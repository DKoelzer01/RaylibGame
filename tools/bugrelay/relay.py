#!/usr/bin/env python3
"""Bug report relay (#169): the game's "Report a bug" window posts here; the relay stores the attachments in the Claude
media store and opens a GitHub issue as the user (with their token, read from a file the user places; agents never
handle it). Tailnet-only for now (bind it behind `tailscale serve`, or to the LAN/tailnet address); public later.

Standard library only. Configuration (environment):
  BUGRELAY_PORT          8765
  BUGRELAY_BIND          0.0.0.0
  BUGRELAY_TOKEN_FILE    /run/secrets/github_token   (a fine-grained token: this repo, Issues read/write)
  BUGRELAY_REPO          DKoelzer01/RaylibGame
  BUGRELAY_MEDIA         /srv/claude-media/bug-reports        (where attachments are written)
  BUGRELAY_SHARE         \\\\DEWEYSERVER\\files\\claude-media\\bug-reports   (how the issue links them)
  BUGRELAY_MAX_MB        64
  BUGRELAY_DRY_RUN       1: store the report and print the issue instead of posting it (tests)

The game sends multipart/form-data: title, body, labels (comma separated), severity, meta (text), and the optional files
screenshot, log, capture_sheet, capture_gif, save. The answer is JSON: {"url": ..., "number": N} or {"error": ...}.
"""
import datetime
import email.parser
import email.policy
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

PORT = int(os.environ.get("BUGRELAY_PORT", "8765"))
BIND = os.environ.get("BUGRELAY_BIND", "0.0.0.0")
TOKEN_FILE = os.environ.get("BUGRELAY_TOKEN_FILE", "/run/secrets/github_token")
REPO = os.environ.get("BUGRELAY_REPO", "DKoelzer01/RaylibGame")
MEDIA = os.environ.get("BUGRELAY_MEDIA", "/srv/claude-media/bug-reports")
SHARE = os.environ.get("BUGRELAY_SHARE", r"\\DEWEYSERVER\files\claude-media\bug-reports")
MAX_BYTES = int(float(os.environ.get("BUGRELAY_MAX_MB", "64")) * 1024 * 1024)
DRY_RUN = os.environ.get("BUGRELAY_DRY_RUN", "") == "1"
FILES = {"screenshot": "screenshot.png", "log": "log.txt", "capture_sheet": "capture_sheet.png", "capture_gif": "capture.gif", "save": "world.sav"}

_last_by_client = {}
_lock = threading.Lock()
_labels_cache = {"at": 0.0, "names": set()}


def token():
    with open(TOKEN_FILE, "r", encoding="utf-8") as f:
        return f.read().strip()


def github(method, path, payload=None):
    req = urllib.request.Request("https://api.github.com" + path, method=method,
                                 data=json.dumps(payload).encode() if payload is not None else None,
                                 headers={"Authorization": "Bearer " + token(), "Accept": "application/vnd.github+json",
                                          "X-GitHub-Api-Version": "2022-11-28", "User-Agent": "raylibgame-bugrelay"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


def repo_labels():
    """The repo's existing labels (cached 10 min): only those are applied; other tags go in the body."""
    if time.time() - _labels_cache["at"] > 600:
        names, page = set(), 1
        while True:
            batch = github("GET", f"/repos/{REPO}/labels?per_page=100&page={page}")
            names.update(l["name"] for l in batch)
            if len(batch) < 100:
                break
            page += 1
        _labels_cache.update(at=time.time(), names=names)
    return _labels_cache["names"]


def parse_multipart(content_type, raw):
    msg = email.parser.BytesParser(policy=email.policy.HTTP).parsebytes(
        b"Content-Type: " + content_type.encode() + b"\r\nMIME-Version: 1.0\r\n\r\n" + raw)
    fields, files = {}, {}
    for part in msg.iter_parts():
        name = part.get_param("name", header="content-disposition")
        if not name:
            continue
        data = part.get_payload(decode=True) or b""
        if part.get_filename() is not None:
            files[name] = data
        else:
            fields[name] = data.decode("utf-8", "replace")
    return fields, files


def slug(s):
    s = re.sub(r"[^a-zA-Z0-9]+", "-", s.lower()).strip("-")
    return s[:40] or "report"


class Handler(BaseHTTPRequestHandler):
    def reply(self, code, obj):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/health":
            return self.reply(200, {"ok": True, "repo": REPO, "token": os.path.exists(TOKEN_FILE), "dry_run": DRY_RUN})
        self.reply(404, {"error": "POST /report"})

    def do_POST(self):
        if self.path != "/report":
            return self.reply(404, {"error": "POST /report"})
        client = self.client_address[0]
        with _lock:  # one report per client every 20 s
            if time.time() - _last_by_client.get(client, 0) < 20:
                return self.reply(429, {"error": "slow down: one report every 20 s"})
            _last_by_client[client] = time.time()
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0 or length > MAX_BYTES:
            return self.reply(413, {"error": f"report size {length} bytes (limit {MAX_BYTES})"})
        ctype = self.headers.get("Content-Type", "")
        if not ctype.startswith("multipart/form-data"):
            return self.reply(400, {"error": "multipart/form-data expected"})
        fields, files = parse_multipart(ctype, self.rfile.read(length))
        title = fields.get("title", "").strip()[:200]
        if not title:
            return self.reply(400, {"error": "a title is needed"})
        # The attachments: <media>/<date>/<HHMMSS>-<slug>/
        now = datetime.datetime.now()
        rel = f"{now:%Y-%m-%d}/{now:%H%M%S}-{slug(title)}"
        folder = os.path.join(MEDIA, *rel.split("/"))
        os.makedirs(folder, exist_ok=True)
        stored = []
        for key, fname in FILES.items():
            if key in files and files[key]:
                with open(os.path.join(folder, fname), "wb") as f:
                    f.write(files[key])
                stored.append((key, fname, len(files[key])))
        for key in ("body", "meta"):
            with open(os.path.join(folder, key + ".txt"), "w", encoding="utf-8") as f:
                f.write(fields.get(key, ""))
        # The issue.
        tags = [t.strip() for t in fields.get("labels", "").split(",") if t.strip()]
        severity = fields.get("severity", "minor").strip()
        try:
            known = repo_labels() if not DRY_RUN else set()
        except Exception:
            known = set()
        wanted = ["bug report", "bug"] + tags
        labels = [l for l in wanted if l in known]
        other = [t for t in tags if t not in known]
        share = SHARE.rstrip("\\") + "\\" + rel.replace("/", "\\")
        lines = [fields.get("body", "").strip() or "_(no description)_", "",
                 f"**Severity:** {severity}" + (f"  \n**Tags:** {', '.join(other)}" if other else ""), "",
                 "**Attachments**" if stored else "_No attachments._"]
        for key, fname, size in stored:
            lines.append(f"- {key}: `{share}\\{fname}` ({size // 1024} KB)")
        if stored:
            lines.append(f"- In a Claude session: `media_list bug-reports/{rel}` / `media_get bug-reports/{rel}/<file>`")
        lines += ["", "<details><summary>Game info</summary>", "", "```", fields.get("meta", "").strip(), "```", "</details>", "",
                  "_Sent from the in-game bug report (#169); triage pending._"]
        payload = {"title": title, "body": "\n".join(lines), "labels": labels}
        if DRY_RUN:
            print(json.dumps(payload, indent=1))
            return self.reply(200, {"url": f"dry-run:{share}", "number": 0})
        try:
            issue = github("POST", f"/repos/{REPO}/issues", payload)
        except FileNotFoundError:
            return self.reply(503, {"error": f"no token at {TOKEN_FILE}; the report is kept at {share}"})
        except urllib.error.HTTPError as e:
            return self.reply(502, {"error": f"GitHub said {e.code}: {e.read()[:200].decode('utf-8', 'replace')}; kept at {share}"})
        print(f"{now:%H:%M:%S} {client} -> #{issue['number']} {title}", flush=True)
        self.reply(200, {"url": issue["html_url"], "number": issue["number"]})

    def log_message(self, fmt, *args):  # quieter: the issue line above is the log
        pass


if __name__ == "__main__":
    print(f"bug relay on {BIND}:{PORT} for {REPO}, media {MEDIA}{' (dry run)' if DRY_RUN else ''}", flush=True)
    ThreadingHTTPServer((BIND, PORT), Handler).serve_forever()
