#!/usr/bin/env python3
"""A small coding agent for the local-model benchmark (#2095), stdlib only.

It works like Claude Code on big files: search, read a range of lines, edit by exact replacement, build and test. The
game's source files are 27-70K tokens, more than a 32K slot holds, so whole-file agents (Aider) can't do most chunks.

    mini_agent.py --task /task --run /run --api http://host:8095/v1 [--model NAME] [--max-tests 4] [--max-steps 60]

Reads /task/task.json and /task/spec.md, works in the current directory (the repo), runs /bench/testcmd.sh for the
`test` tool, writes /run/agent.json (stats) and /run/transcript.jsonl.
"""
import argparse, json, os, re, subprocess, time, urllib.request

SYSTEM = """You are a careful C++ engineer implementing one small, fully specified change in the RaylibGame repository
(a C++17 raylib planet game). Claude designed the change; you implement exactly the spec, nothing more.

Tools: grep (search), read (a range of lines; files are big, read only what you need), edit (replace an exact unique
snippet), test (builds and runs the task's test; a limited number of runs), done (finish with a one-line summary).
Workflow: find the code the spec names with grep, read around it, make the smallest edits that satisfy the spec, call
test, fix what fails, call done when the test passes (or when you can't make progress).
Rules:
- Edit only the files the task allows. Keep the existing style: comments say why, numbers carry units.
- No debug prints, no commented-out code, no unrelated changes. Don't change the test.
- Planet-scale positions stay double (DVec3); GL calls only on the main thread; the factory sim stays deterministic.
- `edit.old` must match the file exactly (indentation included) and only once: read the lines first."""

TOOLS = [
    {"type": "function", "function": {"name": "grep", "description": "Search the repo with an extended regex. Returns up to 60 matching lines as path:line:text.",
        "parameters": {"type": "object", "properties": {"pattern": {"type": "string"}, "path": {"type": "string", "description": "file or folder, default src"}}, "required": ["pattern"]}}},
    {"type": "function", "function": {"name": "read", "description": "Read lines start..end (1-based, at most 250 lines) of a file, numbered.",
        "parameters": {"type": "object", "properties": {"path": {"type": "string"}, "start": {"type": "integer"}, "end": {"type": "integer"}}, "required": ["path", "start", "end"]}}},
    {"type": "function", "function": {"name": "edit", "description": "Replace the exact text `old` (must occur once) with `new` in a file the task allows.",
        "parameters": {"type": "object", "properties": {"path": {"type": "string"}, "old": {"type": "string"}, "new": {"type": "string"}}, "required": ["path", "old", "new"]}}},
    {"type": "function", "function": {"name": "test", "description": "Build the game and run the task's test. Returns build errors or the test output. Limited runs.",
        "parameters": {"type": "object", "properties": {}}}},
    {"type": "function", "function": {"name": "done", "description": "Finish. Give a one-line summary of the change.",
        "parameters": {"type": "object", "properties": {"summary": {"type": "string"}}, "required": ["summary"]}}},
]


class Agent:
    def __init__(self, a):
        self.a = a
        self.task = json.load(open(f"{a.task}/task.json"))
        self.spec = open(f"{a.task}/spec.md").read()
        self.stats = {"steps": 0, "tests": 0, "edits": 0, "edit_errors": 0, "out_of_scope": 0, "tool_errors": 0,
                      "prompt_tokens": 0, "completion_tokens": 0, "model_s": 0.0, "test_s": 0.0, "finished": False,
                      "summary": "", "stop": ""}
        self.log = open(f"{a.run}/transcript.jsonl", "a")

    # ---------------------------------------------------------------- tools
    def t_grep(self, pattern, path="src"):
        p = subprocess.run(["grep", "-rnE", "--include=*.cpp", "--include=*.h", "--include=*.txt", "--include=*.glsl",
                            pattern, path or "src"], capture_output=True, text=True, errors="replace")
        lines = [l[:220] for l in p.stdout.splitlines()]
        return "\n".join(lines[:60]) + (f"\n[{len(lines)} matches, first 60 shown]" if len(lines) > 60 else "") or "no matches"

    def t_read(self, path, start, end):
        start, end = max(1, int(start)), int(end)
        end = min(end, start + 249)
        with open(path, errors="replace", newline="") as f:
            lines = f.read().split("\n")
        part = lines[start - 1:end]
        return "\n".join(f"{start + i}\t{l.rstrip(chr(13))}" for i, l in enumerate(part)) + f"\n[{len(lines)} lines in file]"

    def t_edit(self, path, old, new):
        if os.path.normpath(path) not in [os.path.normpath(f) for f in self.task["files"]]:
            self.stats["out_of_scope"] += 1
            return f"refused: {path} is not one of the files this task allows ({', '.join(self.task['files'])})"
        data = open(path, "rb").read().decode("utf-8", "replace")
        crlf = "\r\n" in data
        if crlf:
            old = old.replace("\r\n", "\n").replace("\n", "\r\n")
            new = new.replace("\r\n", "\n").replace("\n", "\r\n")
        n = data.count(old)
        if n != 1:
            self.stats["edit_errors"] += 1
            return f"edit failed: `old` occurs {n} times (it must occur exactly once). Read the lines and copy them exactly."
        open(path, "wb").write(data.replace(old, new, 1).encode("utf-8"))
        self.stats["edits"] += 1
        return "ok"

    def t_test(self):
        if self.stats["tests"] >= self.a.max_tests:
            return f"no test runs left ({self.a.max_tests} used): call done."
        self.stats["tests"] += 1
        t0 = time.time()
        p = subprocess.run(["bash", "/bench/testcmd.sh", f"test-{self.stats['tests']}"], capture_output=True, text=True, errors="replace")
        self.stats["test_s"] += time.time() - t0
        out = (p.stdout + p.stderr)[-3500:]
        return f"{'PASS' if p.returncode == 0 else 'FAIL'} (run {self.stats['tests']} of {self.a.max_tests})\n{out}"

    # ---------------------------------------------------------------- model
    def chat(self, messages):
        body = {"model": self.a.model, "messages": messages, "tools": TOOLS, "temperature": 0.2, "max_tokens": 2048}
        if self.a.no_think:
            body["chat_template_kwargs"] = {"enable_thinking": False}
        req = urllib.request.Request(f"{self.a.api}/chat/completions", json.dumps(body).encode(),
                                     {"Content-Type": "application/json", "Authorization": "Bearer x"})
        t0 = time.time()
        with urllib.request.urlopen(req, timeout=600) as r:
            j = json.load(r)
        self.stats["model_s"] += time.time() - t0
        u = j.get("usage", {})
        self.stats["prompt_tokens"] += u.get("prompt_tokens", 0)
        self.stats["completion_tokens"] += u.get("completion_tokens", 0)
        return j["choices"][0]["message"]

    def trim(self, messages, budget_chars):
        """Keep the conversation under the slot's context: shorten the oldest tool results first."""
        size = lambda: sum(len(json.dumps(m)) for m in messages)
        i = 2
        while size() > budget_chars and i < len(messages) - 4:
            m = messages[i]
            if m["role"] == "tool" and len(m["content"]) > 200:
                m["content"] = m["content"][:120] + "\n[older output trimmed]"
            i += 1

    def run(self):
        files = ", ".join(self.task["files"])
        user = (f"{self.spec}\n\nFiles you may edit: {files}\nThe test: {self.task['test_cmd']} "
                f"(already in the tree; the `test` tool builds and runs it).")
        messages = [{"role": "system", "content": SYSTEM}, {"role": "user", "content": user}]
        tools = {"grep": self.t_grep, "read": self.t_read, "edit": self.t_edit, "test": self.t_test}
        while self.stats["steps"] < self.a.max_steps:
            self.stats["steps"] += 1
            self.trim(messages, self.a.ctx_chars)
            try:
                msg = self.chat(messages)
            except Exception as e:
                self.stats["stop"] = f"model error: {e}"
                break
            messages.append({k: v for k, v in msg.items() if k in ("role", "content", "tool_calls") and v})
            self.log.write(json.dumps({"step": self.stats["steps"], "assistant": msg}) + "\n")
            calls = msg.get("tool_calls") or []
            if not calls:
                messages.append({"role": "user", "content": "Use the tools: continue, or call done."})
                continue
            for c in calls:
                name = c["function"]["name"]
                try:
                    args = json.loads(c["function"].get("arguments") or "{}")
                except json.JSONDecodeError:
                    args = None
                if name == "done":
                    self.stats["finished"], self.stats["summary"] = True, (args or {}).get("summary", "")
                    self.stats["stop"] = "done"
                    break
                try:
                    out = tools[name](**args) if args is not None else "bad arguments: send valid JSON"
                except Exception as e:
                    self.stats["tool_errors"] += 1
                    out = f"tool error: {e}"
                self.log.write(json.dumps({"step": self.stats["steps"], "tool": name, "args": args, "out": out[:4000]}) + "\n")
                messages.append({"role": "tool", "tool_call_id": c.get("id", name), "content": out})
            if self.stats["finished"]:
                break
        else:
            self.stats["stop"] = "max steps"
        json.dump(self.stats, open(f"{self.a.run}/agent.json", "w"), indent=2)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--task", default="/task"); p.add_argument("--run", default="/run")
    p.add_argument("--api", required=True); p.add_argument("--model", default="local")
    p.add_argument("--max-tests", type=int, default=4); p.add_argument("--max-steps", type=int, default=60)
    p.add_argument("--ctx-chars", type=int, default=80000, help="conversation budget in characters (~4 per token)")
    p.add_argument("--no-think", action="store_true", default=True)
    Agent(p.parse_args()).run()
