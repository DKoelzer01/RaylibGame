#!/usr/bin/env python3
"""Make a benchmark task from a real commit (#2095).

    mktask.py <repo> <sha> <task-id> [--difficulty easy|medium|hard]

Writes tasks/<task-id>/: task.json, test.patch (the commit's tests: files under tests/ and hunks that add a GoogleTest
block), ref.patch (the rest: the reference implementation), and spec.md (to be written by Claude from the commit).
The task's base is the commit's parent; the model starts from base + test.patch and must make the test pass.
Byte-exact (the repo's files are CRLF).
"""
import json, os, re, subprocess, sys

TEST_HUNK = re.compile(rb"^\+.*(TEST(_F)?\(|#ifdef RAYLIBGAME_UNIT_TESTS|SCENARIO\()")


def git(repo, *args) -> bytes:
    return subprocess.run(["git", "-C", repo, *args], check=True, capture_output=True).stdout


def split(diff: bytes):
    """Unified diff -> (test patch, code patch, code files, test names, scenario names)."""
    test, code, code_files, gtests, scenarios = [], [], [], [], []
    files = re.split(rb"(?m)^(?=diff --git )", diff)
    for f in files:
        if not f.startswith(b"diff --git "):
            continue
        parts = re.split(rb"(?m)^(?=@@ )", f)
        header, hunks = parts[0], parts[1:]
        path = re.search(rb"^diff --git a/(\S+) b/", header).group(1).decode()
        for h in hunks:
            for m in re.finditer(rb"^\+.*TEST(?:_F)?\(\s*(\w+)\s*,\s*(\w+)\s*\)", h, re.M):
                gtests.append(f"{m.group(1).decode()}.{m.group(2).decode()}")
            for m in re.finditer(rb'^\+.*SCENARIO\(\s*"([^"]+)"', h, re.M):
                scenarios.append(m.group(1).decode())
        if path.startswith("tests/"):
            test.append(f)
            continue
        t = [h for h in hunks if any(TEST_HUNK.match(l) for l in h.splitlines())]
        c = [h for h in hunks if h not in t]
        if t:
            test.append(header + b"".join(t))
        if c:
            code.append(header + b"".join(c))
            code_files.append(path)
    return b"".join(test), b"".join(code), code_files, gtests, scenarios


def main():
    repo, sha, tid = sys.argv[1], sys.argv[2], sys.argv[3]
    diff_level = sys.argv[sys.argv.index("--difficulty") + 1] if "--difficulty" in sys.argv else "medium"
    full = git(repo, "rev-parse", sha).decode().strip()
    base = git(repo, "rev-parse", f"{full}^").decode().strip()
    diff = git(repo, "show", "--pretty=", "--binary", full)
    test, code, files, gtests, scenarios = split(diff)
    if not test or not code:
        sys.exit(f"{sha}: needs both a test part and a code part (test {len(test)} B, code {len(code)} B)")
    if gtests:
        kind, test_cmd = "unit", "$TESTS --gtest_brief=1 --gtest_filter=" + ":".join(gtests)
    else:
        kind, test_cmd = "scenario", "$GAME --scenario " + ",".join(scenarios) + " --scenario-no-record --scenario-out /run/scenario-out"
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tasks", tid)
    os.makedirs(d, exist_ok=True)
    open(f"{d}/test.patch", "wb").write(test)
    open(f"{d}/ref.patch", "wb").write(code)
    subject = git(repo, "log", "-1", "--pretty=%s", full).decode().strip()
    stat = git(repo, "show", "--numstat", "--pretty=", full).decode().split("\n")
    added = sum(int(l.split()[0]) + int(l.split()[1]) for l in stat if l and l.split()[0].isdigit() and not l.split()[2].startswith("tests/"))
    task = {"id": tid, "commit": full, "base": base, "subject": subject, "kind": kind, "difficulty": diff_level,
            "files": files, "test_cmd": test_cmd, "ref_lines_changed": added,
            "file_tokens_k": {f: round(len(git(repo, "show", f"{base}:{f}")) / 4000, 1) for f in files}}
    json.dump(task, open(f"{d}/task.json", "w"), indent=2)
    if not os.path.exists(f"{d}/spec.md"):
        open(f"{d}/spec.md", "w").write(f"TODO: spec for {subject}\n")
    print(json.dumps(task))


if __name__ == "__main__":
    main()
