#!/usr/bin/env python3
"""Summarize the local-model benchmark (#2095): results.jsonl -> markdown tables.

    report.py [results.jsonl] [--label L ...]   (default: every label except validate; the latest run per task+label)
"""
import json, statistics as st, sys

path = next((a for a in sys.argv[1:] if not a.startswith("--") and a.endswith(".jsonl")), "results.jsonl")
want = [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--label"]
rows = [json.loads(l) for l in open(path) if l.strip()]
val = {r["task"]: r for r in rows if r["label"] == "validate"}
latest = {}
for r in rows:
    if r["label"] != "validate" and (not want or r["label"] in want):
        latest[(r["task"], r["label"])] = r
labels = sorted({l for _, l in latest})
tasks = sorted({t for t, _ in latest}, key=lambda t: ({"easy": 0, "medium": 1, "hard": 2}.get(latest[next(k for k in latest if k[0] == t)]["difficulty"], 3), t))


import os, re
BENCH = os.path.dirname(os.path.abspath(path))


def new_failures(r):
    """Unit tests failing after the model's change that don't fail at the reference (tasks/<id>/hidden_baseline.txt)."""
    try:
        log = open(f"{BENCH}/runs/{r['run']}/hidden.log", errors="replace").read()
    except OSError:
        return None
    failed = set(re.findall(r"^\[  FAILED  \] (\w+\.\w+)", log, re.M))
    try:
        base = set(open(f"{BENCH}/tasks/{r['task']}/hidden_baseline.txt").read().split())
    except OSError:
        return sorted(failed) if failed else []  # no baseline yet: every failure counts
    return sorted(failed - base)


for r in latest.values():
    nf = new_failures(r)
    r["hidden_new"] = nf
    r["hidden_pass"] = nf is not None and not nf


def fmt(x, nd=0):
    return "-" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


print("## Tasks\n")
print("| task | kind | difficulty | reference lines | validated (base fails, ref passes) |")
print("|---|---|---|---|---|")
for t in tasks:
    r = next(latest[k] for k in latest if k[0] == t)
    v = val.get(t)
    ok = "yes" if v and v["base_exit"] not in ("0", "") and v["ref_exit"] == "0" else ("NO" if v else "not run")
    print(f"| {t} | {r['kind']} | {r['difficulty']} | {r['ref_lines']} | {ok} |")

for label in labels:
    try:
        reviews = json.load(open(f"{BENCH}/reviews-{label}.json"))
    except OSError:
        reviews = {}
    rs = [latest[(t, label)] for t in tasks if (t, label) in latest]
    passed = [r for r in rs if r["pass"]]
    print(f"\n## {label}: {len(passed)}/{len(rs)} pass\n")
    print("| task | pass | unit suite (vs reference) | rounds | steps | wall | build+test | model gen | gen tok/s | lines (ref) | scope | CRLF / debug / warnings | edit errors |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rs:
        print(f"| {r['task']} | {'PASS' if r['pass'] else 'fail'} | {'ok' if r['hidden_pass'] else 'new failures: ' + ', '.join(r['hidden_new'] or ['?'])} | {fmt(r['rounds'])} | {fmt(r['steps'])} "
              f"| {r['wall_s'] // 60}m{r['wall_s'] % 60:02d}s | {r['build_test_s']}s | {fmt(r['gen_s'], 0)}s | {fmt(r['gen_tok_s'], 1)} "
              f"| {r['lines_changed']} ({r['ref_lines']}) | {r['scope_violations'] or 'ok'} | {r['crlf_bad_lines']} / {r['debug_lines']} / {r['new_warnings']} "
              f"| {fmt(r['edit_errors'])} |")
    if rs:
        pr = [r["rounds"] for r in passed if r["rounds"]]
        print(f"\n- pass rate {100 * len(passed) / len(rs):.0f}%; within 2 test runs: {sum(1 for r in passed if (r['rounds'] or 9) <= 2)}/{len(rs)}"
              f"; median rounds when passing: {st.median(pr) if pr else '-'}")
        print(f"- median wall {st.median(r['wall_s'] for r in rs) / 60:.1f} min; model share of wall "
              f"{100 * sum((r['prompt_s'] or 0) + (r['gen_s'] or 0) for r in rs) / max(1, sum(r['wall_s'] for r in rs)):.0f}%; "
              f"tokens in/out {sum(r['prompt_tokens'] or 0 for r in rs):.0f} / {sum(r['gen_tokens'] or 0 for r in rs):.0f}")
        bad = [r["task"] for r in rs if r["scope_violations"] or r["crlf_bad_lines"] or r["debug_lines"]]
        print(f"- hygiene problems (scope, line endings, debug output): {', '.join(bad) or 'none'}")
        if reviews:
            rev = [reviews.get(r["task"], {}).get("verdict") for r in passed]
            print(f"- Claude review of the passing diffs: {rev.count('correct')} correct, {rev.count('minor')} minor, "
                  f"{rev.count('defect')} defect (bugs the tests miss), {rev.count(None)} not reviewed")
            for r in rs:
                v = reviews.get(r["task"])
                if v and v["verdict"] != "correct":
                    print(f"  - {r['task']}: **{v['verdict']}**: {v['note']}")
