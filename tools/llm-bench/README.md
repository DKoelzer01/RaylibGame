# Local-model benchmark (#2095)

Measures how well a local model (llama.cpp on deweyserver) implements small, Claude-specified chunks of this repo, and
how fast: the trial in the local-model plan (#1834) before sessions delegate work to it. Results and the reviews so
far are on #2095.

## How it works
- **Tasks** come from real commits (`mktask.py <repo> <sha> <id> [--difficulty easy|medium|hard]`): the task's base is
  the commit's parent; the commit's tests (files under `tests/`, hunks that add a GoogleTest or a `SCENARIO`) become
  `test.patch`, the rest `ref.patch`. Claude writes `spec.md` the way a chunk spec is written: the problem, the files and
  functions, the signatures to add, the behaviour; not the code.
- **Validation** (`run.sh <task> validate`): base + tests must fail, + the reference must pass.
- **A run** (`run.sh <task> <label>`): a fresh checkout at the base with the tests applied; the agent implements the spec;
  then the task's test (`final`) and the whole unit suite (`hidden`). One JSON row per run goes to `results.jsonl`.
- **Agents:** `AGENT=mini` (default, `bin/mini_agent.py`): grep, read a line range, exact edit, build+test, done. Our
  files are 27-70K tokens, more than a 32K slot holds, so it never puts a whole file in context (`CTX_CHARS` = its
  conversation budget, ~4 characters a token). `AGENT=aider` (`bin/aider_agent.sh`): Aider 0.86.2 with its own test loop.
- **Builds** run in the `llm-bench:1` image (`image/Dockerfile`: the CI runner image + Aider) at low priority next to CI
  (4 cpus, cpu-shares 512), on the 2080 SUPER for scenarios, with the build sandbox's ccache. One job at a time
  (`slot1.lock`). Scenario verdicts come from `results.json` (a single-name run exits 0 on FAIL: #2137).
- **Measures:** pass, unit-suite failures that the reference doesn't have (`hidden_baseline.txt`, from `baseline.sh`),
  test runs, steps, wall/build/model time, prompt and generated tokens and tok/s (llama.cpp's /metrics), lines changed
  vs the reference, files outside the spec, line endings, debug output, new warnings in touched files.
- **Review:** a passing test is not a correct change (every configuration has produced passing-but-wrong diffs). Claude
  reviews each passing diff against the spec and the reference: `reviews-<label>.json` (correct / minor / defect).
  `mutant.sh <task> <diff>` checks that a task's test catches a known-bad diff (the `mutant-*.diff` files).

## Use (on deweyserver, as claude)
```
rsync -a tools/llm-bench/ deweyserver:llm/bench/          # deploy (scripts expect ~/llm/bench)
docker build -t llm-bench:1 ~/llm/bench/image              # once
cd ~/llm/bench && ./validate-all.sh                        # check the tasks
CTX_CHARS=80000 ./run-all.sh <label> [task ...]            # our agent, ~20K-token budget (the best setup so far)
AGENT=aider ./run-all.sh <label> [task ...]                # Aider
./baseline.sh <task ...>                                   # unit tests failing with the reference
python3 report.py results.jsonl --label <label>            # tables (+ reviews-<label>.json)
```
Environment assumptions (deweyserver): the model at `localhost:8095` (`~/llm/docker-compose.yml`), the repo clone and
ccache of the build sandbox (`~/rg-sandbox`), the runners' GPU lock (`~/runner-stack/cache-1`), the 2080 SUPER's UUID
in `run.sh`/`mutant.sh`. `repo/`, `work1/`, `out1/`, `runs/` and `results.jsonl` are created there and not kept here.

## Tasks
`*-v2` tasks have stronger tests that catch the defects the first runs got past (an unseen target chased, a centre with
x == 0 ignored). `mining-yield` was dropped (its scenario fails with the reference fix at that commit).
