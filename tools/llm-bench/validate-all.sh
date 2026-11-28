#!/bin/bash
# Validate tasks: base + tests must fail, + the reference must pass (#2095). Sequential; 90 min cap.
#   validate-all.sh [task ...]   (default: every task in tasks/)
B=~/llm/bench; cd $B
tasks=${*:-$(ls tasks)}
: > validate.log
for t in $tasks; do
  echo "== $t $(date +%T)" >> validate.log
  timeout 1200 ./run.sh $t validate > /tmp/val-$t.log 2>&1
  { grep -E "^(base|ref) exit" /tmp/val-$t.log | tr "\n" " "; echo; } >> validate.log
done
echo "all done $(date +%T)" >> validate.log
