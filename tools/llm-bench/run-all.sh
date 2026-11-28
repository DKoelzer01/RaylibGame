#!/bin/bash
# Run every validated task for one model label (#2095). Sequential (run.sh takes the bench lock); 4 h cap.
#   run-all.sh <label> [task ...]
B=~/llm/bench; cd $B; label=${1:?label}; shift
tasks=${*:-"weather-seed noise-seed ship-turn creature-chase trees-clearing fluid-cellkey site-query mining-yield catalog-slots"}
log=$B/run-all-$label.log; : > $log
timeout 14400 bash -c "
for t in $tasks; do
  echo \"== \$t \$(date +%T)\"
  timeout 3300 ./run.sh \$t $label > /tmp/bench-\$t-$label.log 2>&1
  tail -1 /tmp/bench-\$t-$label.log | jq -c '{task, pass, hidden_pass, rounds, steps, wall_s, gen_tok_s, lines_changed, scope_violations}' 2>/dev/null || tail -3 /tmp/bench-\$t-$label.log
done" >> $log 2>&1
echo "all done $(date +%T)" >> $log
