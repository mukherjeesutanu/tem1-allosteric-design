#!/bin/bash
# wait for the allosteric screen to finish, then dock the controls
cd ~/tem1-allosteric-design
while pgrep -f "run_screen.py --workers" > /dev/null 2>&1; do sleep 30; done
echo "allosteric done: $(wc -l < data/generated/docked_allosteric.jsonl) results" >> outputs/chain.log
PYTHONPATH=. ~/anaconda3/envs/vina/bin/python run_screen.py \
  --workers 22 --exhaustiveness 8 \
  --sites orthosteric,surface --top-from-allosteric 400 >> outputs/controls.log 2>&1
echo "controls done" >> outputs/chain.log
