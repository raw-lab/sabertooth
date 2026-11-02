#!/bin/bash
# A/B: hit-loop posting prefetch on/off, 68M spaced index, default-like config, alternating runs
Q=/home/claude/bench_data/ds_scale/q_tune.fa
IDX=/home/claude/exp_dbs/st_scale_sp9.stidx
for i in 1 2 3; do
  for v in pf nopf; do
    /home/claude/bin/st_$v search $Q $IDX --evalue 10 --max-hits 4000 --threads 2 --spaced 111001101 --kmer-score 22 --max-kmers 512 --min-diag-score 41 \
      --out /dev/null --stats-json /tmp/claude-0/-home-claude/6d52089c-4cc5-5fb7-b131-8bb4426efce6/scratchpad/ab_$v.$i.json 2>/dev/null
    python3 -c "
import json,sys
s=json.load(open(sys.argv[1])); c=s['counters']
print('%-5s run %s wall %.2f seeds %.2f hits %.2f ev %.2f gap %.2f tb %.2f' % (sys.argv[2], sys.argv[3], s['wall_s']['search'], c['cpu_s_seed_generation'], c['cpu_s_hit_collection_sort'], c['cpu_s_diagonal_evidence'], c['cpu_s_gapped_scoring'], c['cpu_s_traceback']))
" /tmp/claude-0/-home-claude/6d52089c-4cc5-5fb7-b131-8bb4426efce6/scratchpad/ab_$v.$i.json $v $i
  done
done
