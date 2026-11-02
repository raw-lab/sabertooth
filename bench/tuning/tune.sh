#!/bin/bash
# tune.sh GRID TARGET EXTRA_IDS|- BASEGATE "T:C T:C ..." [extra sabertooth args]
# Environment: SABERTOOTH (binary), DS (make_datasets.py output with q_tune.fa),
# TUNE_OUT (output root), GHI (highest replayed gate, default 46).
set -e
DB=$1; TARGET=$2; EXTRA=$3; G0=$4; GRID=$5; shift 5
HERE=$(cd "$(dirname "$0")" && pwd)
N=${SABERTOOTH:-sabertooth}
DS=${DS:?set DS to the dataset directory}
OUT=${TUNE_OUT:-tune_runs}/$DB
mkdir -p $OUT
for tc in $GRID; do
  T=${tc%:*}; C=${tc#*:}
  tag=t${T}_c${C}
  [ -f $OUT/$tag.gates.json ] && continue
  SABERTOOTH_TRACE_CANDIDATES=$OUT/$tag.cand $N search $DS/q_tune.fa $TARGET --evalue 10 --max-hits 4000 --threads 2 \
     --out $OUT/$tag.m8 --stats-json $OUT/$tag.stats.json --kmer-score $T --max-kmers $C --min-diag-score $G0 "$@" 2>/dev/null
  EX=""; [ "$EXTRA" != "-" ] && EX="--db-extra $EXTRA"
  python3 $HERE/../tune_gate.py --trace $OUT/$tag.cand --hits $OUT/$tag.m8 \
     --labels $DS/scope40.labels.tsv --queries $DS/q_tune.labels.tsv $EX --gates $G0-${GHI:-46} --json $OUT/$tag.gates.json > $OUT/$tag.gates.txt
  rm -f $OUT/$tag.cand
  python3 - $OUT/$tag.stats.json <<'PY'
import json,sys
s=json.load(open(sys.argv[1])); c=s['counters']
print('%s wall %.2f cpu seeds %.2f hits %.2f ev %.2f gap %.2f tb %.2f kmers/q %.0f' % (sys.argv[1].split('/')[-1], s['wall_s']['search'], c['cpu_s_seed_generation'], c['cpu_s_hit_collection_sort'], c['cpu_s_diagonal_evidence'], c['cpu_s_gapped_scoring'], c['cpu_s_traceback'], c['kmers_generated']/c['queries']))
PY
done
