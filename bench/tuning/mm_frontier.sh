#!/bin/bash
# mm_frontier.sh NAME MMDB EXTRA_IDS|- "s1 s2 ..."  (tune queries)
# Environment: MMSEQS (binary), DS (dataset directory), TUNE_OUT (output root).
NAME=$1; DB=$2; EXTRA=$3; SLIST=$4
HERE=$(cd "$(dirname "$0")" && pwd)
MM=${MMSEQS:-mmseqs}
DS=${DS:?set DS to the dataset directory}
OUT=${TUNE_OUT:-tune_runs}/mm_$NAME
mkdir -p $OUT
for s in $SLIST; do
  tag=s$s
  [ -f $OUT/$tag.eval.json ] && continue
  rm -rf $OUT/tmp
  /usr/bin/time -f "%e %U %S %M" -o $OUT/$tag.time $MM easy-search $DS/q_tune.fa $DB $OUT/$tag.m8 $OUT/tmp -e 10 --max-seqs 4000 --threads 2 -s $s -v 1
  rm -rf $OUT/tmp
  EX=""; [ "$EXTRA" != "-" ] && EX="--db-extra $EXTRA"
  python3 $HERE/../eval_scop.py --labels $DS/scope40.labels.tsv --queries $DS/q_tune.labels.tsv $EX \
     --hits $OUT/$tag.m8 --name $tag --json $OUT/$tag.eval.json --recall-levels 0.05,0.06,0.07,0.08 > /dev/null
  python3 - $OUT/$tag <<'PY'
import json,sys
p=sys.argv[1]; e=json.load(open(p+'.eval.json')); t=open(p+'.time').read().split()
th={x['evalue']:x for x in e['thresholds']}
print('%s wall %s user %s sys %s rssKB %s AUC1 %.4f TP/FP@1e-3 %d/%d FP@rec %s' % (p.split('/')[-1], t[0], t[1], t[2], t[3], e['mean_auc1'], th[1e-3]['TP'], th[1e-3]['FP'], [(m['recall'], m['FP']) for m in e['matched_recall']]))
PY
  rm -f $OUT/$tag.m8
done
