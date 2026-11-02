#!/usr/bin/env python3
"""Evaluate every prefilter gate >= the run's own gate from ONE traced search.

A stricter ungapped gate keeps a subset of the prefilter candidates and never
changes their gapped scores or E-values, so the hit list a stricter gate would
produce is exactly the traced run's hit list restricted to candidates whose
prefilter score passes it. This script replays that for a range of gates and
reports, per gate: candidates per query, DP cells per query (gapped-alignment
work), mean AUC1, TP/FP at E <= 1e-3 and FP at matched recall.

  SABERTOOTH_TRACE_CANDIDATES=run.cand sabertooth search Q DB --min-diag-score 26 ... --out run.m8
  tune_gate.py --trace run.cand --hits run.m8 --labels scope40.labels.tsv \
               --queries q.labels.tsv [--db-extra extra.ids] --gates 26-44 [--json out.json]

Trace format (one line per candidate): query, target, prefilter score, cells.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eval_scop as ev  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trace', required=True)
    ap.add_argument('--hits', required=True)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--queries', required=True)
    ap.add_argument('--db-extra', default=None)
    ap.add_argument('--gates', default='26-44')
    ap.add_argument('--recall-levels', default='0.05,0.06,0.07,0.08')
    ap.add_argument('--json', default=None)
    a = ap.parse_args()
    g0, g1 = (int(x) for x in a.gates.split('-'))
    gates = list(range(g0, g1 + 1))

    lab, _ = ev.read_labels(a.labels)
    ev.WINSHUF_LAB.update(lab)
    queries = ev.read_query_labels(a.queries, lab)
    sf_count = Counter(ev.sf(s) for s in lab.values())
    q_tp = {}
    tp_total = 0
    for q, (srcs, _qlen, src_ids) in queries.items():
        n = sum(sf_count[x] for x in {ev.sf(s) for s in srcs}) - sum(1 for d in src_ids if d in lab)
        q_tp[q] = n
        tp_total += n
    hits = ev.read_hits(a.hits)

    # prefilter score of every (query, target) candidate; per-gate totals
    score = {}
    n_cand = Counter()
    n_cells = Counter()
    with open(a.trace) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                continue
            s = int(p[2])
            c = int(p[3])
            score[(p[0], p[1])] = max(s, score.get((p[0], p[1]), -10**9))
            # cumulative: a candidate counts for every gate it passes
            n_cand[s] += 1
            n_cells[s] += c
    smax = max(n_cand) if n_cand else 0
    cum_cand = {}
    cum_cells = {}
    acc_c = acc_x = 0
    for s in range(max(smax, g1), g0 - 1, -1):
        acc_c += n_cand.get(s, 0)
        acc_x += n_cells.get(s, 0)
        cum_cand[s] = acc_c
        cum_cells[s] = acc_x
    # candidates scoring below g0 (if the run used a lower gate) never matter
    nq = len(queries)
    levels = [float(x) for x in a.recall_levels.split(',')]
    out = []
    for g in gates:
        auc = []
        lh = []
        for q, (srcs, _qlen, src_ids) in queries.items():
            seen = 0
            first_fp = False
            for e, _b, t in hits.get(q, []):
                if t in src_ids:
                    continue
                if score.get((q, t), -10**9) < g:
                    continue
                c = ev.classify(srcs, lab.get(t), t)
                lh.append((e, c == 'TP', c == 'FP'))
                if c == 'TP' and not first_fp:
                    seen += 1
                if c == 'FP':
                    first_fp = True
            if q_tp[q] > 0:
                auc.append(seen / q_tp[q])
        lh.sort()
        tp3 = sum(1 for e, t, f in lh if e <= 1e-3 and t)
        fp3 = sum(1 for e, t, f in lh if e <= 1e-3 and f)
        fp_at = {}
        for lvl in levels:
            need = lvl * tp_total
            tp = fp = 0
            fp_at[lvl] = None
            for e, t, f in lh:
                tp += t
                fp += f
                if tp >= need:
                    fp_at[lvl] = fp
                    break
        row = {'gate': g, 'cand_per_q': cum_cand.get(g, 0) / nq, 'cells_per_q': cum_cells.get(g, 0) / nq,
               'auc1': sum(auc) / len(auc) if auc else 0.0, 'tp_1e3': tp3, 'fp_1e3': fp3,
               'tp_10': sum(1 for e, t, f in lh if t), 'fp_at_recall': {str(k): v for k, v in fp_at.items()}}
        out.append(row)
        print('gate %3d  cand/q %8.1f  Mcells/q %7.2f  AUC1 %.4f  TP/FP@1e-3 %5d/%-3d  TP@10 %5d  FP@rec %s' % (
            g, row['cand_per_q'], row['cells_per_q'] / 1e6, row['auc1'], tp3, fp3, row['tp_10'],
            ' '.join('%.2f:%s' % (k, v) for k, v in fp_at.items())), flush=True)
    if a.json:
        with open(a.json, 'w') as f:
            json.dump({'trace': a.trace, 'hits': a.hits, 'queries': nq, 'rows': out}, f, indent=1)


if __name__ == '__main__':
    main()
