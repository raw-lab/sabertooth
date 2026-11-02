#!/usr/bin/env python3
"""Cost-model frontier of (kmer score, cap, gate) from tuning runs.

frontier.py DIR c_kmer c_posting c_diag c_cell c_hit
CPU ns per unit; prints, per configuration, modeled CPU s per 1000 queries and
AUC1, then the Pareto frontier.
"""
import json, glob, sys, os

d = sys.argv[1]
ck, cp, cd, cc, ch = (float(x) for x in sys.argv[2:7])
pts = []
for f in sorted(glob.glob(os.path.join(d, '*.gates.json'))):
    tag = os.path.basename(f).replace('.gates.json', '')
    T = int(tag.split('_')[0][1:]); C = int(tag.split('_')[1][1:])
    st = json.load(open(os.path.join(d, tag + '.stats.json')))['counters']
    nq = st['queries']
    km, po, dg, hits = st['kmers_generated'] / nq, st['postings_visited'] / nq, st['diagonals_scored'] / nq, st['hits'] / nq
    for r in json.load(open(f))['rows']:
        g = r['gate']
        cost = (km * ck + po * cp + dg * cd + r['cells_per_q'] * cc + hits * ch) * 1e-9 * 1000  # CPU s per 1000 q
        pts.append({'T': T, 'C': C, 'G': g, 'cpu': cost, 'auc1': r['auc1'], 'km': km, 'po': po, 'dg': dg,
                    'cand': r['cand_per_q'], 'mcell': r['cells_per_q'] / 1e6, 'fp08': r['fp_at_recall'].get('0.08')})
pts.sort(key=lambda p: p['cpu'])
front = []
best = -1
for p in pts:
    if p['auc1'] > best + 1e-9:
        front.append(p)
        best = p['auc1']
print('Pareto frontier (modeled CPU s / 1000 queries vs AUC1):')
for p in front:
    print('  T%2d C%5d G%2d  cpu %7.1f  AUC1 %.4f  kmers/q %7.0f post/q %8.0f diag/q %8.0f cand/q %7.1f Mcell/q %6.1f FP@.08 %s' % (
        p['T'], p['C'], p['G'], p['cpu'], p['auc1'], p['km'], p['po'], p['dg'], p['cand'], p['mcell'], p['fp08']))
