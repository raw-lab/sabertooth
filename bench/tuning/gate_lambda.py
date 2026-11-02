#!/usr/bin/env python3
"""Optimal gate per DB size under  max AUC1 - lam*ln(cost)  for given (T,C).

gate_lambda.py LAMBDAS T:C[,T:C...] [ROOT]   (ROOT: directory holding the grids,
default ../results/v0.4/tuning)
"""
import json, os, math, sys
consts = {  # c_k, c_p, c_d, c_cell, c_tb (ns), uncontended estimates per size
    'scope40_sp9': (70, 50, 120, 0.131, 92000),
    's17_sp9': (70, 54, 126, 0.129, 91000),
    'scale_sp9': (70, 54, 114, 0.131, 90000),
    's265_sp9': (90, 54, 161, 0.137, 100000),
}
size = {'scope40_sp9': 2.63e6, 's17_sp9': 17.19e6, 'scale_sp9': 68.31e6, 's265_sp9': 265.1e6}
lams = [float(x) for x in sys.argv[1].split(',')]
tcs = sys.argv[2].split(',')
ROOT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'v0.4', 'tuning')
for tc in tcs:
    T, C = tc.split(':')
    for lam in lams:
        row = []
        for d in consts:
            tag = 't%s_c%s' % (T, C)
            p = os.path.join(ROOT, d, tag + '.stats.json')
            if not os.path.exists(p):
                row.append('%s:-' % d); continue
            st = json.load(open(p))['counters']; nq = st['queries']
            ck, cp, cd, cc, ch = consts[d]
            km, po, dg, hits = st['kmers_generated']/nq, st['postings_visited']/nq, st['diagonals_scored']/nq, st['hits']/nq
            best = None
            for r in json.load(open(os.path.join(ROOT, d, tag + '.gates.json')))['rows']:
                cost = (km*ck + po*cp + dg*cd + r['cells_per_q']*cc + hits*ch) * 1e-6
                v = r['auc1'] - lam*math.log(cost)
                if best is None or v > best[0]:
                    best = (v, r['gate'], cost, r['auc1'])
            row.append('%5.1fM G%d (%.1fs %.4f)' % (size[d]/1e6, best[1], best[2], best[3]))
        print('T%s C%s lam %.3f: ' % (T, C, lam) + '  '.join(row))
