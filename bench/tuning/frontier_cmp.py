#!/usr/bin/env python3
"""frontier_cmp.py c_k c_p c_d c_c c_h DIR1 DIR2 ... : best AUC1 (and config) within CPU budgets."""
import json, glob, sys, os
ck, cp, cd, cc, ch = (float(x) for x in sys.argv[1:6])
dirs = sys.argv[6:]
def pts(d):
    out = []
    for f in sorted(glob.glob(os.path.join(d, '*.gates.json'))):
        tag = os.path.basename(f).replace('.gates.json', '')
        T = int(tag.split('_')[0][1:]); C = int(tag.split('_')[1][1:])
        st = json.load(open(os.path.join(d, tag + '.stats.json')))['counters']
        nq = st['queries']
        km, po, dg, hits = st['kmers_generated']/nq, st['postings_visited']/nq, st['diagonals_scored']/nq, st['hits']/nq
        for r in json.load(open(f))['rows']:
            cost = (km*ck + po*cp + dg*cd + r['cells_per_q']*cc + hits*ch) * 1e-6  # CPU s per 1000 q
            out.append((cost, r['auc1'], 'T%dC%dG%d' % (T, C, r['gate']), r['fp_at_recall'].get('0.08')))
    return out
P = {d: pts(d) for d in dirs}
budgets = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100]
print('budget(s/1000q) ' + ' '.join('%-26s' % d for d in dirs))
for b in budgets:
    row = []
    for d in dirs:
        best = max((p for p in P[d] if p[0] <= b), key=lambda p: p[1], default=None)
        row.append('%-26s' % ('%.4f %s(%.1f)' % (best[1], best[2], best[0]) if best else '-'))
    print('%8.1f        ' % b + ' '.join(row))
