#!/usr/bin/env python3
"""gate_opt.py c_k c_p c_d c_c c_h DIR [DIR...]: per (T,C), gate closest to the DB's frontier.

For each configuration and gate: modeled CPU (s/1000 q), AUC1, and the deficit
to the frontier (best AUC1 of any configuration at no more CPU). Prints per
(T,C) the gate band within 0.002 of the minimum deficit and the best gate."""
import json, glob, sys, os
ck, cp, cd, cc, ch = (float(x) for x in sys.argv[1:6])
def load(d):
    pts = {}
    for f in sorted(glob.glob(os.path.join(d, '*.gates.json'))):
        tag = os.path.basename(f).replace('.gates.json', '')
        T = int(tag.split('_')[0][1:]); C = int(tag.split('_')[1][1:])
        st = json.load(open(os.path.join(d, tag + '.stats.json')))['counters']
        nq = st['queries']
        km, po, dg, hits = st['kmers_generated']/nq, st['postings_visited']/nq, st['diagonals_scored']/nq, st['hits']/nq
        for r in json.load(open(f))['rows']:
            cost = (km*ck + po*cp + dg*cd + r['cells_per_q']*cc + hits*ch) * 1e-6
            pts.setdefault((T, C), []).append((r['gate'], cost, r['auc1']))
    return pts
res = {}
for d in sys.argv[6:]:
    P = load(d)
    allp = sorted((c, a) for v in P.values() for (_, c, a) in v)
    fr = []  # running max
    m = -1
    for c, a in allp:
        m = max(m, a); fr.append((c, m))
    import bisect
    cs = [c for c, _ in fr]
    def front(c):
        i = bisect.bisect_right(cs, c) - 1
        return fr[i][1] if i >= 0 else 0
    res[d] = {}
    for tc, v in sorted(P.items()):
        de = [(front(c) - a, g, c, a) for g, c, a in v]
        best = min(de)
        band = [g for (x, g, c, a) in de if x <= best[0] + 0.002]
        res[d][tc] = (best[1], min(band), max(band), best[2], best[3], best[0])
dirs = sys.argv[6:]
print('T  C     ' + '   '.join('%s: best[band] cpu auc def' % d for d in dirs))
for tc in sorted(set().union(*[set(r) for r in res.values()])):
    row = []
    for d in dirs:
        r = res[d].get(tc)
        row.append('G%2d[%2d-%2d] %6.1f %.4f %.4f' % r if r else '%-32s' % '-')
    print('T%2d C%4d  ' % tc + '   '.join(row))
