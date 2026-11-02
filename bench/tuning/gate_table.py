#!/usr/bin/env python3
"""gate_table.py c_k c_p c_d c_c c_h DIR T:C[,T:C...] [GLO-GHI]: per-gate cost/AUC1/FP@0.08 for configs."""
import json, sys, os
ck, cp, cd, cc, ch = (float(x) for x in sys.argv[1:6])
d = sys.argv[6]
glo, ghi = (int(x) for x in (sys.argv[8].split('-') if len(sys.argv) > 8 else ('26', '46')))
for tc in sys.argv[7].split(','):
    T, C = tc.split(':')
    tag = 't%s_c%s' % (T, C)
    st = json.load(open(os.path.join(d, tag + '.stats.json')))['counters']
    nq = st['queries']
    km, po, dg, hits = st['kmers_generated']/nq, st['postings_visited']/nq, st['diagonals_scored']/nq, st['hits']/nq
    rows = json.load(open(os.path.join(d, tag + '.gates.json')))['rows']
    print('%s T%s C%s:' % (d, T, C), '  '.join('G%d %.1fs %.4f' % (r['gate'], (km*ck + po*cp + dg*cd + r['cells_per_q']*cc + hits*ch)*1e-6, r['auc1']) for r in rows if glo <= r['gate'] <= ghi))
