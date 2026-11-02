#!/usr/bin/env python3
"""Markdown tables of the v0.4 benchmarks (README.md / bench/results/v0.4).

  readme_tables_v04.py --v04 SUITES/v0.4 --v03 SUITES/v0.3 --series SERIES [--fill README.md]

SUITES/<version>/<suite>/<run>__<workload>.json come from run_suite.py
(run_v04_suites.sh); SERIES is run_series.py's output. Prints every table,
or rewrites the `<!-- table:NAME -->` … `<!-- /table -->` blocks of a Markdown
file in place.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import series_tables as ST  # noqa: E402

A = None


def rec(ver, suite, run, wl):
    root = A.v04 if ver == 4 else A.v03
    p = os.path.join(root, suite, '%s__%s.json' % (run, wl))
    if not os.path.exists(p):
        return None
    r = json.load(open(p))
    return None if 'failed' in r else r


def thr(r, e):
    return next((t for t in r['evaluation']['thresholds'] if t['evalue'] == e), None)


def fp_at(r, lvl):
    for m in r['evaluation']['matched_recall']:
        if abs(m['recall'] - lvl) < 1e-9:
            return m['FP'] if m['reached'] else None
    return None


def fmt(x, f='%d'):
    return '—' if x is None else f % x


def med(r):
    return r['measure']['wall_s']['median']


def wall(r, f='%.2f'):
    m = r['measure']['wall_s']
    if r['measure']['repeat'] > 1:
        return (f + ' (' + f + '–' + f + ')') % (m['median'], m['min'], m['max'])
    return f % m['median']


def rss(r):
    return r['measure']['peak_rss_bytes']['max'] / 1048576


def null_row(r, e):
    return next((t for t in r['null']['thresholds'] if t['evalue'] == e), None)


def cg_peak(r):
    peaks = [x['cgroup']['max_usage_bytes'] or 0 for x in r['measure']['runs'] if x.get('cgroup')]
    return '%.0f' % (max(peaks) / 1048576) if peaks else '—'


def ratio(o, n, word):
    return '%.1f× %s' % (o / n, word) if o and n else ''


def pct(o, n):
    if not o:
        return ''
    d = 100.0 * (n - o) / o
    return (('%+.1f%%' if abs(d) < 10 else '%+.0f%%') % d).replace('-', '−')


def fewer_more(o, n):
    """Change of a false-positive count (regressions in bold)."""
    if o is None or n is None:
        return ''
    if n == o:
        return 'same'
    if n < o:
        return ('%.1f× fewer' % (o / n)) if n and o / n >= 1.5 else '%d%% fewer' % round(100.0 * (o - n) / o)
    return '**%d%% more**' % round(100.0 * (n - o) / o)


def series_run(size, name):
    d = os.path.join(A.series, '%dM' % size, name + '.json')
    if not os.path.exists(d):
        return None
    r = json.load(open(d))
    return r if r.get('evaluation') else None


def series_build(size, tool):
    d = os.path.join(A.series, '%dM' % size, 'build_%s.json' % tool)
    return json.load(open(d)) if os.path.exists(d) else None


def before_after():
    rows = ['| | v0.3.0 | v0.4.0 | change |', '|---|---|---|---|']

    def add(label, old, new, change=''):
        rows.append('| %s | %s | %s | %s |' % (label, old, new, change))

    o, n = rec(3, 'eval', 'st_s5.7', 'eval'), rec(4, 'eval', 'st_s5.7', 'eval')
    if o and n:
        rows.append('| *held-out SCOPe40 (2.6 M residues): 4,000 queries, 2 threads* | | | |')
        add('wall time, s', '%.2f' % med(o), '%.2f' % med(n), ratio(med(o), med(n), 'faster'))
        add('peak RSS, MiB', '%.0f' % rss(o), '%.0f' % rss(n), ratio(rss(o), rss(n), 'less'))
        ao, an = o['evaluation']['mean_auc1'], n['evaluation']['mean_auc1']
        add('sensitivity (AUC1)', '%.4f' % ao, '%.4f' % an, pct(ao, an))
        to, tn = thr(o, 1e-3), thr(n, 1e-3)
        add('true / false positives at E ≤ 1e-3', '{:,} / {}'.format(to['TP'], to['FP']), '{:,} / {}'.format(tn['TP'], tn['FP']),
            pct(to['TP'], tn['TP']) + ' TP')
        for lvl in (0.08, 0.10):
            fo, fn = fp_at(o, lvl), fp_at(n, lvl)
            add('false positives at matched recall %.2f' % lvl, fmt(fo), fmt(fn), fewer_more(fo, fn))
    o, n = rec(3, 'eval', 'st_s5.7', 'lc'), rec(4, 'eval', 'st_s5.7', 'lc')
    if o and n:
        add('low-complexity-insert queries: TP / FP at E ≤ 1e-3', '%d / %d' % (thr(o, 1e-3)['TP'], thr(o, 1e-3)['FP']),
            '%d / %d' % (thr(n, 1e-3)['TP'], thr(n, 1e-3)['FP']))
    o, n = rec(3, 'null', 'st_s5.7', 'null_q'), rec(4, 'null', 'st_s5.7', 'null_q')
    if o and n:
        add('shuffled queries: chance hits per query at E ≤ 1 (≤ 1 = calibrated)', '%.3f' % null_row(o, 1.0)['per_query'],
            '%.3f' % null_row(n, 1.0)['per_query'])
    for size, label in ((68, '68 M residues (374,323 sequences)'), (265, '265 M residues (1,454,323 sequences)')):
        o, n = series_run(size, 'st03_s5.7'), series_run(size, 'st_s5.7')
        if not (o and n):
            continue
        rows.append('| *%s: 1,000 held-out queries* | | | |' % label)
        add('wall time, s', '%.1f' % med(o), '%.1f' % med(n), ratio(med(o), med(n), 'faster'))
        add('peak RSS, MiB', '%.0f' % rss(o), '%.0f' % rss(n), ratio(rss(o), rss(n), 'less'))
        ao, an = o['evaluation']['mean_auc1'], n['evaluation']['mean_auc1']
        add('sensitivity (AUC1)', '%.4f' % ao, '%.4f' % an, pct(ao, an))
        add('true / false positives at E ≤ 1e-3', '%d / %d' % (thr(o, 1e-3)['TP'], thr(o, 1e-3)['FP']),
            '%d / %d' % (thr(n, 1e-3)['TP'], thr(n, 1e-3)['FP']))
        fo, fn = fp_at(o, 0.10), fp_at(n, 0.10)
        add('false positives at matched recall 0.10', fmt(fo), fmt(fn), fewer_more(fo, fn))
        bo, bn = series_build(size, 'st03'), series_build(size, 'st')
        if bo and bn:
            add('makedb + createindex (single build): s / peak RSS MiB / on disk MiB',
                '%.1f / %.0f / %.0f' % (bo['wall_s'], bo['peak_rss_bytes'] / 1048576, bo['bytes'] / 1048576),
                '%.1f / %.0f / %.0f' % (bn['wall_s'], bn['peak_rss_bytes'] / 1048576, bn['bytes'] / 1048576),
                ratio(bo['bytes'], bn['bytes'], 'smaller'))
        if size == 68:
            o1, n1 = series_run(68, 'st03_s5.7_mem1g'), series_run(68, 'st_s5.7_mem1g')
            if o1 and n1:
                add('on a 1 GiB machine (cgroup): wall time, s', '%.1f' % med(o1), '%.1f' % med(n1), ratio(med(o1), med(n1), 'faster'))
    return '\n'.join(rows)


def head_to_head():
    """Default vs default (and -s 7.5 vs -s 7.5) at every database size."""
    sizes = sorted(int(os.path.basename(d)[:-1]) for d in os.listdir(A.series) if d.endswith('M') and os.path.isdir(os.path.join(A.series, d)))
    rows = ['| database | Sabertooth 0.4 default: wall s / RSS MiB / AUC1 | MMseqs2 `-s 5.7`: wall s / RSS MiB / AUC1 | speed-up | AUC1 | '
            'Sabertooth 0.4 `-s 7.5` | MMseqs2 `-s 7.5` | speed-up | AUC1 |', '|---|---|---|---|---|---|---|---|---|']
    rows[0] = rows[0].replace('| speed-up | AUC1 |', '| Sabertooth speed-up | Sabertooth AUC1 |')
    for s in sizes:
        cells = []
        for a, b in (('st_s5.7', 'mm_s5.7'), ('st_s7.5', 'mm_s7.5')):
            x, y = series_run(s, a), series_run(s, b)
            if x and y:
                cells += ['%.2f / %.0f / %.4f' % (med(x), rss(x), x['evaluation']['mean_auc1']),
                          '%.2f / %.0f / %.4f' % (med(y), rss(y), y['evaluation']['mean_auc1']),
                          '%.2f×' % (med(y) / med(x)), pct(y['evaluation']['mean_auc1'], x['evaluation']['mean_auc1'])]
            else:
                cells += ['', '', '', '']
        label = '2.6 M (SCOPe40)' if s == 3 else '%d M' % s
        rows.append('| %s | %s |' % (label, ' | '.join(cells)))
    return '\n'.join(rows)


HEAD = [('st03', 'st_s5.7', 'v0.3 default'), (4, 'st_s2', 'v0.4 `-s 2`'), (4, 'st_s5.7', '**v0.4 default**'),
        (4, 'st_s7.5', 'v0.4 `-s 7.5`'), (4, 'mm_s5.7', 'MMseqs2 `-s 5.7`'), (4, 'mm_s7.5', 'MMseqs2 `-s 7.5`'),
        (4, 'dm_sensitive', 'DIAMOND `--sensitive`'), (4, 'dm_ultra', 'DIAMOND `--ultra-sensitive`'),
        (4, 'st_exhaustive', 'exhaustive SW')]


def headline():
    cols = [(3 if v == 'st03' else 4, k, h) for v, k, h in HEAD]
    cols = [c for c in cols if rec(c[0], 'eval', c[1], 'eval')]
    rows = ['| | ' + ' | '.join(h for _, _, h in cols) + ' |', '|---' * (len(cols) + 1) + '|']

    def line(label, f):
        vals = []
        for v, k, _ in cols:
            try:
                vals.append(f(v, k))
            except (TypeError, KeyError, AttributeError):
                vals.append('—')
        rows.append('| %s | %s |' % (label, ' | '.join(vals)))

    e = lambda v, k: rec(v, 'eval', k, 'eval')
    line('*held-out SCOPe40: 4,000 queries × 14,323 domains*', lambda v, k: '')
    line('wall time, s (median of 3; exhaustive SW: 1 run)', lambda v, k: '%.2f' % med(e(v, k)))
    line('peak RSS, MiB', lambda v, k: '%.0f' % rss(e(v, k)))
    line('AUC1 (sensitivity to 1st FP)', lambda v, k: '%.4f' % e(v, k)['evaluation']['mean_auc1'])
    line('true positives at E ≤ 1e-3', lambda v, k: '%d' % thr(e(v, k), 1e-3)['TP'])
    line('false positives at E ≤ 1e-3', lambda v, k: '%d' % thr(e(v, k), 1e-3)['FP'])
    for lvl in (0.07, 0.08, 0.10):
        line('FP at matched recall %.2f' % lvl, lambda v, k, lvl=lvl: fmt(fp_at(e(v, k), lvl)))
    lc = lambda v, k: rec(v, 'eval', k, 'lc')
    line('*low-complexity control (1,000 LC-insert queries × SCOPe40 + decoys)*', lambda v, k: '')
    line('TP / FP at E ≤ 1e-3', lambda v, k: '%d / %d' % (thr(lc(v, k), 1e-3)['TP'], thr(lc(v, k), 1e-3)['FP']))
    nq = lambda v, k: rec(v, 'null', k, 'null_q')
    line('*null: 1,000 shuffled queries (expected chance hits per query = E)*', lambda v, k: '')
    line('chance hits per query at E ≤ 1', lambda v, k: '%.3f' % null_row(nq(v, k), 1.0)['per_query'])
    line('chance hits per query at E ≤ 0.01', lambda v, k: '%.4f' % null_row(nq(v, k), 0.01)['per_query'])
    return '\n'.join(rows)


WL = {'eval': ('held-out homology (4,000 queries; median of 3 runs after a warm-up, exhaustive SW one run)', (0.06, 0.07, 0.08, 0.10)),
      'lc': ('low-complexity inserts vs SCOPe40 + decoys (single runs)', (0.04, 0.06, 0.07, 0.08)),
      'frag': ('fragments, 30–60% of the domain (single runs)', (0.03, 0.04, 0.05, 0.06)),
      'multi': ('two-domain fusions (single runs)', (0.06, 0.07, 0.08, 0.10))}
WL_RUNS = [(3, 'st_s5.7', 'v0.3 st_s5.7'), (3, 'st_s7.5', 'v0.3 st_s7.5'), (4, 'st_s2', 'st_s2'), (4, 'st_s5.7', 'st_s5.7'),
           (4, 'st_s7.5', 'st_s7.5'), (4, 'mm_s5.7', 'mm_s5.7'), (4, 'mm_s7.5', 'mm_s7.5'), (4, 'dm_sensitive', 'dm_sensitive'),
           (4, 'dm_ultra', 'dm_ultra'), (4, 'st_exhaustive', 'st_exhaustive')]


def workload_table(wl):
    label, levels = WL[wl]
    out = ['**%s**' % label, '', '| run | wall s | RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | '
           + ' | '.join('FP@rec %.2f' % x for x in levels) + ' |', '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for v, k, name in WL_RUNS:
        r = rec(v, 'eval', k, wl)
        if not r:
            continue
        t = thr(r, 1e-3)
        out.append('| %s | %s | %.0f | %.4f | %d | %d | %.4f | %d | %s |' % (
            name, wall(r), rss(r), r['evaluation']['mean_auc1'], t['TP'], t['FP'], t['FDP'], thr(r, 10.0)['FP'],
            ' | '.join(fmt(fp_at(r, x)) for x in levels)))
    return '\n'.join(out)


def null_table():
    out = ['| run | workload | hits E≤1e-3 | E≤0.01 | E≤0.1 | E≤1 | observed/expected at E≤0.1 | at E≤1 |', '|---|---|---|---|---|---|---|---|']
    for wl in ('null_q', 'null_db'):
        for v, k, name in [(3, 'st_s5.7', 'v0.3 st_s5.7'), (3, 'st_s7.5', 'v0.3 st_s7.5'), (4, 'st_s5.7', 'st_s5.7'),
                           (4, 'st_s7.5', 'st_s7.5'), (4, 'st_s5.7_plain', 'st_s5.7_plain'), (4, 'st_exhaustive', 'st_exhaustive'),
                           (4, 'st_exhaustive_plain', 'st_exhaustive_plain'), (4, 'mm_s5.7', 'mm_s5.7'), (4, 'dm_sensitive', 'dm_sensitive')]:
            r = rec(v, 'null', k, wl)
            if not r or 'null' not in r:
                continue
            g = lambda e: null_row(r, e)
            out.append('| %s | %s | %d | %d | %d | %d | %.3f | %.3f |' % (
                name, wl, g(1e-3)['hits'], g(0.01)['hits'], g(0.1)['hits'], g(1.0)['hits'],
                g(0.1)['observed_over_expected'], g(1.0)['observed_over_expected']))
    return '\n'.join(out)


PRESETS = {'1': '32 / 64', '2': '30 / 128', '4': '26 / 256', '5.7': '22 / 512', '6.5': '20 / 512', '7.5': '18 / 1024',
           '8.5': '18 / 2048', '9': '16 / 4096'}


def presets():
    rows = ['| `-s` | v0.4 k-mer score / cap | SCOPe40: wall s | AUC1 | TP / FP@1e-3 | v0.3 same `-s`: wall s | AUC1 | '
            '68 M: wall s | AUC1 | TP / FP@1e-3 | v0.3 same `-s`: wall s | AUC1 |', '|---' * 12 + '|']
    for s, p in PRESETS.items():
        cells = []
        for wl in ('tune', 'scale'):
            n, o = rec(4, 'presets', 'st_s' + s, wl), rec(3, 'presets', 'st_s' + s, wl)
            if n:
                t = thr(n, 1e-3)
                cells += ['%.2f' % med(n), '%.4f' % n['evaluation']['mean_auc1'], '%d / %d' % (t['TP'], t['FP'])]
            else:
                cells += ['', '', '']
            cells += ['%.2f' % med(o), '%.4f' % o['evaluation']['mean_auc1']] if o else ['', '']
        rows.append('| %s%s | %s | %s |' % (s, ' (default)' if s == '5.7' else '', p, ' | '.join(cells)))
    return '\n'.join(rows)


def scale_limits():
    out = ['| run (68 M residues) | machine | shards | wall s | peak RSS MiB | cgroup peak MiB | AUC1 |', '|---|---|---|---|---|---|---|']
    items = [(4, 'st_s5.7', 'v0.4 default'), (3, 'st_s5.7', 'v0.3 default')]
    for v, k, name in items:
        for run, machine in (('st_s5.7_cold', 'cold page cache'), ('st_s5.7_mem512m', '512 MiB RAM, `--memory-limit 384M`'),
                             ('st_s5.7_mem256m', '256 MiB RAM, `--memory-limit 192M`')):
            r = rec(v, 'scale', run, 'scale')
            if not r:
                continue
            sh = (r.get('sabertooth_stats') or {}).get('shards', {}).get('count', 1)
            out.append('| %s | %s | %s | %s | %.0f | %s | %.4f |' % (name, machine, sh, wall(r, '%.1f'), rss(r), cg_peak(r),
                                                                    r['evaluation']['mean_auc1']))
    return '\n'.join(out)


def mem1g():
    """The 1 GiB machine at 68 M residues (database-size series, host A, held-out queries)."""
    out = ['| run (68 M residues, 1 GiB RAM) | memory option | shards | wall s | peak RSS MiB | cgroup peak MiB | AUC1 |',
           '|---|---|---|---|---|---|---|']
    for name, label in (('st_s5.7_mem1g', 'v0.4 default'), ('st03_s5.7_mem1g', 'v0.3 default'), ('mm_s5.7_mem1g', 'MMseqs2 `-s 5.7`'),
                        ('dm_sensitive_mem1g', 'DIAMOND `--sensitive`')):
        r = series_run(68, name)
        if not r:
            p = os.path.join(A.series, '68M', name + '.json')
            if os.path.exists(p):
                out.append('| %s | | | failed | | | |' % label)
            continue
        v = ST.row_vals(r)
        sh = '—'
        if name.startswith('st'):
            m = re.search(r'searched as (\d+) residue-balanced shard', r['measure']['runs'][-1]['stderr_tail'])
            sh = m.group(1) if m else '1'
        out.append('| %s | %s | %s | %.1f (%.1f–%.1f) | %.0f | %s | %.4f |' % (
            label, '`--memory-limit 1G`' if name.startswith('st') else 'none (see text)', sh, v['wall'], v['wmin'], v['wmax'],
            v['rss'], '%.0f' % v['cgroup_peak'] if v['cgroup_peak'] else '—', v['auc1']))
    return '\n'.join(out)


def tables():
    t = {'before_after_v04': before_after(), 'head_to_head': head_to_head(), 'headline_v04': headline(),
         'null_v04': null_table(), 'presets_v04': presets(), 'scale_limits_v04': scale_limits(), 'mem1g_v04': mem1g()}
    for wl in WL:
        t['eval_%s_v04' % wl] = workload_table(wl)
    md = []
    old = sys.stdout
    try:
        import io
        buf = io.StringIO()
        sys.stdout = buf
        sys.argv = ['series_tables.py', A.series, '--level', '4']
        ST.main()
        md = buf.getvalue()
    finally:
        sys.stdout = old
    t['series_v04'] = md.split('\n', 1)[1].strip() if md else ''
    return t


def fill(path, t):
    text = open(path).read()

    def rep(m):
        name = m.group(1)
        return '<!-- table:%s -->\n%s\n<!-- /table -->' % (name, t[name]) if name in t else m.group(0)
    text = re.sub(r'<!-- table:(\S+) -->\n.*?<!-- /table -->', rep, text, flags=re.S)
    open(path, 'w').write(text)


def main():
    global A
    ap = argparse.ArgumentParser()
    ap.add_argument('--v04', required=True)
    ap.add_argument('--v03', required=True)
    ap.add_argument('--series', required=True)
    ap.add_argument('--fill', nargs='*', default=[])
    A = ap.parse_args()
    t = tables()
    if A.fill:
        for p in A.fill:
            fill(p, t)
    else:
        for k, v in t.items():
            print('<!-- table:%s -->\n%s\n<!-- /table -->\n' % (k, v))


if __name__ == '__main__':
    main()
