#!/usr/bin/env python3
"""Markdown tables for README.md from run_suite.py results.

  readme_tables.py RESULTS_DIR > tables.md
  readme_tables.py RESULTS_DIR --fill README.md   # refresh the embedded tables

Reads RESULTS_DIR/<suite>/<run>__<workload>.json for the tune, eval, null,
scale and build suites (whatever exists) and prints every table (before/after,
tool comparison, presets, per-workload, null, scale, builds, ablations), or
rewrites the `<!-- table:NAME -->` … `<!-- /table -->` blocks of a Markdown
file in place.
"""
import json
import os
import sys

R = sys.argv[1]


def rec(suite, run, wl):
    p = os.path.join(R, suite, '%s__%s.json' % (run, wl))
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def thr(r, e):
    for t in r['evaluation']['thresholds']:
        if t['evalue'] == e:
            return t
    return None


def fp_at(r, lvl):
    for m in r['evaluation']['matched_recall']:
        if abs(m['recall'] - lvl) < 1e-9:
            return m['FP'] if m['reached'] else None
    return None


def wall(r):
    m = r['measure']['wall_s']
    if r['measure']['repeat'] > 1:
        return '%.2f (%.2f–%.2f)' % (m['median'], m['min'], m['max'])
    return '%.2f' % m['median']


def cg_peak(m):
    """Largest cgroup peak (MiB) over the repetitions of a measurement."""
    peaks = [x['cgroup']['max_usage_bytes'] or 0 for x in m['runs'] if x.get('cgroup')]
    return '%.0f' % (max(peaks) / 1048576) if peaks else '—'


def wall1(r):
    m = r['measure']['wall_s']
    if r['measure']['repeat'] > 1:
        return '%.1f (%.1f–%.1f)' % (m['median'], m['min'], m['max'])
    return '%.1f' % m['median']


def rss(r):
    return '%.0f' % (r['measure']['peak_rss_bytes']['max'] / 1048576)


def fmt(x, f='%d'):
    return '—' if x is None else f % x


def null_ratio(r, e):
    for t in r['null']['thresholds']:
        if t['evalue'] == e:
            return t
    return None


COLS = [('v02_default', 'v0.2'), ('st_s5.7', '**v0.3 default** (`-s 5.7`)'), ('st_s2', 'v0.3 `-s 2`'),
        ('st_s7.5', 'v0.3 `-s 7.5`'), ('mm_s5.7', 'MMseqs2 `-s 5.7`'), ('mm_s7.5', 'MMseqs2 `-s 7.5`'),
        ('dm_sensitive', 'DIAMOND `--sensitive`'), ('dm_ultra', 'DIAMOND `--ultra-sensitive`')]


def headline():
    rows = []
    cols = [(k, h) for k, h in COLS if rec('eval', k, 'eval')]
    if not cols:
        return ''
    head = '| | ' + ' | '.join(h for _, h in cols) + ' |'
    rows.append(head)
    rows.append('|---' * (len(cols) + 1) + '|')

    def line(label, f):
        vals = []
        for k, _ in cols:
            try:
                vals.append(f(k))
            except (TypeError, KeyError):
                vals.append('—')
        rows.append('| %s | %s |' % (label, ' | '.join(vals)))

    e = lambda k: rec('eval', k, 'eval')
    line('*held-out SCOPe40: 4,000 queries × 14,323 domains*', lambda k: '')
    line('wall time, s (median of 3)', lambda k: '%.2f' % e(k)['measure']['wall_s']['median'])
    line('peak RSS, MiB', lambda k: rss(e(k)))
    line('AUC1 (sensitivity to 1st FP)', lambda k: '%.4f' % e(k)['evaluation']['mean_auc1'])
    line('true positives at E ≤ 1e-3', lambda k: '%d' % thr(e(k), 1e-3)['TP'])
    line('false positives at E ≤ 1e-3', lambda k: '%d' % thr(e(k), 1e-3)['FP'])
    line('FP at matched recall 0.07', lambda k: fmt(fp_at(e(k), 0.07)))
    line('FP at matched recall 0.08', lambda k: fmt(fp_at(e(k), 0.08)))
    line('FP at matched recall 0.10', lambda k: fmt(fp_at(e(k), 0.10)))
    lc = lambda k: rec('eval', k, 'lc')
    line('*low-complexity control (1,000 LC-insert queries × SCOPe40 + decoys)*', lambda k: '')
    line('TP / FP at E ≤ 1e-3', lambda k: '%d / %d' % (thr(lc(k), 1e-3)['TP'], thr(lc(k), 1e-3)['FP']))
    nq = lambda k: rec('null', k, 'null_q')
    line('*null: 1,000 shuffled queries (expected chance hits per query = E)*', lambda k: '')
    line('chance hits per query at E ≤ 1', lambda k: '%.3f' % null_ratio(nq(k), 1.0)['per_query'])
    line('chance hits per query at E ≤ 0.01', lambda k: '%.4f' % null_ratio(nq(k), 0.01)['per_query'])
    sc = lambda k: rec('scale', k, 'scale')
    line('*68 M residues (374,323 seqs): 1,000 queries*', lambda k: '')
    line('wall time, s', lambda k: '%.1f' % sc(k)['measure']['wall_s']['median'])
    line('peak RSS, MiB', lambda k: rss(sc(k)))
    line('AUC1', lambda k: '%.4f' % sc(k)['evaluation']['mean_auc1'])
    line('TP / FP at E ≤ 1e-3', lambda k: '%d / %d' % (thr(sc(k), 1e-3)['TP'], thr(sc(k), 1e-3)['FP']))
    return '\n'.join(rows)


def before_after():
    """Concise v0.2 → v0.3 (default `-s 5.7`) table with the change."""
    rows = ['| | v0.2 | v0.3 (default) | change |', '|---|---|---|---|']
    E = lambda wl, k: rec('eval', k, wl)
    S = lambda k: rec('scale', k, 'scale')
    N = lambda k: rec('null', k, 'null_q')

    def add(label, old, new, change):
        rows.append('| %s | %s | %s | %s |' % (label, old, new, change))

    def ratio(o, n, more='faster', fewer=False):
        if o is None or n is None or n == 0:
            return ''
        return '%.1f× %s' % (o / n, more)

    def pct(o, n):
        return ('%+.0f%%' % (100.0 * (n - o) / o)).replace('-', '−') if o else ''

    o, n = E('eval', 'v02_default'), E('eval', 'st_s5.7')
    if o and n:
        rows.append('| *held-out SCOPe40: 4,000 queries × 14,323 domains (2 threads)* | | | |')
        wo, wn = o['measure']['wall_s']['median'], n['measure']['wall_s']['median']
        add('wall time, s', '%.1f' % wo, '%.1f' % wn, ratio(wo, wn))
        ro, rn = o['measure']['peak_rss_bytes']['max'] / 1048576, n['measure']['peak_rss_bytes']['max'] / 1048576
        add('peak RSS, MiB', '%.0f' % ro, '%.0f' % rn, ratio(ro, rn, 'less'))
        ao, an = o['evaluation']['mean_auc1'], n['evaluation']['mean_auc1']
        add('sensitivity (AUC1)', '%.3f' % ao, '%.3f' % an, pct(ao, an))
        to, tn = thr(o, 1e-3), thr(n, 1e-3)
        add('true positives at E ≤ 1e-3', '{:,}'.format(to['TP']), '{:,}'.format(tn['TP']), pct(to['TP'], tn['TP']))
        add('false positives at E ≤ 1e-3', '%d' % to['FP'], '%d' % tn['FP'], pct(to['FP'], tn['FP']))
        fo, fn = fp_at(o, 0.08), fp_at(n, 0.08)
        add('false positives at matched recall 0.08', '—' if fo is None else '{:,}'.format(fo),
            fmt(fn), ('%.0f× fewer' % (fo / max(fn, 1))) if fo is not None and fn is not None else '')
    o, n = E('lc', 'v02_default'), E('lc', 'st_s5.7')
    if o and n:
        add('low-complexity-insert queries: FP at E ≤ 1e-3', '{:,}'.format(thr(o, 1e-3)['FP']), '%d' % thr(n, 1e-3)['FP'], '')
    o, n = N('v02_default'), N('st_s5.7')
    if o and n:
        add('shuffled queries: chance hits per query at E ≤ 1 (≤ 1 = calibrated)', '%.2f' % null_ratio(o, 1.0)['per_query'],
            '%.2f' % null_ratio(n, 1.0)['per_query'], 'overconfident → conservative')
    o, n = S('v02_default'), S('st_s5.7')
    if o and n:
        rows.append('| *68 M residues (374,323 sequences): 1,000 queries* | | | |')
        wo, wn = o['measure']['wall_s']['median'], n['measure']['wall_s']['median']
        add('wall time, s', '%.1f' % wo, '%.1f' % wn, ratio(wo, wn))
        ro, rn = o['measure']['peak_rss_bytes']['max'] / 1048576, n['measure']['peak_rss_bytes']['max'] / 1048576
        add('peak RSS, MiB (v0.3: mostly the memory-mapped index)', '%.0f' % ro, '%.0f' % rn, ratio(ro, rn, 'less'))
        ao, an = o['evaluation']['mean_auc1'], n['evaluation']['mean_auc1']
        add('sensitivity (AUC1)', '%.3f' % ao, '%.3f' % an, pct(ao, an))
        add('true / false positives at E ≤ 1e-3', '%d / %d' % (thr(o, 1e-3)['TP'], thr(o, 1e-3)['FP']),
            '%d / %d' % (thr(n, 1e-3)['TP'], thr(n, 1e-3)['FP']), '')
        for k_old, k_new, label in (('v02_default_mem1g', 'st_s5.7_mem1g', 'on a 1 GiB machine: wall time, s'),
                                    (None, 'st_s5.7_mem512m', 'on a 512 MiB machine: wall time, s'),
                                    (None, 'st_s5.7_mem256m', 'on a 256 MiB machine: wall time, s')):
            rn = S(k_new)
            ro = S(k_old) if k_old else None
            if not rn:
                continue
            v02_1g = S('v02_default_mem1g')
            if ro and 'failed' in ro:
                old = 'killed (out of memory)'
            elif ro:
                old = '%.1f' % ro['measure']['wall_s']['median']
            elif v02_1g and 'failed' in v02_1g:
                old = '— (killed at 1 GiB)'
            else:
                old = 'not run'
            new = '%.1f' % rn['measure']['wall_s']['median'] if 'failed' not in rn else 'failed'
            sh = (rn.get('sabertooth_stats') or {}).get('shards', {}).get('count', 1)
            add(label, old, new + (' (%d shards)' % sh if sh > 1 else ''), '')
        b = n.get('db_build') or {}
        if b.get('createindex'):
            ci = b['createindex']
            add('seed index', 'rebuilt in every search', 'built once: %.0f s, %.0f MiB on disk' % (
                ci['wall_s']['median'], (b.get('index_bytes') or 0) / 1048576), '')
    return '\n'.join(rows)


def presets():
    rows = ['| `-s` | k-mer score / cap / gate | wall s | AUC1 | TP@1e-3 | FP@1e-3 | FP@E10 |', '|---|---|---|---|---|---|---|']
    table = {'1': '26 / 192 / 30', '2': '25 / 256 / 30', '4': '24 / 320 / 30', '5.7': '22 / 512 / 33',
             '6.5': '21 / 640 / 33', '7.5': '19 / 1024 / 34', '8.5': '17 / 1536 / 36', '9': '16 / 2048 / 38'}
    for s, p in table.items():
        r = rec('tune', 'st_s%s' % s, 'tune')
        if not r:
            continue
        t = thr(r, 1e-3)
        rows.append('| %s%s | %s | %s | %.4f | %d | %d | %d |' % (
            s, ' (default)' if s == '5.7' else '', p, wall(r), r['evaluation']['mean_auc1'], t['TP'], t['FP'],
            thr(r, 10.0)['FP']))
    return '\n'.join(rows)


def workload_table(suite, wl, runs, label, levels=(0.06, 0.07, 0.08, 0.10)):
    out = ['**%s**' % label, '', '| run | wall s | RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | '
           + ' | '.join('FP@rec %.2f' % x for x in levels) + ' |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for k in runs:
        r = rec(suite, k, wl)
        if not r:
            continue
        if 'failed' in r:
            out.append('| %s | failed (exit %s) | | | | | | | | | | |' % (k, r['failed']['exit']))
            continue
        t = thr(r, 1e-3)
        out.append('| %s | %s | %s | %.4f | %d | %d | %.4f | %d | %s | %s | %s | %s |' % (
            k, wall(r), rss(r), r['evaluation']['mean_auc1'], t['TP'], t['FP'], t['FDP'], thr(r, 10.0)['FP'],
            *[fmt(fp_at(r, x)) for x in levels]))
    return '\n'.join(out)


def null_table():
    out = ['| run | workload | hits E≤1e-3 | E≤0.01 | E≤0.1 | E≤1 | observed/expected at E≤0.1 | at E≤1 |', '|---|---|---|---|---|---|---|---|']
    for wl in ('null_q', 'null_db'):
        for k in ('v02_default', 'st_s5.7', 'st_s7.5', 'st_s5.7_plain', 'st_exhaustive', 'st_exhaustive_plain', 'mm_s5.7', 'dm_sensitive'):
            r = rec('null', k, wl)
            if not r or 'null' not in r:
                continue
            g = lambda e: null_ratio(r, e)
            out.append('| %s | %s | %d | %d | %d | %d | %.3f | %.3f |' % (
                k, wl, g(1e-3)['hits'], g(0.01)['hits'], g(0.1)['hits'], g(1.0)['hits'],
                g(0.1)['observed_over_expected'], g(1.0)['observed_over_expected']))
    return '\n'.join(out)


def scale_table():
    out = ['| run | threads | machine | shards | wall s | peak RSS MiB | cgroup peak MiB | AUC1 | TP@1e-3 | FP@1e-3 | FP@rec 0.06 | FP@rec 0.07 |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for k in ('st_s5.7', 'st_s5.7_t1', 'st_s5.7_cold', 'st_s5.7_mem1g', 'st_s5.7_mem1g_noshard', 'st_s5.7_mem512m',
              'st_s5.7_mem512m_noshard', 'st_s5.7_mem256m', 'st_s2', 'st_s7.5', 'st_s5.7_occ200', 'st_s5.7_occ50',
              'mm_s5.7', 'mm_s5.7_t1', 'mm_s5.7_noidx', 'mm_s5.7_mem1g', 'mm_s5.7_mem1g_split', 'dm_sensitive', 'dm_sensitive_t1',
              'dm_sensitive_mem1g', 'v02_default', 'v02_default_mem1g'):
        r = rec('scale', k, 'scale')
        if not r:
            continue
        opt = r.get('options', {})
        mem = ('%d MiB RAM' % (opt['cgroup'] / 1048576)) if opt.get('cgroup') else ('cold cache' if opt.get('evict_db') else '8 GiB, warm')
        if '--memory-limit' in r.get('args', []):
            mem += ', `--memory-limit %s`' % r['args'][r['args'].index('--memory-limit') + 1]
        if '--split-memory-limit' in r.get('args', []):
            mem += ', `--split-memory-limit %s`' % r['args'][r['args'].index('--split-memory-limit') + 1]
        if opt.get('no_index'):
            mem += ', no precomputed index'
        sh = (r.get('sabertooth_stats') or {}).get('shards', {}).get('count', 1) if r['tool'] == 'sabertooth' else '—'
        cgp = cg_peak(r['measure'])
        if 'failed' in r:
            out.append('| %s | %s | %s | %s | failed (exit %s) after %.1f | %s | %s | | | | | |' % (
                k, r.get('threads'), mem, sh, r['failed']['exit'], r['measure']['wall_s']['median'], rss(r), cgp))
            continue
        t = thr(r, 1e-3)
        out.append('| %s | %s | %s | %s | %s | %s | %s | %.4f | %d | %d | %s | %s |' % (
            k, r.get('threads'), mem, sh, wall1(r), rss(r), cgp, r['evaluation']['mean_auc1'], t['TP'], t['FP'],
            fmt(fp_at(r, 0.06)), fmt(fp_at(r, 0.07))))
    return '\n'.join(out)


def build_budgets():
    p = os.path.join(R, 'build', 'summary.json')
    if not os.path.exists(p):
        return ''
    d = json.load(open(p))
    out = ['| `createindex` (68 M residues) | machine | wall s | peak RSS MiB | cgroup peak MiB | partitions | index MiB |',
           '|---|---|---|---|---|---|---|']
    for r in d['runs']:
        lim = ' '.join(r['args'][1:]) if r['args'] else 'default'
        machine = ('%d MiB RAM' % (r['options']['cgroup'] >> 20)) if r['options'].get('cgroup') else '8 GiB'
        m = r['measure']
        cgp = cg_peak(m)
        if 'failed' in r:
            out.append('| `--memory-limit %s` | %s | failed (exit %s) | %.0f | %s | | |' % (lim, machine, r['failed']['exit'],
                                                                                       m['peak_rss_bytes']['max'] / 1048576, cgp))
            continue
        b = r['createindex_stats']['build']
        out.append('| %s | %s | %s | %.0f | %s | %d | %.0f |' % (
            '`--memory-limit %s`' % lim if r['args'] else 'default (1 GiB partition buffer)', machine, wall1(r),
            m['peak_rss_bytes']['max'] / 1048576, cgp, b['partitions'], r['index_bytes'] / 1048576))
    return '\n'.join(out)


def builds():
    out = ['| tool | database | build wall s | build peak RSS MiB | size on disk MiB |', '|---|---|---|---|---|']
    for suite, wl in (('eval', 'eval'), ('scale', 'scale')):
        for k in ('st_s5.7', 'mm_s5.7', 'dm_sensitive'):
            r = rec(suite, k, wl)
            if not r or not r.get('db_build'):
                continue
            b = r['db_build']
            steps = [s for s in ('makedb', 'createdb', 'createindex') if s in b]
            w = sum(b[s]['wall_s']['median'] for s in steps)
            m = max(b[s]['peak_rss_bytes']['max'] for s in steps) / 1048576
            out.append('| %s | %s | %.1f | %.0f | %.0f |' % (r['tool'], b['db'], w, m, (b.get('db_bytes') or 0) / 1048576))
    return '\n'.join(out)


EVAL_RUNS = ['v02_default', 'st_s2', 'st_s5.7', 'st_s7.5', 'mm_s5.7', 'mm_s7.5', 'dm_sensitive', 'dm_ultra', 'st_exhaustive']
EVAL_WL = {'eval': ('held-out homology (4,000 queries)', (0.06, 0.07, 0.08, 0.10)),
           'lc': ('low-complexity inserts vs SCOPe40 + decoys', (0.04, 0.06, 0.07, 0.08)),
           'frag': ('fragments (30–60% of the domain)', (0.03, 0.04, 0.05, 0.06)),
           'multi': ('two-domain fusions', (0.06, 0.07, 0.08, 0.10))}
ABLATIONS = ['st_s5.7', 'st_s5.7_nomaskscore', 'st_s5.7_nocbc', 'st_s5.7_nomask', 'st_s5.7_plain',
             'st_s5.7_kadane', 'st_s5.7_diag1', 'st_v02seeds']


def tables():
    """name -> Markdown table (the blocks README.md embeds)."""
    t = {'before_after': before_after(), 'headline': headline(), 'presets': presets(), 'null': null_table(),
         'scale': scale_table(), 'builds': builds(), 'build_budgets': build_budgets()}
    for wl, (label, lv) in EVAL_WL.items():
        t['eval_' + wl] = workload_table('eval', wl, EVAL_RUNS, label, lv)
    t['ablation_tune'] = workload_table('tune', 'tune', ABLATIONS, 'tuning split (1,000 queries vs SCOPe40)',
                                        (0.06, 0.07, 0.08, 0.10))
    t['ablation_lc'] = workload_table('tune', 'lc_tune', ABLATIONS, 'low-complexity tuning queries vs SCOPe40 + decoys',
                                      (0.04, 0.06, 0.07, 0.08))
    return t


def fill(path, t):
    """Replace every `<!-- table:NAME -->` … `<!-- /table -->` block of `path`."""
    import re
    text = open(path).read()

    def rep(m):
        name = m.group(1)
        if name not in t:
            sys.exit('unknown table %s in %s' % (name, path))
        return '<!-- table:%s -->\n%s\n<!-- /table -->' % (name, t[name])
    new, n = re.subn(r'<!-- table:([a-z_]+) -->.*?<!-- /table -->', rep, text, flags=re.S)
    with open(path, 'w') as f:
        f.write(new)
    print('filled %d table(s) in %s' % (n, path), file=sys.stderr)


if __name__ == '__main__':
    t = tables()
    if len(sys.argv) > 3 and sys.argv[2] == '--fill':
        for p in sys.argv[3:]:
            fill(p, t)
    else:
        for name, body in t.items():
            print('## %s\n\n%s\n' % (name, body))
