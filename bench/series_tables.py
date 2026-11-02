#!/usr/bin/env python3
"""Render the database-size series (run_series.py output) as Markdown tables.

  series_tables.py OUT [--md OUT/series.md]

Per size: every search configuration (median wall of the measured runs, CPU =
user + system, peak RSS of the largest process, AUC1, TP/FP at E <= 1e-3, FP at
matched recall), the database builds, a head-to-head line per sensitivity
level, and the scaling exponent of wall time with database size (least-squares
slope of log wall vs log residues over the sizes >= 17 M).
"""
import argparse
import glob
import json
import math
import os

ORDER = ['st_s1', 'st_s2', 'st_s4', 'st_s5.7', 'st_s7.5', 'st03_s5.7', 'mm_s1', 'mm_s4', 'mm_s5.7', 'mm_s7.5',
         'dm_sensitive', 'st_s5.7_t1', 'mm_s5.7_t1', 'st_s5.7_mem1g', 'st03_s5.7_mem1g', 'mm_s5.7_mem1g',
         'dm_sensitive_mem1g']
LABEL = {
    'st_s1': 'Sabertooth 0.4 `-s 1`', 'st_s2': 'Sabertooth 0.4 `-s 2`', 'st_s4': 'Sabertooth 0.4 `-s 4`',
    'st_s5.7': 'Sabertooth 0.4 default (`-s 5.7`)', 'st_s7.5': 'Sabertooth 0.4 `-s 7.5`',
    'st03_s5.7': 'Sabertooth 0.3 default', 'mm_s1': 'MMseqs2 `-s 1`', 'mm_s4': 'MMseqs2 `-s 4`',
    'mm_s5.7': 'MMseqs2 default (`-s 5.7`)', 'mm_s7.5': 'MMseqs2 `-s 7.5`', 'dm_sensitive': 'DIAMOND `--sensitive`',
    'st_s5.7_t1': 'Sabertooth 0.4 default, 1 thread', 'mm_s5.7_t1': 'MMseqs2 default, 1 thread',
    'st_s5.7_mem1g': 'Sabertooth 0.4 default, 1 GiB RAM', 'st03_s5.7_mem1g': 'Sabertooth 0.3 default, 1 GiB RAM',
    'mm_s5.7_mem1g': 'MMseqs2 default, 1 GiB RAM', 'dm_sensitive_mem1g': 'DIAMOND `--sensitive`, 1 GiB RAM',
}
TOOL = {'st': 'Sabertooth 0.4', 'st03': 'Sabertooth 0.3', 'mm': 'MMseqs2', 'dm': 'DIAMOND'}


def load(out):
    sizes = {}
    for d in sorted(glob.glob(os.path.join(out, '*M'))):
        size = int(os.path.basename(d)[:-1])
        runs = {}
        for p in glob.glob(os.path.join(d, '*.json')):
            b = os.path.basename(p)[:-5]
            if b.startswith('build_') or b.endswith('.measure') or b.endswith('.eval'):
                continue
            runs[b] = json.load(open(p))
        builds = {}
        for p in glob.glob(os.path.join(d, 'build_*.json')):
            builds[os.path.basename(p)[6:-5]] = json.load(open(p))
        sizes[size] = (runs, builds)
    return sizes


def row_vals(r):
    m = r.get('measure')
    if not m or m['runs'][-1]['exit'] != 0 or 'evaluation' not in r:
        return None
    ev = r['evaluation']
    t = {x['evalue']: x for x in ev['thresholds']}
    cpu = sorted(x['user_s'] + x['sys_s'] for x in m['runs'])[len(m['runs']) // 2]
    mr = {x['recall']: (x['FP'] if x['reached'] else None) for x in ev['matched_recall']}
    cg = [x.get('cgroup', {}).get('max_usage_bytes') for x in m['runs'] if x.get('cgroup')]
    return {
        'wall': m['wall_s']['median'], 'wmin': m['wall_s']['min'], 'wmax': m['wall_s']['max'], 'n': m['repeat'],
        'cpu': cpu, 'rss': m['peak_rss_bytes']['max'] / 1048576, 'auc1': ev['mean_auc1'],
        'tp3': t[1e-3]['TP'], 'fp3': t[1e-3]['FP'], 'fp_rec': mr, 'cgroup_peak': max(cg) / 1048576 if cg and all(cg) else None,
    }


def curve(runs, names):
    pts = []
    for n in names:
        if n in runs:
            v = row_vals(runs[n])
            if v:
                pts.append((v['auc1'], v['wall']))
    pts.sort()
    return pts


def interp(pts, a):
    """(text, seconds or None) for the wall time at AUC1 `a` along `pts`."""
    if a <= pts[0][0]:
        return '≤ %.2f' % pts[0][1], None
    if a > pts[-1][0]:
        return '> %.2f' % pts[-1][1], None
    for (a0, w0), (a1, w1) in zip(pts, pts[1:]):
        if a0 <= a <= a1:
            f = (a - a0) / (a1 - a0) if a1 > a0 else 0.0
            w = math.exp(math.log(w0) + f * (math.log(w1) - math.log(w0)))
            return '%.2f' % w, w
    return '–', None


def index_mib(b):
    """Size of Sabertooth's .stidx from createindex's report (other tools: –)."""
    import re
    st = b['steps'].get('createindex')
    if b['tool'] in ('st', 'st03') and st:
        m = re.search(r', ([0-9.]+) MiB\) in ', st['runs'][-1]['stderr_tail'])
        if m:
            return m.group(1)
    return '–'


def slabel(size):
    # series directories are named by rounded millions; the smallest is SCOPe40 alone (2.6 M)
    return '2.6 M' if size == 3 else '%d M' % size


def fmt_fp(v):
    return '–' if v is None else str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--md')
    ap.add_argument('--level', type=int, default=2, help='markdown level of the section headings (README: 4)')
    a = ap.parse_args()
    h = '#' * a.level
    S = load(a.out)
    md = []
    md.append('# Database-size series: Sabertooth vs MMseqs2 (and DIAMOND)\n')
    prov = os.path.join(a.out, 'provenance.json')
    if os.path.exists(prov):
        p = json.load(open(prov))
        md.append('Queries: `%s` (sha256 %s…), %s threads unless stated, median of %s measured runs after %s warm-up run(s); peak RSS: the largest of the measured runs.\n' % (
            os.path.basename(p['queries']['fasta']), p['queries']['sha256'][:12], p['threads'], p['repeat'], p['warmup']))
        md.append('| binary | version | sha256 |\n|---|---|---|')
        for k, b in sorted(p.get('binaries', {}).items()):
            md.append('| %s | %s | `%s…` |' % (TOOL.get(k, k), b['version'], b['sha256'][:16]))
        md.append('')
    for size in sorted(S):
        runs, builds = S[size]
        res = None
        for r in runs.values():
            ev = r.get('evaluation')
            if ev:
                res = ev.get('db_size')
                break
        md.append('%s %s\n' % (h, '2.6 M residues (SCOPe40 only)' if size == 3 else '%d M residues (SCOPe40 + filler)' % size))
        md.append('| configuration | wall s (min–max) | CPU s | peak RSS MiB | AUC1 | TP / FP at E ≤ 1e-3 | FP at recall 0.06 / 0.08 / 0.10 |')
        md.append('|---|---|---|---|---|---|---|')
        for name in ORDER:
            if name not in runs:
                continue
            v = row_vals(runs[name])
            if v is None:
                md.append('| %s | failed | | | | | |' % LABEL.get(name, name))
                continue
            rss = '%.0f' % v['rss'] + (' (cgroup peak %.0f)' % v['cgroup_peak'] if v['cgroup_peak'] else '')
            md.append('| %s | %.2f (%.2f–%.2f) | %.1f | %s | %.4f | %d / %d | %s / %s / %s |' % (
                LABEL.get(name, name), v['wall'], v['wmin'], v['wmax'], v['cpu'], rss, v['auc1'], v['tp3'], v['fp3'],
                fmt_fp(v['fp_rec'].get(0.06)), fmt_fp(v['fp_rec'].get(0.08)), fmt_fp(v['fp_rec'].get(0.1))))
        if builds:
            md.append('\n| database build (single run) | wall s | peak RSS MiB | on disk MiB (all files) | seed index MiB | note |')
            md.append('|---|---|---|---|---|---|')
            for t in ['st', 'st03', 'mm', 'dm']:
                if t in builds:
                    b = builds[t]
                    steps = ' + '.join('%s %.1f s' % (k, v['wall_s']['median']) for k, v in b['steps'].items() if v)
                    md.append('| %s | %.1f | %.0f | %.0f | %s | %s%s |' % (TOOL[t], b['wall_s'], b['peak_rss_bytes'] / 1048576,
                                                                     b['bytes'] / 1048576, index_mib(b), steps,
                                                                     '; index build failed' if b.get('index_failed') else ''))
        # head-to-head per sensitivity level
        pairs = [('st_s1', 'mm_s1'), ('st_s4', 'mm_s4'), ('st_s5.7', 'mm_s5.7'), ('st_s7.5', 'mm_s7.5')]
        lines = []
        for x, y in pairs:
            if x in runs and y in runs:
                vx, vy = row_vals(runs[x]), row_vals(runs[y])
                if vx and vy:
                    lines.append('| `-s %s` | %.2f× | %+.1f%% | %.2f× |' % (
                        x.split('_s')[1], vy['wall'] / vx['wall'], 100 * (vx['auc1'] / vy['auc1'] - 1), vy['rss'] / vx['rss']))
        if lines:
            md.append('\n| same `-s` | Sabertooth speed-up | Sabertooth AUC1 vs MMseqs2 | MMseqs2 / Sabertooth peak RSS |')
            md.append('|---|---|---|---|')
            md.extend(lines)
        md.append('')
    # matched sensitivity: time each tool needs to reach the other's AUC1
    # (log-linear interpolation of wall time along a tool's own -s presets)
    md.append('%s Matched sensitivity\n' % h)
    md.append('Wall time each tool needs to reach the other default\'s AUC1, interpolated (log wall time linear in AUC1) '
              'between that tool\'s measured `-s` settings; “>” = not reached by its most sensitive measured setting, '
              '“≤” = already exceeded by its fastest.\n')
    md.append('| database | MMseqs2 `-s 5.7`: AUC1 / wall s | Sabertooth 0.4 at that AUC1, s | | '
              'Sabertooth 0.4 default: AUC1 / wall s | MMseqs2 at that AUC1, s | |')
    md.append('|---|---|---|---|---|---|---|')
    for size in sorted(S):
        runs = S[size][0]
        st = curve(runs, ['st_s1', 'st_s2', 'st_s4', 'st_s5.7', 'st_s7.5'])
        mm = curve(runs, ['mm_s1', 'mm_s4', 'mm_s5.7', 'mm_s7.5'])
        ra, rb = runs.get('mm_s5.7'), runs.get('st_s5.7')
        va, vb = (row_vals(ra) if ra else None), (row_vals(rb) if rb else None)
        if not (va and vb and st and mm):
            continue
        t1 = interp(st, va['auc1'])
        t2 = interp(mm, vb['auc1'])
        md.append('| %s | %.4f / %.2f | %s | %s | %.4f / %.2f | %s | %s |' % (
            slabel(size), va['auc1'], va['wall'], t1[0],
            ('Sabertooth %.1f× faster' % (va['wall'] / t1[1])) if t1[1] else '',
            vb['auc1'], vb['wall'], t2[0], ('MMseqs2 %.1f× slower' % (t2[1] / vb['wall'])) if t2[1] else ''))
    md.append('')
    # scaling
    big = [s for s in sorted(S) if s >= 17]
    if len(big) >= 2:
        md.append('%s Scaling with database size\n' % h)
        md.append('Wall time (s) per size; exponent b of wall ∝ residues^b (least squares over %s M).\n' % ', '.join(map(str, big)))
        md.append('| configuration | ' + ' | '.join(slabel(s) for s in sorted(S)) + ' | b |')
        md.append('|---|' + '---|' * (len(S) + 1))
        for name in ORDER:
            vals = {}
            for s in sorted(S):
                r = S[s][0].get(name)
                v = row_vals(r) if r else None
                if v:
                    vals[s] = v['wall']
            if len(vals) < 2:
                continue
            pts = [(math.log(s), math.log(vals[s])) for s in big if s in vals]
            b = ''
            if len(pts) >= 2:
                mx = sum(p[0] for p in pts) / len(pts)
                my = sum(p[1] for p in pts) / len(pts)
                b = '%.2f' % (sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts))
            md.append('| %s | ' % LABEL.get(name, name) + ' | '.join('%.2f' % vals[s] if s in vals else '' for s in sorted(S)) + ' | %s |' % b)
        md.append('')
        md.append('| configuration | ' + ' | '.join(slabel(s) for s in sorted(S)) + ' |')
        md.append('|---|' + '---|' * len(S))
        for metric, key, f in [('AUC1', 'auc1', '%.4f'), ('peak RSS MiB', 'rss', '%.0f')]:
            for name in ['st_s5.7', 'st03_s5.7', 'mm_s5.7', 'st_s7.5', 'mm_s7.5', 'dm_sensitive']:
                vals = {}
                for s in sorted(S):
                    r = S[s][0].get(name)
                    v = row_vals(r) if r else None
                    if v:
                        vals[s] = v[key]
                if vals:
                    md.append('| %s — %s | ' % (LABEL.get(name, name), metric) + ' | '.join(f % vals[s] if s in vals else '' for s in sorted(S)) + ' |')
        md.append('')
    text = '\n'.join(md)
    if a.md:
        with open(a.md, 'w') as f:
            f.write(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
