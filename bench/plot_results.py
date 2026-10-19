#!/usr/bin/env python3
"""Figures and compact tables from run_suite.py results.

  plot_results.py RESULTS_DIR OUT_DIR

For every suite found in RESULTS_DIR (tune, eval, null, scale) writes
  OUT_DIR/<suite>_<workload>_roc.png   true positives vs false positives
                                       (log x) from each run's PR sweep
  OUT_DIR/<suite>.json                 compact per-run metrics (machine-readable)
and prints Markdown tables. Requires matplotlib for the figures (tables are
written without it).
"""
import json
import os
import sys


def load(results, suite):
    p = os.path.join(results, suite, 'summary.json')
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def compact(rec):
    m = rec['measure']
    out = {'run': rec['run'], 'tool': rec['tool'], 'workload': rec['workload'], 'threads': rec.get('threads'),
           'options': rec.get('options', {}), 'wall_s_median': m['wall_s']['median'],
           'wall_s_min': m['wall_s']['min'], 'wall_s_max': m['wall_s']['max'], 'repeats': m['repeat'],
           'peak_rss_mib': m['peak_rss_bytes']['max'] / 1048576}
    cgs = [x['cgroup'] for x in m['runs'] if x.get('cgroup')]
    if cgs:
        out['cgroup_limit_mib'] = cgs[0]['limit_bytes'] / 1048576
        out['cgroup_peak_mib'] = max((c['max_usage_bytes'] or 0) for c in cgs) / 1048576
    if 'failed' in rec:
        out['failed_exit'] = rec['failed']['exit']
        return out
    if rec.get('db_build'):
        b = rec['db_build']
        steps = [k for k in ('makedb', 'createdb', 'createindex') if k in b]
        out['db_build_wall_s'] = sum(b[k]['wall_s']['median'] for k in steps)
        out['db_build_peak_rss_mib'] = max(b[k]['peak_rss_bytes']['max'] for k in steps) / 1048576
        out['db_bytes'] = b.get('db_bytes')
    if 'evaluation' in rec:
        e = rec['evaluation']
        t = {r['evalue']: r for r in e['thresholds']}
        out.update({'auc1': e['mean_auc1'], 'denominators': e['denominators'],
                    'at_1e-3': {k: t[1e-3][k] for k in ('TP', 'FP', 'recall', 'recall_superfamily_only', 'precision', 'FDP', 'FPR')},
                    'at_10': {k: t[10.0][k] for k in ('TP', 'FP', 'recall', 'precision', 'FDP', 'FPR')},
                    'fp_at_matched_recall': {str(x['recall']): x['FP'] for x in e['matched_recall']},
                    'auc1_by_query_length': e['by_query_length']})
    if 'null' in rec:
        out['null'] = rec['null']
    if 'sabertooth_stats' in rec:
        c = rec['sabertooth_stats']['counters']
        out['stage_cpu_s'] = {k[6:]: c[k] for k in c if k.startswith('cpu_s_')}
        out['candidates_per_query'] = c['candidates'] / max(1, c['queries'])
    return out


def plot(suite, data, outdir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available: no figures', file=sys.stderr)
        return []
    files = []
    by_wl = {}
    for r in data['runs']:
        if 'evaluation' in r:
            by_wl.setdefault(r['workload'], []).append(r)
    for wl, rs in by_wl.items():
        fig, ax = plt.subplots(figsize=(7.5, 5))
        for r in rs:
            pr = r['evaluation']['pr_curve']
            xs = [max(p['FP'], 0.8) for p in pr]
            ys = [p['TP'] for p in pr]
            style = '-' if r['tool'] == 'sabertooth' else ('--' if r['tool'] in ('mmseqs', 'diamond') else ':')
            ax.plot(xs, ys, style, label='%s (AUC1 %.3f)' % (r['run'], r['evaluation']['mean_auc1']), linewidth=1.4)
        ax.set_xscale('log')
        ax.set_xlabel('false positives (cumulative, E-value sweep)')
        ax.set_ylabel('true positives')
        d = rs[0]['evaluation']['denominators']
        ax.set_title('%s / %s  (TP pairs %d, FP pairs %.3g)' % (suite, wl, d['TP_pairs'], d['FP_pairs']))
        ax.grid(True, which='both', alpha=0.3)
        ax.legend(fontsize=7, loc='lower right')
        fig.tight_layout()
        p = os.path.join(outdir, '%s_%s_roc.png' % (suite, wl))
        fig.savefig(p, dpi=130)
        plt.close(fig)
        files.append(p)
    return files


def table(suite, rows):
    lines = ['### %s' % suite, '']
    wls = []
    for r in rows:
        if r['workload'] not in wls:
            wls.append(r['workload'])
    for wl in wls:
        rs = [r for r in rows if r['workload'] == wl]
        lines.append('**%s**' % wl)
        lines.append('')
        if any('auc1' in r for r in rs):
            lines.append('| run | wall s | peak RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FP@E10 | FP@rec .05 | FP@rec .06 | FP@rec .07 |')
            lines.append('|---|---|---|---|---|---|---|---|---|---|')
            for r in rs:
                if 'failed_exit' in r:
                    lines.append('| %s | failed (exit %s) | %.0f | | | | | | | |' % (r['run'], r['failed_exit'], r['peak_rss_mib']))
                    continue
                f = r['fp_at_matched_recall']
                g = lambda k: '—' if f.get(k) is None else str(f.get(k))
                lines.append('| %s | %.2f | %.0f | %.4f | %d | %d | %d | %s | %s | %s |' % (
                    r['run'], r['wall_s_median'], r['peak_rss_mib'], r['auc1'], r['at_1e-3']['TP'], r['at_1e-3']['FP'],
                    r['at_10']['FP'], g('0.05'), g('0.06'), g('0.07')))
        else:
            lines.append('| run | wall s | hits E≤0.01 | E≤0.1 | E≤1 | obs/exp E≤1 |')
            lines.append('|---|---|---|---|---|---|')
            for r in rs:
                if 'null' not in r:
                    continue
                t = {x['evalue']: x for x in r['null']['thresholds']}
                lines.append('| %s | %.2f | %d | %d | %d | %.3f |' % (r['run'], r['wall_s_median'], t[0.01]['hits'], t[0.1]['hits'],
                                                                     t[1.0]['hits'], t[1.0]['observed_over_expected']))
        lines.append('')
    return '\n'.join(lines)


def main():
    results, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    md = []
    for suite in ('tune', 'eval', 'null', 'scale'):
        data = load(results, suite)
        if data is None:
            continue
        rows = [compact(r) for r in data['runs']]
        with open(os.path.join(outdir, '%s.json' % suite), 'w') as f:
            json.dump({'suite': suite, 'versions': data['versions'], 'runs': rows}, f, indent=1)
        plot(suite, data, outdir)
        md.append(table(suite, rows))
    text = '\n'.join(md)
    with open(os.path.join(outdir, 'tables.md'), 'w') as f:
        f.write(text + '\n')
    print(text)


if __name__ == '__main__':
    main()
