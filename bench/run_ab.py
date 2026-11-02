#!/usr/bin/env python3
"""Alternating A/B runs of search variants on one database (v0.4 ablations).

Each variant is a (binary, extra arguments, environment) triple. Every variant
runs once unmeasured (warm page cache), then the variants run in alternating
rounds (A B C, A B C, ...) so that drift of the machine affects all of them
alike. Each run is measured by measure.py (wall time, rusage) and, for
Sabertooth, records its --stats-json (CPU time per stage). Hits of each
variant are kept once and evaluated with eval_scop.py when labels are given;
variants that must give identical output are checked by SHA-256.

  run_ab.py --set NAME --out OUT --queries Q.fa --target DB [--labels Q.labels.tsv]
            [--db-extra EXTRA.ids] [--bin BIN] [--dev-pf BIN --dev-nopf BIN --dev-target IDX]
            [--rounds 3] [--threads 2]

Sets (VARIANT_SETS below):
  simd      default kernels (AVX2) vs SABERTOOTH_SIMD=avx512 — identical output
  evidence  default windowed segment evidence (--diag-window 48) vs the whole
            diagonal (--diag-window 0) vs v0.3's X-drop extension (--ungapped-xdrop 30)
  prefetch  two development binaries that differ only in the hit-loop posting
            prefetch (PF_POSTINGS = 8 vs 0 in src/prefilter.rs), as in the v0.4
            development A/B (identical output)

Writes OUT/<set>/<variant>.<round>.json (measure record + stats) and
OUT/<set>/summary.json.
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MEASURE = os.path.join(HERE, 'measure.py')
EVAL = os.path.join(HERE, 'eval_scop.py')

# name -> list of (variant, binary key, extra args, env)
VARIANT_SETS = {
    'simd': [
        ('avx2', 'bin', [], {'SABERTOOTH_SIMD': 'avx2'}),
        ('avx512', 'bin', [], {'SABERTOOTH_SIMD': 'avx512'}),
    ],
    'evidence': [
        ('window48', 'bin', [], {}),
        ('whole_diagonal', 'bin', ['--diag-window', '0'], {}),
        ('xdrop30', 'bin', ['--ungapped-xdrop', '30'], {}),
    ],
    'prefetch': [
        # the development A/B used explicit settings equal to the default at 68 M
        ('prefetch', 'dev_pf', ['--spaced', '111001101', '--kmer-score', '22', '--max-kmers', '512',
                                '--min-diag-score', '41'], {}),
        ('no_prefetch', 'dev_nopf', ['--spaced', '111001101', '--kmer-score', '22', '--max-kmers', '512',
                                     '--min-diag-score', '41'], {}),
    ],
}
IDENTICAL = {'simd', 'prefetch'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def measured(cmd, env, path, label):
    e = dict(os.environ)
    e.update(env)
    r = subprocess.run([sys.executable, MEASURE, '--json', path, '--repeat', '1', '--label', label, '--'] + cmd,
                       env=e, capture_output=True, text=True)
    if not os.path.exists(path):
        sys.exit('measure.py failed for %s: %s' % (label, r.stderr[-2000:]))
    m = json.load(open(path))
    if m['runs'][-1]['exit'] != 0:
        sys.exit('run failed: %s\n%s' % (' '.join(cmd), m['runs'][-1]['stderr_tail']))
    return m


def host():
    model = ''
    for line in open('/proc/cpuinfo'):
        if line.startswith('model name'):
            model = line.split(':', 1)[1].strip()
            break
    return {'cpu': model, 'cpus': os.cpu_count(), 'kernel': platform.release()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True, choices=sorted(VARIANT_SETS))
    ap.add_argument('--out', required=True)
    ap.add_argument('--queries', required=True)
    ap.add_argument('--target', help='database prefix with its default .stidx (bin variants)')
    ap.add_argument('--labels', help='query labels (eval_scop.py --queries)')
    ap.add_argument('--scop-labels', help='SCOPe40 labels (eval_scop.py --labels)')
    ap.add_argument('--db-extra', help='identifiers of filler sequences (false positives)')
    ap.add_argument('--bin')
    ap.add_argument('--dev-pf')
    ap.add_argument('--dev-nopf')
    ap.add_argument('--dev-target', help='index used by the development binaries')
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--threads', type=int, default=2)
    a = ap.parse_args()
    bins = {'bin': a.bin, 'dev_pf': a.dev_pf, 'dev_nopf': a.dev_nopf}
    out = os.path.join(a.out, a.set)
    os.makedirs(out, exist_ok=True)
    variants = VARIANT_SETS[a.set]

    def cmd(v, m8, stats):
        name, bkey, extra, env = v
        target = a.dev_target if bkey.startswith('dev') else a.target
        return [bins[bkey], 'search', a.queries, target, '--evalue', '10', '--max-hits', '4000',
                '--threads', str(a.threads), '--out', m8, '--stats-json', stats] + extra

    summary = {'set': a.set, 'queries': {'fasta': a.queries, 'sha256': sha256(a.queries)},
               'target': a.dev_target if a.set == 'prefetch' else a.target, 'rounds': a.rounds,
               'threads': a.threads, 'host': host(), 'binaries': {}, 'variants': {}}
    for k, b in bins.items():
        if b and any(v[1] == k for v in variants):
            summary['binaries'][k] = {'path': b, 'sha256': sha256(b)}
    hits = {}
    # warm-up: every variant once, unmeasured; its hits are the variant's output
    for v in variants:
        m8 = os.path.join(out, v[0] + '.m8')
        stats = os.path.join(out, v[0] + '.warmup.stats.json')
        measured(cmd(v, m8, stats), v[3], os.path.join(out, 'warmup.tmp.json'), v[0] + ' warm-up')
        os.remove(os.path.join(out, 'warmup.tmp.json'))
        os.remove(stats)
        hits[v[0]] = sha256(m8)
    runs = {v[0]: [] for v in variants}
    for rnd in range(1, a.rounds + 1):
        for v in variants:
            stats = os.path.join(out, '%s.%d.stats.json' % (v[0], rnd))
            m8 = os.path.join(out, '%s.%d.m8' % (v[0], rnd))
            path = os.path.join(out, '%s.%d.json' % (v[0], rnd))
            m = measured(cmd(v, m8, stats), v[3], path + '.tmp', '%s round %d' % (v[0], rnd))
            os.remove(path + '.tmp')
            same = sha256(m8) == hits[v[0]]
            os.remove(m8)
            s = json.load(open(stats))
            os.remove(stats)
            rec = {'variant': v[0], 'round': rnd, 'command': cmd(v, '<hits>', '<stats>'), 'env': v[3],
                   'measure': m, 'stats': s, 'same_hits_as_warmup': same}
            with open(path, 'w') as f:
                json.dump(rec, f, indent=1)
            c = s.get('counters', {})
            runs[v[0]].append({'wall_s': m['runs'][0]['wall_s'], 'peak_rss_bytes': m['runs'][0]['peak_rss_bytes'],
                               'cpu_s': {k[6:]: c[k] for k in c if k.startswith('cpu_s_')}})
            print('%-8s %-15s round %d  wall %.2f s  %s' % (a.set, v[0], rnd, m['runs'][0]['wall_s'],
                  ' '.join('%s %.2f' % (k[6:], c[k]) for k in sorted(c) if k.startswith('cpu_s_'))), flush=True)
    for v in variants:
        name = v[0]
        rs = runs[name]
        walls = sorted(r['wall_s'] for r in rs)
        stages = sorted(rs[0]['cpu_s'])
        summ = {'command': cmd(v, '<hits>', '<stats>'), 'env': v[3], 'hits_sha256': hits[name],
                'wall_s': {'median': walls[len(walls) // 2], 'min': walls[0], 'max': walls[-1]},
                'peak_rss_mib_max': max(r['peak_rss_bytes'] for r in rs) / 1048576,
                'cpu_s': {k: {'min': min(r['cpu_s'][k] for r in rs), 'max': max(r['cpu_s'][k] for r in rs)} for k in stages}}
        m8 = os.path.join(out, name + '.m8')
        if a.labels and a.scop_labels:
            ej = os.path.join(out, name + '.eval.json')
            ecmd = [sys.executable, EVAL, '--labels', a.scop_labels, '--queries', a.labels, '--hits', m8,
                    '--name', name, '--json', ej, '--recall-levels', '0.04,0.05,0.06,0.07,0.08,0.10']
            if a.db_extra:
                ecmd += ['--db-extra', a.db_extra]
            subprocess.run(ecmd, check=True, capture_output=True)
            e = json.load(open(ej))
            t = {x['evalue']: x for x in e['thresholds']}
            summ['auc1'] = e['mean_auc1']
            summ['tp_fp_1e-3'] = [t[1e-3]['TP'], t[1e-3]['FP']]
            summ['fp_at_recall'] = {str(x['recall']): x['FP'] for x in e['matched_recall']}
        os.remove(m8)
        summary['variants'][name] = summ
    if a.set in IDENTICAL:
        summary['identical_hits'] = len(set(hits.values())) == 1
    with open(os.path.join(out, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    for name, s in summary['variants'].items():
        print('%-8s %-15s wall %.2f (%.2f–%.2f) s  rss %.0f MiB  %s' % (
            a.set, name, s['wall_s']['median'], s['wall_s']['min'], s['wall_s']['max'], s['peak_rss_mib_max'],
            ('AUC1 %.4f TP/FP %s FP@rec %s' % (s['auc1'], s['tp_fp_1e-3'], s['fp_at_recall'])) if 'auc1' in s else ''))
    if a.set in IDENTICAL:
        print('identical hits across variants:', summary['identical_hits'])


if __name__ == '__main__':
    main()
