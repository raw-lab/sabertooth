#!/usr/bin/env python3
"""Reproducible benchmark driver.

Runs search tools on the SCOPe-derived workloads made by make_datasets.py,
evaluates every run with eval_scop.py (or the null-calibration check) and
writes machine-readable results plus a Markdown summary.

  run_suite.py --data DS --out DIR --suite SUITE --sabertooth BIN
               [--baseline BIN] [--mmseqs BIN] [--diamond BIN]
               [--threads 2] [--repeat 1] [--warmup 0] [--only RUN ...]
               [--skip-existing]

Suites
  tune     sabertooth -s presets and ablations on the tuning split
  eval     held-out split: v0.2 baseline, sabertooth (fast/default/sensitive),
           MMseqs2 (-s 5.7 / 7.5), DIAMOND (sensitive / ultra-sensitive),
           exhaustive Smith-Waterman
  null     E-value calibration on null data (shuffled queries vs the real
           database; real queries vs a shuffled database)
  scale    68 M-residue database (--data DS_SCALE): presets, 1 vs 2 threads,
           cold cache, seed caps, cgroup-limited machines (1 GiB / 512 MiB /
           256 MiB; out-of-core shards), every tool
  build    createindex time and peak memory under the same limits

Every tool searches its own prebuilt database (built once per suite and timed
separately); the v0.2 baseline reads FASTA (its database path was broken).
All searches report E <= 10 with up to 4000 hits per query so precision-recall
curves cover the same range for every tool.

Output (DIR/SUITE/):
  dbs/                       prebuilt databases + build measurements
  <run>.m8                   raw hits
  <run>.json                 command, timing (measure.py), evaluation
  summary.json / summary.md  all runs
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MEASURE = os.path.join(HERE, 'measure.py')
EVAL = os.path.join(HERE, 'eval_scop.py')

# workload -> (query fasta, query labels, database key, extra-id file or None)
WORKLOADS = {
    'tune':       ('q_tune.fa', 'q_tune.labels.tsv', 'scope40', None),
    'lc_tune':    ('q_lc_tune.fa', 'q_lc_tune.labels.tsv', 'decoy', 'db_decoy.extra.ids'),
    'frag_tune':  ('q_frag_tune.fa', 'q_frag_tune.labels.tsv', 'scope40', None),
    'multi_tune': ('q_multi_tune.fa', 'q_multi_tune.labels.tsv', 'scope40', None),
    'null_tune':  ('q_null_tune.fa', 'q_null_tune.labels.tsv', 'scope40', None),
    'eval':       ('q_eval.fa', 'q_eval.labels.tsv', 'scope40', None),
    'lc':         ('q_lc.fa', 'q_lc.labels.tsv', 'decoy', 'db_decoy.extra.ids'),
    'frag':       ('q_frag.fa', 'q_frag.labels.tsv', 'scope40', None),
    'multi':      ('q_multi.fa', 'q_multi.labels.tsv', 'scope40', None),
    'null_q':     ('q_null.fa', 'q_null.labels.tsv', 'scope40', None),
    'null_db':    ('q_eval.fa', 'q_eval.labels.tsv', 'null', None),
    # scaling: SCOPe40 + order-2 Markov filler (make_scale_db.sh); filler hits are FPs
    'scale':      ('q_tune.fa', 'q_tune.labels.tsv', 'scale', 'scale.extra.ids'),
}
DATABASES = {'scope40': 'scope40.fa', 'decoy': 'db_decoy.fa', 'null': 'db_null.fa', 'scale': 'scale_db.fa'}
NULL_WORKLOADS = {'null_tune', 'null_q', 'null_db'}

ST_PRESETS = {
    'st_s1': ['-s', '1'], 'st_s2': ['-s', '2'], 'st_s4': ['-s', '4'], 'st_s5.7': [],
    'st_s6.5': ['-s', '6.5'], 'st_s7.5': ['-s', '7.5'], 'st_s8.5': ['-s', '8.5'], 'st_s9': ['-s', '9'],
}


def suite_runs(suite):
    """(run name, tool, extra args, workload) tuples of a suite."""
    runs = []
    if suite == 'tune':
        for name, a in ST_PRESETS.items():
            runs.append((name, 'sabertooth', a, 'tune'))
        for wl in ('lc_tune', 'frag_tune', 'multi_tune', 'null_tune'):
            for name in ('st_s2', 'st_s5.7', 'st_s7.5'):
                runs.append((name, 'sabertooth', ST_PRESETS[name], wl))
        # ablations at the default operating point
        abl = {
            'st_s5.7_nomaskscore': ['--mask-scoring', '0'],
            'st_s5.7_nocbc': ['--comp-bias-corr', '0'],
            'st_s5.7_nomask': ['--mask', '0'],
            'st_s5.7_plain': ['--mask', '0', '--comp-bias-corr', '0', '--mask-scoring', '0'],
            'st_v02seeds': ['--kmer-score', '25', '--max-kmers', '256', '--min-diag-score', '15'],
            # prefilter heuristics: whole-diagonal best segment instead of X-drop,
            # one diagonal per target (v0.2) instead of four + chaining
            'st_s5.7_kadane': ['--ungapped-xdrop', '0'],
            'st_s5.7_diag1': ['--prefilter-diags', '1'],
        }
        for name, a in abl.items():
            for wl in ('tune', 'lc_tune'):
                runs.append((name, 'sabertooth', a, wl))
    elif suite == 'eval':
        tools = [
            ('v02_default', 'baseline', []),
            ('st_s2', 'sabertooth', ST_PRESETS['st_s2']),
            ('st_s5.7', 'sabertooth', ST_PRESETS['st_s5.7']),
            ('st_s7.5', 'sabertooth', ST_PRESETS['st_s7.5']),
            ('mm_s5.7', 'mmseqs', ['-s', '5.7']),
            ('mm_s7.5', 'mmseqs', ['-s', '7.5']),
            ('dm_sensitive', 'diamond', ['--sensitive']),
            ('dm_ultra', 'diamond', ['--ultra-sensitive']),
        ]
        for wl in ('eval', 'lc', 'frag', 'multi'):
            for name, tool, a in tools:
                runs.append((name, tool, a, wl))
        # sensitivity ceiling of the scoring scheme: exhaustive Smith-Waterman
        runs.append(('st_exhaustive', 'sabertooth_exh', [], 'eval'))
    elif suite == 'null':
        tools = [
            ('v02_default', 'baseline', []),
            ('st_s5.7', 'sabertooth', []),
            ('st_s7.5', 'sabertooth', ST_PRESETS['st_s7.5']),
            ('st_s5.7_plain', 'sabertooth', ['--mask', '0', '--comp-bias-corr', '0', '--mask-scoring', '0']),
            ('mm_s5.7', 'mmseqs', ['-s', '5.7']),
            ('dm_sensitive', 'diamond', ['--sensitive']),
            ('st_exhaustive', 'sabertooth_exh', []),
            ('st_exhaustive_plain', 'sabertooth_exh', ['--mask', '0', '--comp-bias-corr', '0', '--mask-scoring', '0']),
        ]
        for wl in ('null_q', 'null_db'):
            for name, tool, a in tools:
                runs.append((name, tool, a, wl))
    elif suite == 'scale':
        gib = 1 << 30
        runs += [
            ('st_s5.7', 'sabertooth', [], 'scale'),
            ('st_s2', 'sabertooth', ST_PRESETS['st_s2'], 'scale'),
            ('st_s7.5', 'sabertooth', ST_PRESETS['st_s7.5'], 'scale'),
            ('st_s5.7_t1', 'sabertooth', [], 'scale', {'threads': 1}),
            ('st_s5.7_cold', 'sabertooth', [], 'scale', {'evict_db': True}),
            # a machine with 1 GiB / 512 MiB of RAM (page cache included): the
            # index is larger than RAM, so pages are evicted and re-read
            # (the database files are evicted first so their pages are charged to
            # the limited cgroup instead of an earlier run's)
            # --memory-limit = 3/4 of the machine; when the index does not fit next
            # to the working memory it is searched as shards one at a time
            # (out-of-core); the shards are built once in a separate measured
            # step ('prepare': the same command with one query, same cgroup)
            ('st_s5.7_mem1g', 'sabertooth', ['--memory-limit', '768M'], 'scale',
             {'cgroup': gib, 'evict_db': True, 'prepare': True}),
            ('st_s5.7_mem512m', 'sabertooth', ['--memory-limit', '384M'], 'scale',
             {'cgroup': gib // 2, 'evict_db': True, 'prepare': True}),
            ('st_s5.7_mem256m', 'sabertooth', ['--memory-limit', '192M'], 'scale',
             {'cgroup': gib // 4, 'evict_db': True, 'prepare': True}),
            # the same machines without shards: on 1 GiB the 832 MiB index still
            # fits in the page cache; on 512 MiB it is paged in and out
            ('st_s5.7_mem1g_noshard', 'sabertooth', ['--memory-limit', '768M', '--shards', '1'], 'scale',
             {'cgroup': gib, 'evict_db': True}),
            ('st_s5.7_mem512m_noshard', 'sabertooth', ['--memory-limit', '384M', '--shards', '1'], 'scale',
             {'cgroup': gib // 2, 'evict_db': True}),
            ('mm_s5.7', 'mmseqs', ['-s', '5.7'], 'scale'),
            ('mm_s5.7_t1', 'mmseqs', ['-s', '5.7'], 'scale', {'threads': 1}),
            ('dm_sensitive', 'diamond', ['--sensitive'], 'scale'),
            ('dm_sensitive_t1', 'diamond', ['--sensitive'], 'scale', {'threads': 1}),
            ('mm_s5.7_mem1g', 'mmseqs', ['-s', '5.7'], 'scale', {'cgroup': gib, 'evict_db': True}),
            # MMseqs2 configured for the limit: its precomputed index cannot be
            # split ("Cannot fit databases into 767M"), so search the sequence
            # database and let the prefilter build split index tables itself
            ('mm_s5.7_mem1g_split', 'mmseqs', ['-s', '5.7', '--split-memory-limit', '768M'], 'scale',
             {'cgroup': gib, 'evict_db': True, 'no_index': True}),
            # ... and the same without a memory limit (cost of the on-the-fly index)
            ('mm_s5.7_noidx', 'mmseqs', ['-s', '5.7'], 'scale', {'no_index': True}),
            ('dm_sensitive_mem1g', 'diamond', ['--sensitive'], 'scale', {'cgroup': gib, 'evict_db': True}),
            ('v02_default', 'baseline', [], 'scale'),
            # v0.2 on a 1 GiB machine (it holds a hash-map index in anonymous memory)
            ('v02_default_mem1g', 'baseline', [], 'scale', {'cgroup': gib, 'evict_db': True}),
            # frequent-seed caps (off by default): speed vs recall
            ('st_s5.7_occ200', 'sabertooth', ['--max-seed-occ', '200'], 'scale'),
            ('st_s5.7_occ50', 'sabertooth', ['--max-seed-occ', '50'], 'scale'),
        ]
    else:
        sys.exit('unknown suite %s' % suite)
    return [r if len(r) == 5 else r + ({},) for r in runs]


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit('command failed (%d): %s\n%s' % (r.returncode, ' '.join(cmd), r.stderr[-3000:]))
    return r


def measure(json_path, label, cmd, repeat, stdout='-', evict=(), cgroup=0, must_succeed=True, warmup=0):
    mcmd = [sys.executable, MEASURE, '--json', json_path, '--repeat', str(repeat), '--label', label, '--stdout', stdout]
    if warmup:
        mcmd += ['--warmup', str(warmup)]
    if evict:
        mcmd += ['--evict-prefix'] + list(evict)
    if cgroup:
        mcmd += ['--cgroup-limit', str(int(cgroup))]
    r = subprocess.run(mcmd + ['--'] + cmd, capture_output=True, text=True)
    if not os.path.exists(json_path):
        sys.exit('measure.py failed for %s: %s' % (label, r.stderr[-2000:]))
    with open(json_path) as f:
        m = json.load(f)
    if r.returncode != 0 and must_succeed:
        sys.exit('command failed (%s): %s\n%s' % (m['runs'][-1]['exit'], ' '.join(cmd), m['runs'][-1]['stderr_tail']))
    return m


def tool_version(tool, binpath):
    try:
        if tool in ('sabertooth', 'baseline'):
            out = subprocess.run([binpath, 'version'], capture_output=True, text=True).stdout
            for line in out.splitlines():
                if 'version' in line.lower() or 'sabertooth' in line.lower():
                    return line.strip()
            return out.strip().splitlines()[0] if out.strip() else '?'
        if tool == 'mmseqs':
            return 'MMseqs2 ' + subprocess.run([binpath, 'version'], capture_output=True, text=True).stdout.strip()
        if tool == 'diamond':
            return subprocess.run([binpath, 'version'], capture_output=True, text=True).stdout.strip()
    except OSError:
        return '?'
    return '?'


def build_db(tool, bins, data, dbdir, key, threads):
    """Prebuild the database of `tool` for DATABASES[key]; returns
    (search target, build record)."""
    fasta = os.path.join(data, DATABASES[key])
    os.makedirs(dbdir, exist_ok=True)
    rec_path = os.path.join(dbdir, '%s_%s.build.json' % (tool, key))
    if tool in ('baseline', 'sabertooth_exh'):
        return fasta, None
    if tool == 'sabertooth':
        prefix = os.path.join(dbdir, 'st_%s' % key)
        if not os.path.exists(rec_path):
            for f in os.listdir(dbdir):
                if f.startswith('st_%s' % key) and not f.endswith('.json'):
                    fp = os.path.join(dbdir, f)
                    shutil.rmtree(fp) if os.path.isdir(fp) else os.remove(fp)
            m1 = measure(rec_path + '.makedb', 'makedb', [bins['sabertooth'], 'makedb', fasta, prefix], 1)
            m2 = measure(rec_path + '.index', 'createindex',
                         [bins['sabertooth'], 'createindex', prefix, '--threads', str(threads)], 1)
            idx = prefix + '.stidx'
            rec = {'tool': tool, 'db': key, 'makedb': m1, 'createindex': m2,
                   'index_bytes': os.path.getsize(idx),
                   'db_bytes': sum(os.path.getsize(os.path.join(dbdir, f)) for f in os.listdir(dbdir)
                                   if f.startswith('st_%s' % key) and not f.endswith('.json')
                                   and os.path.isfile(os.path.join(dbdir, f)))}
            with open(rec_path, 'w') as f:
                json.dump(rec, f, indent=1)
        return prefix, json.load(open(rec_path))
    if tool == 'mmseqs':
        prefix = os.path.join(dbdir, 'mm_%s' % key)
        tmp = os.path.join(dbdir, 'mm_tmp_%s' % key)
        if not os.path.exists(rec_path):
            m1 = measure(rec_path + '.createdb', 'createdb', [bins['mmseqs'], 'createdb', fasta, prefix, '-v', '1'], 1)
            m2 = measure(rec_path + '.index', 'createindex',
                         [bins['mmseqs'], 'createindex', prefix, tmp, '--threads', str(threads), '-v', '1'], 1)
            shutil.rmtree(tmp, ignore_errors=True)
            rec = {'tool': tool, 'db': key, 'createdb': m1, 'createindex': m2,
                   'db_bytes': sum(os.path.getsize(os.path.join(dbdir, f)) for f in os.listdir(dbdir)
                                   if f.startswith('mm_%s' % key) and os.path.isfile(os.path.join(dbdir, f)))}
            with open(rec_path, 'w') as f:
                json.dump(rec, f, indent=1)
        return prefix, json.load(open(rec_path))
    if tool == 'diamond':
        prefix = os.path.join(dbdir, 'dm_%s' % key)
        if not os.path.exists(rec_path):
            m1 = measure(rec_path + '.makedb', 'makedb',
                         [bins['diamond'], 'makedb', '--in', fasta, '-d', prefix, '--threads', str(threads), '--quiet'], 1)
            rec = {'tool': tool, 'db': key, 'makedb': m1, 'db_bytes': os.path.getsize(prefix + '.dmnd')}
            with open(rec_path, 'w') as f:
                json.dump(rec, f, indent=1)
        return prefix, json.load(open(rec_path))
    raise ValueError(tool)


def search_cmd(tool, bins, target, query, out_m8, args, threads, tmpdir, stats_json):
    if tool == 'sabertooth':
        return [bins['sabertooth'], 'search', query, target, '--evalue', '10', '--max-hits', '4000',
                '--threads', str(threads), '--out', out_m8, '--stats-json', stats_json] + args
    if tool == 'sabertooth_exh':
        return [bins['sabertooth'], 'easysearch', query, target, '--evalue', '10', '--max-hits', '4000',
                '--threads', str(threads), '--out', out_m8, '--stats-json', stats_json] + args
    if tool == 'baseline':
        return [bins['baseline'], 'search', query, target, '--evalue', '10', '--max-hits', '4000',
                '--threads', str(threads), '--out', out_m8] + args
    if tool == 'mmseqs':
        return [bins['mmseqs'], 'easy-search', query, target, out_m8, tmpdir, '-e', '10', '--max-seqs', '4000',
                '--threads', str(threads), '-v', '1'] + args
    if tool == 'diamond':
        return [bins['diamond'], 'blastp', '-q', query, '-d', target, '-o', out_m8, '-e', '10', '-k', '4000',
                '--threads', str(threads), '--quiet'] + args
    raise ValueError(tool)


def null_eval(m8, n_queries):
    """Observed chance hits per query at E-value thresholds (expected = x)."""
    ev = []
    with open(m8) as f:
        for line in f:
            if not line.strip() or line.startswith('query\t') or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            try:
                ev.append(float(p[10]))
            except (IndexError, ValueError):
                pass
    rows = []
    for x in (1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0):
        obs = sum(1 for e in ev if e <= x)
        rows.append({'evalue': x, 'hits': obs, 'per_query': obs / n_queries, 'expected_per_query': x,
                     'observed_over_expected': (obs / n_queries) / x})
    return {'queries': n_queries, 'thresholds': rows}


def count_fasta(path):
    with open(path) as f:
        return sum(1 for line in f if line.startswith('>'))


BUILD_RUNS = [
    # seed-index construction on the scale database under memory limits
    # (`--memory-limit` = 3/4 of the machine; partitions sized to half of it)
    ('st_createindex', [], {}),
    ('st_createindex_mem1g', ['--memory-limit', '768M'], {'cgroup': 1 << 30}),
    ('st_createindex_mem512m', ['--memory-limit', '384M'], {'cgroup': 1 << 29}),
    ('st_createindex_mem256m', ['--memory-limit', '192M'], {'cgroup': 1 << 28}),
]


def build_suite(a, bins, out):
    """`--suite build`: time and peak memory of `createindex` under limits."""
    dbdir = os.path.join(out, 'dbs')
    os.makedirs(dbdir, exist_ok=True)
    key = 'scale' if os.path.exists(os.path.join(a.data, DATABASES['scale'])) else 'scope40'
    prefix = os.path.join(dbdir, 'st_%s' % key)
    if not os.path.exists(prefix + '.index'):
        run([bins['sabertooth'], 'makedb', os.path.join(a.data, DATABASES[key]), prefix])
    version = tool_version('sabertooth', bins['sabertooth'])
    records = []
    for name, args, opt in BUILD_RUNS:
        if a.only and name not in a.only:
            p = os.path.join(out, name + '.json')
            if os.path.exists(p):
                records.append(json.load(open(p)))
            continue
        idx = os.path.join(out, 'tmp_%s.stidx' % name)
        sj = os.path.join(out, name + '.stats.json')
        cmd = [bins['sabertooth'], 'createindex', prefix, idx, '--threads', str(a.threads), '--stats-json', sj] + args
        meas = measure(os.path.join(out, name + '.measure'), name, cmd, a.repeat, evict=[prefix],
                       cgroup=opt.get('cgroup', 0), must_succeed=False)
        rec = {'run': name, 'tool': 'sabertooth', 'version': version, 'database': key, 'args': args, 'options': opt,
               'command': cmd, 'threads': a.threads, 'measure': meas}
        if meas['runs'][-1]['exit'] != 0:
            rec['failed'] = {'exit': meas['runs'][-1]['exit'], 'stderr_tail': meas['runs'][-1]['stderr_tail'][-1500:]}
        else:
            rec['createindex_stats'] = json.load(open(sj))
            rec['index_bytes'] = os.path.getsize(idx)
        for f in (idx, sj, os.path.join(out, name + '.measure')):
            if os.path.exists(f):
                os.remove(f)
        with open(os.path.join(out, name + '.json'), 'w') as f:
            json.dump(rec, f, indent=1)
        records.append(rec)
        print('%-24s build      %s' % (name, build_line(rec)), flush=True)
    lines = ['# Benchmark suite `build`', '', 'Version: `%s`; database: `%s`.' % (version, key), '',
             '| run | --memory-limit | machine (cgroup) | wall s | peak RSS MiB | cgroup peak MiB | partitions | index MiB |',
             '|---|---|---|---|---|---|---|---|']
    for r in records:
        lines.append('| %s | %s | %s | %s |' % (r['run'], ' '.join(r['args'][1:]) or 'default (1 GiB partitions)',
                                              ('%d MiB' % (r['options']['cgroup'] >> 20)) if r['options'].get('cgroup') else '—',
                                              build_line(r, md=True)))
    with open(os.path.join(out, 'summary.json'), 'w') as f:
        json.dump({'suite': 'build', 'versions': {'sabertooth': version}, 'runs': records}, f, indent=1)
    with open(os.path.join(out, 'summary.md'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('wrote %s' % os.path.join(out, 'summary.md'))


def max_cgroup(m):
    """The cgroup record with the largest peak over a measurement's repetitions."""
    cgs = [x['cgroup'] for x in m['runs'] if x.get('cgroup')]
    return max(cgs, key=lambda c: c['max_usage_bytes'] or 0) if cgs else None


def build_line(r, md=False):
    m = r['measure']
    cg = max_cgroup(m)
    cgp = '%.0f' % ((cg['max_usage_bytes'] or 0) / 1048576) if cg else '—'
    if 'failed' in r:
        return ('failed (exit %s) | %.0f | %s | | ' if md else 'failed (exit %s) rss %.0f MiB cgroup %s') % (
            r['failed']['exit'], m['peak_rss_bytes']['max'] / 1048576, cgp)
    b = r['createindex_stats']['build']
    if md:
        return '%.1f | %.0f | %s | %d | %.0f' % (m['wall_s']['median'], m['peak_rss_bytes']['max'] / 1048576, cgp,
                                                b['partitions'], r['index_bytes'] / 1048576)
    return 'wall %6.2fs rss %7.1f MiB cgroup peak %s MiB partitions %d index %.0f MiB' % (
        m['wall_s']['median'], m['peak_rss_bytes']['max'] / 1048576, cgp, b['partitions'], r['index_bytes'] / 1048576)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--suite', required=True)
    ap.add_argument('--sabertooth', required=True)
    ap.add_argument('--baseline')
    ap.add_argument('--mmseqs')
    ap.add_argument('--diamond')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--repeat', type=int, default=1)
    ap.add_argument('--warmup', type=int, default=0,
                    help='unmeasured runs before each warm-cache measurement (not for --evict / cgroup runs)')
    ap.add_argument('--only', nargs='*', default=[])
    ap.add_argument('--skip-existing', action='store_true')
    a = ap.parse_args()
    bins = {'sabertooth': os.path.abspath(a.sabertooth)}
    for t in ('baseline', 'mmseqs', 'diamond'):
        if getattr(a, t):
            bins[t] = os.path.abspath(getattr(a, t))
    out = os.path.join(a.out, a.suite)
    dbdir = os.path.join(out, 'dbs')
    os.makedirs(out, exist_ok=True)
    if a.suite == 'build':
        build_suite(a, bins, out)
        return
    versions = {t: tool_version(t, b) for t, b in bins.items()}
    versions['sabertooth_exh'] = versions['sabertooth'] + ' (easysearch: exhaustive Smith-Waterman, no prefilter)'
    records = []
    for name, tool, args, wl, opt in suite_runs(a.suite):
        if (tool if tool != 'sabertooth_exh' else 'sabertooth') not in bins:
            print('skip %s (%s not given)' % (name, tool), file=sys.stderr)
            continue
        rid = '%s__%s' % (name, wl)
        rec_path = os.path.join(out, rid + '.json')
        if a.only and name not in a.only and rid not in a.only:
            # keep the suite summary complete: reuse the existing record
            if os.path.exists(rec_path):
                records.append(json.load(open(rec_path)))
            continue
        if a.skip_existing and os.path.exists(rec_path):
            records.append(json.load(open(rec_path)))
            continue
        qfa, qlab, dbkey, extra = WORKLOADS[wl]
        target, build = build_db(tool, bins, a.data, dbdir, dbkey, a.threads)
        if tool == 'mmseqs' and opt.get('no_index'):
            # the same sequence database without the precomputed index
            target = os.path.join(dbdir, 'mmnoidx_%s' % dbkey)
            if not os.path.exists(target + '.dbtype'):
                run([bins['mmseqs'], 'createdb', os.path.join(a.data, DATABASES[dbkey]), target, '-v', '1'])
        m8 = os.path.join(out, rid + '.m8')
        stats = os.path.join(out, rid + '.stats.json')
        tmp = os.path.join(out, 'tmp_' + rid)
        threads = opt.get('threads', a.threads)
        cmd = search_cmd(tool, bins, target, os.path.join(a.data, qfa), m8, args, threads, tmp, stats)
        evict = [target] if opt.get('evict_db') else []
        prep = None
        if opt.get('prepare'):
            # build what the run reuses (out-of-core shards) with one query, under
            # the same memory limit, and record it separately
            q1 = os.path.join(out, 'prep_' + rid + '.fa')
            with open(os.path.join(a.data, qfa)) as f, open(q1, 'w') as g:
                g.write(f.read().split('>', 2)[1].join(['>', '']))
            pcmd = search_cmd(tool, bins, target, q1, os.devnull, args, threads, tmp, stats + '.prepare')
            prep = measure(rec_path + '.prepare', rid + ' prepare', pcmd, 1, evict=evict, cgroup=opt.get('cgroup', 0))
            if os.path.exists(stats + '.prepare'):
                prep['sabertooth_stats'] = json.load(open(stats + '.prepare'))
                os.remove(stats + '.prepare')
            os.remove(rec_path + '.prepare')
            os.remove(q1)
        cold = bool(evict) or bool(opt.get('cgroup'))
        meas = measure(rec_path + '.measure', rid, cmd, a.repeat, evict=evict, cgroup=opt.get('cgroup', 0),
                       must_succeed=False, warmup=0 if cold else a.warmup)
        shutil.rmtree(tmp, ignore_errors=True)
        rec = {'run': name, 'tool': tool, 'version': versions.get(tool), 'workload': wl, 'database': dbkey,
               'args': args, 'options': opt, 'command': cmd, 'threads': threads, 'measure': meas, 'db_build': build}
        if prep is not None:
            rec['prepare'] = prep
        failed = meas['runs'][-1]['exit'] != 0
        if failed:
            rec['failed'] = {'exit': meas['runs'][-1]['exit'], 'stderr_tail': meas['runs'][-1]['stderr_tail'][-1500:]}
            os.remove(rec_path + '.measure')
            with open(rec_path, 'w') as f:
                json.dump(rec, f, indent=1)
            records.append(rec)
            print('%-24s %-10s FAILED (exit %s)' % (name, wl, rec['failed']['exit']), flush=True)
            continue
        if os.path.exists(stats):
            rec['sabertooth_stats'] = json.load(open(stats))
        if wl in NULL_WORKLOADS:
            rec['null'] = null_eval(m8, count_fasta(os.path.join(a.data, qfa)))
        else:
            ej = rec_path + '.eval'
            ecmd = [sys.executable, EVAL, '--labels', os.path.join(a.data, 'scope40.labels.tsv'),
                    '--queries', os.path.join(a.data, qlab), '--hits', m8, '--name', rid, '--json', ej,
                    '--recall-levels', '0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.10']
            if extra:
                ecmd += ['--db-extra', os.path.join(a.data, extra)]
            run(ecmd)
            rec['evaluation'] = json.load(open(ej))
            os.remove(ej)
        os.remove(rec_path + '.measure')
        with open(rec_path, 'w') as f:
            json.dump(rec, f, indent=1)
        records.append(rec)
        print(summarize(rec), flush=True)
    with open(os.path.join(out, 'summary.json'), 'w') as f:
        json.dump({'suite': a.suite, 'versions': versions, 'runs': records}, f, indent=1)
    with open(os.path.join(out, 'summary.md'), 'w') as f:
        f.write(markdown(a.suite, versions, records))
    print('wrote %s' % os.path.join(out, 'summary.md'))


def metric_row(rec):
    m = rec['measure']
    if 'failed' in rec:
        return {'run': rec['run'], 'workload': rec['workload'], 'failed': rec['failed']['exit'],
                'wall_s': m['wall_s']['median'], 'rss_mib': m['peak_rss_bytes']['max'] / 1048576}
    row = {'run': rec['run'], 'workload': rec['workload'], 'wall_s': m['wall_s']['median'],
           'wall_min': m['wall_s']['min'], 'wall_max': m['wall_s']['max'], 'n': m['repeat'],
           'rss_mib': m['peak_rss_bytes']['max'] / 1048576}
    if 'evaluation' in rec:
        e = rec['evaluation']
        t = {r['evalue']: r for r in e['thresholds']}
        mr = {round(x['recall'], 3): x for x in e['matched_recall']}
        row.update({'auc1': e['mean_auc1'], 'tp_1e3': t[1e-3]['TP'], 'fp_1e3': t[1e-3]['FP'],
                    'recall_1e3': t[1e-3]['recall'], 'fdp_1e3': t[1e-3]['FDP'],
                    'tp_10': t[10.0]['TP'], 'fp_10': t[10.0]['FP'],
                    'fp_at_recall': {k: v['FP'] for k, v in mr.items()}})
    if 'null' in rec:
        row['null'] = {r['evalue']: r['observed_over_expected'] for r in rec['null']['thresholds']}
        row['null_hits'] = {r['evalue']: r['hits'] for r in rec['null']['thresholds']}
    return row


def summarize(rec):
    r = metric_row(rec)
    s = '%-24s %-10s wall %6.2fs rss %7.1f MiB' % (r['run'], r['workload'], r['wall_s'], r['rss_mib'])
    cg = max_cgroup(rec['measure'])
    if cg:
        s += ' cgroup peak %.0f MiB (limit %.0f)' % ((cg['max_usage_bytes'] or 0) / 1048576, cg['limit_bytes'] / 1048576)
    if 'auc1' in r:
        s += '  AUC1 %.4f  TP/FP@1e-3 %d/%d  FP@10 %d' % (r['auc1'], r['tp_1e3'], r['fp_1e3'], r['fp_10'])
    if 'null' in r:
        s += '  null hits E<=1e-3/1e-2/0.1/1: %s' % '/'.join(str(r['null_hits'][x]) for x in (1e-3, 1e-2, 0.1, 1.0))
    return s


def markdown(suite, versions, records):
    lines = ['# Benchmark suite `%s`' % suite, '', 'Versions:', '']
    for t, v in versions.items():
        lines.append('- %s: `%s`' % (t, v))
    lines.append('')
    rows = [metric_row(r) for r in records]
    by_wl = {}
    for r in rows:
        by_wl.setdefault(r['workload'], []).append(r)
    for wl, rs in by_wl.items():
        lines.append('## %s' % wl)
        lines.append('')
        for r in [r for r in rs if 'failed' in r]:
            lines.append('- **%s failed** (exit %s after %.1f s, peak RSS %.0f MiB)' % (r['run'], r['failed'], r['wall_s'], r['rss_mib']))
        rs = [r for r in rs if 'failed' not in r]
        if not rs:
            continue
        lines.append('')
        if 'auc1' in rs[0]:
            levels = sorted(rs[0]['fp_at_recall'])
            lines.append('| run | wall s (min–max, n) | peak RSS MiB | AUC1 | TP@1e-3 | FP@1e-3 | FDP@1e-3 | FP@E10 | '
                         + ' | '.join('FP@rec %.2f' % l for l in levels) + ' |')
            lines.append('|' + '---|' * (8 + len(levels)))
            for r in rs:
                fps = ' | '.join(str(r['fp_at_recall'].get(l)) if r['fp_at_recall'].get(l) is not None else 'n/r'
                                 for l in levels)
                lines.append('| %s | %.2f (%.2f–%.2f, %d) | %.0f | %.4f | %d | %d | %.4f | %d | %s |' % (
                    r['run'], r['wall_s'], r['wall_min'], r['wall_max'], r['n'], r['rss_mib'], r['auc1'],
                    r['tp_1e3'], r['fp_1e3'], r['fdp_1e3'], r['fp_10'], fps))
        else:
            lines.append('| run | wall s | peak RSS MiB | hits E≤1e-3 | E≤1e-2 | E≤0.1 | E≤1 | obs/exp E≤0.1 | obs/exp E≤1 |')
            lines.append('|---|---|---|---|---|---|---|---|---|')
            for r in rs:
                h = r['null_hits']
                lines.append('| %s | %.2f | %.0f | %d | %d | %d | %d | %.3f | %.3f |' % (
                    r['run'], r['wall_s'], r['rss_mib'], h[1e-3], h[1e-2], h[0.1], h[1.0], r['null'][0.1], r['null'][1.0]))
        lines.append('')
    lines.append('FP@rec x = false positives at the first E-value threshold reaching recall x (n/r = not reached).')
    lines.append('')
    return '\n'.join(lines)


if __name__ == '__main__':
    main()
