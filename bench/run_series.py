#!/usr/bin/env python3
"""Head-to-head scaling benchmark over a database-size series.

For every database size: build each tool's database/index (timed, peak memory,
size on disk), run every search configuration on the held-out queries
(warm-up runs, then repeated measured runs), evaluate against SCOP truth
(filler hits are false positives), record everything, then delete the
databases (disk) unless --keep-dbs.

  run_series.py --data DS_SCALE --series SERIES_DIR --out OUT \\
      --sabertooth BIN [--sabertooth-v03 BIN] [--mmseqs BIN] [--diamond BIN] \\
      [--sizes 3,17,34,68,134,265] [--queries q_eval1k] [--threads 2] \\
      [--repeat 3] [--warmup 1] [--only RUN ...] [--keep-dbs]

Size 3 is SCOPe40 itself (2.6 M residues); other sizes read
SERIES_DIR/series_<N>M.fa (make_scale_series.py). Output: OUT/<N>M/ with one
JSON record per run (`<run>.json`) and per build (`build_<tool>.json`), and
OUT/summary.json (every record plus provenance: binary and input hashes).
`series_tables.py OUT` renders the tables.

Search configurations (all report E <= 10, up to 4,000 hits per query, every
tool on its own prebuilt index; RUNS below):
  st_s1 .. st_s7.5           Sabertooth v0.4 presets
  st03_s5.7                  Sabertooth v0.3.0 default (its own format-1 index)
  mm_s1 .. mm_s7.5           MMseqs2 -s 1 / 4 / 5.7 / 7.5
  dm_sensitive               DIAMOND --sensitive
  *_t1                       1 thread (265 M only)
  *_mem1g                    inside a 1 GiB memory cgroup, database files evicted from
                             the page cache before each measured run (68 M only)
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MEASURE = os.path.join(HERE, 'measure.py')
EVAL = os.path.join(HERE, 'eval_scop.py')
GIB = 1 << 30

# name, tool, extra args, options (threads: override; mem: cgroup bytes; sizes: restrict)
RUNS = [
    ('st_s1', 'st', ['-s', '1'], {}),
    ('st_s2', 'st', ['-s', '2'], {}),
    ('st_s4', 'st', ['-s', '4'], {}),
    ('st_s5.7', 'st', [], {}),
    ('st_s7.5', 'st', ['-s', '7.5'], {}),
    ('st03_s5.7', 'st03', [], {}),
    ('mm_s1', 'mm', ['-s', '1'], {}),
    ('mm_s4', 'mm', ['-s', '4'], {}),
    ('mm_s5.7', 'mm', ['-s', '5.7'], {}),
    ('mm_s7.5', 'mm', ['-s', '7.5'], {}),
    ('dm_sensitive', 'dm', ['--sensitive'], {}),
    ('st_s5.7_t1', 'st', [], {'threads': 1, 'sizes': [265]}),
    ('mm_s5.7_t1', 'mm', ['-s', '5.7'], {'threads': 1, 'sizes': [265]}),
    ('st_s5.7_mem1g', 'st', ['--memory-limit', '1G'], {'mem': GIB, 'sizes': [68]}),
    ('st03_s5.7_mem1g', 'st03', ['--memory-limit', '1G'], {'mem': GIB, 'sizes': [68]}),
    ('mm_s5.7_mem1g', 'mm', ['-s', '5.7'], {'mem': GIB, 'sizes': [68]}),
    ('dm_sensitive_mem1g', 'dm', ['--sensitive'], {'mem': GIB, 'sizes': [68]}),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def sh(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('command failed (%d): %s\n%s' % (r.returncode, ' '.join(cmd), r.stderr[-3000:]))
    return r


def measure(path, label, cmd, repeat, warmup=0, mem=0, allow_fail=False, evict=None):
    mcmd = [sys.executable, MEASURE, '--json', path, '--repeat', str(repeat), '--label', label]
    if warmup:
        mcmd += ['--warmup', str(warmup)]
    if mem:
        mcmd += ['--cgroup-limit', str(mem)]
    if evict:
        # page-cache pages already cached outside the memory cgroup are not
        # charged to it: evict the database files before every measured run so
        # the limit covers them (a machine with that much RAM, cold cache)
        mcmd += ['--evict-prefix', evict]
    r = subprocess.run(mcmd + ['--'] + cmd, capture_output=True, text=True)
    if not os.path.exists(path):
        if allow_fail:
            return None
        sys.exit('measure.py failed for %s: %s' % (label, r.stderr[-3000:]))
    m = json.load(open(path))
    os.remove(path)
    if m['runs'][-1]['exit'] != 0 and not allow_fail:
        sys.exit('command failed: %s\n%s' % (' '.join(cmd), m['runs'][-1]['stderr_tail']))
    return m


def files_bytes(prefix):
    d, b = os.path.dirname(prefix), os.path.basename(prefix)
    return sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d)
               if f.startswith(b) and os.path.isfile(os.path.join(d, f)))


def remove_prefix(prefix):
    d, b = os.path.dirname(prefix), os.path.basename(prefix)
    for f in os.listdir(d):
        if f.startswith(b):
            p = os.path.join(d, f)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)


def build(tool, bins, fasta, dbdir, threads, rec_path):
    """Build `tool`'s database for `fasta`; returns (search target, record)."""
    if os.path.exists(rec_path):
        rec = json.load(open(rec_path))
        if any(f.startswith(os.path.basename(rec['target'])) for f in os.listdir(dbdir)):
            return rec['target'], rec
    os.makedirs(dbdir, exist_ok=True)
    pre = os.path.join(dbdir, tool + '_db')
    steps = []
    if tool in ('st', 'st03'):
        steps.append(('makedb', [bins[tool], 'makedb', fasta, pre], False))
        steps.append(('createindex', [bins[tool], 'createindex', pre, '--threads', str(threads)], False))
    elif tool == 'mm':
        steps.append(('createdb', [bins['mm'], 'createdb', fasta, pre, '-v', '1'], False))
        # a failed index build (e.g. not enough RAM) leaves MMseqs2 searching
        # without a precomputed index — recorded, not fatal
        steps.append(('createindex', [bins['mm'], 'createindex', pre, os.path.join(dbdir, 'mm_tmp'),
                                      '--threads', str(threads), '-v', '1'], True))
    elif tool == 'dm':
        steps.append(('makedb', [bins['dm'], 'makedb', '--in', fasta, '-d', pre, '--threads', str(threads),
                                 '--quiet'], False))
    else:
        raise ValueError(tool)
    rec = {'tool': tool, 'fasta': fasta, 'target': pre, 'steps': {}}
    for name, cmd, may_fail in steps:
        m = measure(rec_path + '.m', '%s %s' % (tool, name), cmd, 1, allow_fail=may_fail)
        rec['steps'][name] = m
        if m is None or m['runs'][-1]['exit'] != 0:
            rec['index_failed'] = True
    shutil.rmtree(os.path.join(dbdir, 'mm_tmp'), ignore_errors=True)
    rec['bytes'] = files_bytes(pre)
    ok = [s for s in rec['steps'].values() if s is not None]
    rec['wall_s'] = sum(s['wall_s']['median'] for s in ok)
    rec['peak_rss_bytes'] = max(s['peak_rss_bytes']['max'] for s in ok)
    with open(rec_path, 'w') as f:
        json.dump(rec, f, indent=1)
    return pre, rec


def search_cmd(tool, bins, target, query, out_m8, args, threads, tmp):
    if tool in ('st', 'st03'):
        return [bins[tool], 'search', query, target, '--evalue', '10', '--max-hits', '4000',
                '--threads', str(threads), '--out', out_m8] + args
    if tool == 'mm':
        return [bins['mm'], 'easy-search', query, target, out_m8, tmp, '-e', '10', '--max-seqs', '4000',
                '--threads', str(threads), '-v', '1'] + args
    if tool == 'dm':
        return [bins['dm'], 'blastp', '-q', query, '-d', target, '-o', out_m8, '-e', '10', '-k', '4000',
                '--threads', str(threads), '--quiet'] + args
    raise ValueError(tool)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--series', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--sabertooth', required=True)
    ap.add_argument('--sabertooth-v03')
    ap.add_argument('--mmseqs')
    ap.add_argument('--diamond')
    ap.add_argument('--sizes', default='3,17,34,68,134,265')
    ap.add_argument('--queries', default='q_eval1k')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--repeat', type=int, default=3)
    ap.add_argument('--warmup', type=int, default=1)
    ap.add_argument('--only', nargs='*', default=[])
    ap.add_argument('--keep-dbs', action='store_true')
    a = ap.parse_args()
    bins = {'st': os.path.abspath(a.sabertooth)}
    for k, v in (('st03', a.sabertooth_v03), ('mm', a.mmseqs), ('dm', a.diamond)):
        if v:
            bins[k] = os.path.abspath(v)
    qfa = os.path.join(a.series, a.queries + '.fa')
    qlab = os.path.join(a.series, a.queries + '.labels.tsv')
    if not os.path.exists(qfa):
        qfa = os.path.join(a.data, a.queries + '.fa')
        qlab = os.path.join(a.data, a.queries + '.labels.tsv')
    os.makedirs(a.out, exist_ok=True)
    prov_path = os.path.join(a.out, 'provenance.json')
    prov = json.load(open(prov_path)) if os.path.exists(prov_path) else {}
    for k, b in bins.items():
        r = subprocess.run([b, 'version'], capture_output=True, text=True)
        lines = [x.strip() for x in (r.stdout + '\n' + r.stderr).splitlines() if x.strip()]
        ver = next((x for x in lines if 'sabertooth v' in x.lower() or 'version' in x.lower()), lines[0] if lines else '?')
        prov.setdefault('binaries', {})[k] = {'path': b, 'sha256': sha256(b), 'version': ver}
    prov['queries'] = {'fasta': qfa, 'sha256': sha256(qfa), 'labels': qlab, 'labels_sha256': sha256(qlab)}
    prov['threads'] = a.threads
    prov['repeat'] = a.repeat
    prov['warmup'] = a.warmup
    sizes = [int(x) for x in a.sizes.split(',')]
    for size in sizes:
        if size == 3:
            fasta = os.path.join(a.data, 'scope40.fa')
            extra = None
        else:
            fasta = os.path.join(a.series, 'series_%dM.fa' % size)
            extra = os.path.join(a.series, 'series_%dM.extra.ids' % size)
        if 'db_%d' % size not in prov.get('databases', {}):
            prov.setdefault('databases', {})['db_%d' % size] = {'fasta': fasta, 'sha256': sha256(fasta)}
        with open(prov_path, 'w') as f:
            json.dump(prov, f, indent=1)
        out = os.path.join(a.out, '%dM' % size)
        dbdir = os.path.join(out, 'dbs')
        os.makedirs(dbdir, exist_ok=True)
        builds = {}
        for name, tool, args, opt in RUNS:
            if tool not in bins or (a.only and name not in a.only):
                continue
            want = opt.get('sizes')
            if want == 'largest' and size != max(sizes):
                continue
            if isinstance(want, list) and size not in want:
                continue
            rid = os.path.join(out, name)
            if os.path.exists(rid + '.json'):
                continue
            if tool not in builds:
                builds[tool] = build(tool, bins, fasta, dbdir, a.threads, os.path.join(out, 'build_%s.json' % tool))
            target, brec = builds[tool]
            threads = opt.get('threads', a.threads)
            mem = opt.get('mem', 0)
            m8 = rid + '.m8'
            tmp = os.path.join(out, 'tmp_' + name)
            cmd = search_cmd(tool, bins, target, qfa, m8, args, threads, tmp)
            # memory-limited runs: the warm-up (outside the cgroup) builds what the
            # tool reuses (v0.3's out-of-core shards); the measured runs start cold
            meas = measure(rid + '.measure', '%dM %s' % (size, name), cmd, a.repeat, a.warmup, mem, allow_fail=bool(mem),
                           evict=target if mem else None)
            shutil.rmtree(tmp, ignore_errors=True)
            rec = {'size_m': size, 'run': name, 'tool': tool, 'args': args, 'command': cmd, 'threads': threads,
                   'memory_limit_bytes': mem, 'cold_cache': bool(mem), 'measure': meas,
                   'database': {'fasta': fasta, 'build': {k: v for k, v in brec.items() if k != 'steps'},
                                'build_steps': {k: (None if v is None else {'wall_s': v['wall_s']['median'],
                                                                            'peak_rss_bytes': v['peak_rss_bytes']['max'],
                                                                            'exit': v['runs'][-1]['exit']})
                                                for k, v in brec['steps'].items()}}}
            ok = meas is not None and meas['runs'][-1]['exit'] == 0 and os.path.exists(m8)
            if ok:
                ej = rid + '.eval.json'
                ecmd = [sys.executable, EVAL, '--labels', os.path.join(a.data, 'scope40.labels.tsv'), '--queries', qlab,
                        '--hits', m8, '--name', name, '--json', ej, '--recall-levels', '0.05,0.06,0.07,0.08,0.10']
                if extra:
                    ecmd += ['--db-extra', extra]
                sh(ecmd)
                rec['evaluation'] = json.load(open(ej))
                os.remove(ej)
            if os.path.exists(m8):
                os.remove(m8)
            with open(rid + '.json', 'w') as f:
                json.dump(rec, f, indent=1)
            if ok:
                ev = rec['evaluation']
                t = {x['evalue']: x for x in ev['thresholds']}
                print('%4dM %-18s wall %7.2f s cpu %7.1f s rss %6.0f MiB  AUC1 %.4f TP/FP@1e-3 %d/%d' % (
                    size, name, meas['wall_s']['median'],
                    sum(r['user_s'] + r['sys_s'] for r in meas['runs']) / len(meas['runs']),
                    meas['peak_rss_bytes']['max'] / 1048576, ev['mean_auc1'], t[1e-3]['TP'], t[1e-3]['FP']), flush=True)
            else:
                print('%4dM %-18s FAILED (%s)' % (size, name, 'no record' if meas is None else meas['runs'][-1]['stderr_tail'][-300:]), flush=True)
        if not a.keep_dbs:
            for tool in list(builds):
                remove_prefix(builds[tool][0])
    recs = []
    for size in sizes:
        d = os.path.join(a.out, '%dM' % size)
        for name, _, _, _ in RUNS:
            p = os.path.join(d, name + '.json')
            if os.path.exists(p):
                recs.append(json.load(open(p)))
    with open(os.path.join(a.out, 'summary.json'), 'w') as f:
        json.dump({'provenance': prov, 'runs': recs}, f, indent=1)


if __name__ == '__main__':
    main()
