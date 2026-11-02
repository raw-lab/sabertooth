#!/usr/bin/env python3
"""Seed-index layout comparison: sparse k-mer table vs the default (`auto`).

For each database FASTA: makedb once, then one index per layout
(`createindex --index-layout sparse|auto`) in its own directory, and the same
search (tuning queries, default settings) against each: 1 unmeasured warm-up
and N measured runs. Records wall time, peak RSS (rusage of the child), index
size, the layout `indexinfo` reports and the sha256 of the hit file (the
layouts must give identical output).

  run_layouts.py --bin BIN --queries Q.fa --out OUT_DIR NAME=DB.fa [NAME=DB.fa ...]

Writes OUT_DIR/<NAME>__<layout>.json and OUT_DIR/summary.json.
"""
import argparse, hashlib, json, os, platform, shutil, subprocess, sys, time


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def run_measured(cmd):
    import resource
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, status, ru = os.wait4(p.pid, 0)
    wall = time.perf_counter() - t0
    err = p.stderr.read().decode()
    p.stdout.close()
    p.stderr.close()
    code = os.waitstatus_to_exitcode(status)
    if code != 0:
        sys.exit('FAILED (%d): %s\n%s' % (code, ' '.join(cmd), err[-2000:]))
    return {'wall_s': wall, 'peak_rss_bytes': ru.ru_maxrss * 1024, 'user_s': ru.ru_utime, 'sys_s': ru.ru_stime,
            'major_faults': ru.ru_majflt, 'stderr_tail': err[-300:]}


def host():
    model, mhz = '', ''
    try:
        for line in open('/proc/cpuinfo'):
            if line.startswith('model name') and not model:
                model = line.split(':', 1)[1].strip()
    except OSError:
        pass
    return {'cpu': model, 'cpus': os.cpu_count(), 'kernel': platform.release(),
            'mem_total_mib': os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') // 1048576}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', required=True)
    ap.add_argument('--queries', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--work', help='scratch directory for databases (default OUT/work, removed afterwards)')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--repeat', type=int, default=3)
    ap.add_argument('--layouts', default='sparse,auto')
    ap.add_argument('dbs', nargs='+')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    work = a.work or os.path.join(a.out, 'work')
    os.makedirs(work, exist_ok=True)
    vout = subprocess.run([a.bin, '--version'], capture_output=True, text=True).stdout
    version = next((l.strip() for l in vout.splitlines() if 'Sabertooth v' in l), vout.strip())
    summary = {'binary': {'path': a.bin, 'sha256': sha256(a.bin), 'version': version},
               'queries': {'fasta': a.queries, 'sha256': sha256(a.queries)}, 'threads': a.threads,
               'repeat': a.repeat, 'warmup': 1, 'host': host(), 'runs': []}
    for spec in a.dbs:
        name, fasta = spec.split('=', 1)
        hits = {}
        for layout in a.layouts.split(','):
            d = os.path.join(work, '%s_%s' % (name, layout))
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d)
            db = os.path.join(d, 'db')
            subprocess.run([a.bin, 'makedb', fasta, db], check=True, capture_output=True)
            b = run_measured([a.bin, 'createindex', db, '--index-layout', layout, '--threads', str(a.threads)])
            info = subprocess.run([a.bin, 'indexinfo', db + '.stidx'], capture_output=True, text=True).stdout
            m8 = os.path.join(d, 'hits.m8')
            cmd = [a.bin, 'search', a.queries, db, '--threads', str(a.threads), '--evalue', '10',
                   '--max-hits', '4000', '--out', m8, '--require-index']
            run_measured(cmd)  # warm-up (page cache, not recorded)
            runs = [run_measured(cmd) for _ in range(a.repeat)]
            walls = sorted(r['wall_s'] for r in runs)
            rec = {'name': name, 'layout': layout, 'fasta': fasta, 'fasta_sha256': sha256(fasta),
                   'index_bytes': os.path.getsize(db + '.stidx'), 'indexinfo': info,
                   'build': b, 'command': cmd, 'runs': runs,
                   'wall_s': {'median': walls[len(walls) // 2], 'min': walls[0], 'max': walls[-1]},
                   'peak_rss_bytes_median': sorted(r['peak_rss_bytes'] for r in runs)[len(runs) // 2],
                   'hits_sha256': sha256(m8)}
            hits[layout] = rec['hits_sha256']
            with open(os.path.join(a.out, '%s__%s.json' % (name, layout)), 'w') as f:
                json.dump(rec, f, indent=1)
            summary['runs'].append({k: rec[k] for k in ('name', 'layout', 'index_bytes', 'wall_s',
                                                         'peak_rss_bytes_median', 'hits_sha256')})
            print('%-8s %-7s index %7.1f MiB  wall %.2f (%.2f-%.2f) s  peak RSS %.0f MiB' % (
                name, layout, rec['index_bytes'] / 1048576, rec['wall_s']['median'], rec['wall_s']['min'],
                rec['wall_s']['max'], rec['peak_rss_bytes_median'] / 1048576), flush=True)
            shutil.rmtree(d, ignore_errors=True)
        same = len(set(hits.values())) == 1
        summary.setdefault('identical_hits', {})[name] = same
        print('%-8s identical hits across layouts: %s' % (name, same), flush=True)
    with open(os.path.join(a.out, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    if not a.work:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
