#!/usr/bin/env python3
"""Run a command and record wall time, peak RSS, CPU time and page faults.

  measure.py --json out.json [--repeat N] [--warmup N] [--evict FILE ...] [--label X]
             [--cgroup-limit BYTES] -- cmd args...

Each repetition runs in a fresh child Python process so that
RUSAGE_CHILDREN reflects only that run (peak RSS of the largest descendant —
MMseqs2's workflow modules run sequentially, so this is the peak of the
heaviest module). --evict drops the given files from the page cache before each
run (posix_fadvise DONTNEED) to measure cold starts. --cgroup-limit runs each
repetition inside a fresh cgroup-v1 memory cgroup (a child of this process's
memory cgroup) whose limit counts anonymous memory *and* page cache — i.e. a
machine with that much RAM; the cgroup's peak usage and failcnt are recorded.
--warmup N first runs the command N times unmeasured (no eviction, no cgroup)
so that warm-cache measurements do not depend on what ran before. Output JSON
contains every repetition plus min/median/max of wall time and peak RSS.
"""
import argparse
import json
import os
import subprocess
import sys
import statistics

CHILD = r'''
import os, resource, subprocess, sys, time, json
cmd = json.loads(sys.argv[1])
stdout = sys.argv[2]
cg = sys.argv[3] if len(sys.argv) > 3 else ''
if cg:
    with open(os.path.join(cg, 'cgroup.procs'), 'w') as f:
        f.write(str(os.getpid()))
t = time.perf_counter()
with open(stdout, 'wb') if stdout != '-' else open(os.devnull, 'wb') as so:
    p = subprocess.run(cmd, stdout=so, stderr=subprocess.PIPE)
wall = time.perf_counter() - t
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
print(json.dumps({"wall_s": wall, "exit": p.returncode, "peak_rss_bytes": ru.ru_maxrss * 1024,
                  "user_s": ru.ru_utime, "sys_s": ru.ru_stime, "minor_faults": ru.ru_minflt,
                  "major_faults": ru.ru_majflt, "stderr_tail": p.stderr.decode(errors="replace")[-2000:]}))
'''


def evict(paths):
    for p in paths:
        for root, dirs, files in (os.walk(p) if os.path.isdir(p) else [(os.path.dirname(p) or '.', [], [os.path.basename(p)])]):
            for f in files:
                fp = os.path.join(root, f)
                try:
                    fd = os.open(fp, os.O_RDONLY)
                    os.fsync(fd) if False else None
                    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                    os.close(fd)
                except OSError:
                    pass


def memory_cgroup_root():
    for line in open('/proc/self/cgroup'):
        p = line.rstrip('\n').split(':')
        if len(p) == 3 and 'memory' in p[1].split(','):
            return '/sys/fs/cgroup/memory' + p[2]
    raise SystemExit('no cgroup-v1 memory controller for this process')


def make_cgroup(limit):
    path = os.path.join(memory_cgroup_root(), 'measure_%d_%d' % (os.getpid(), int(limit)))
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, 'memory.limit_in_bytes'), 'w') as f:
        f.write(str(int(limit)))
    return path


def read_int(path):
    try:
        with open(path) as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', required=True)
    ap.add_argument('--repeat', type=int, default=1)
    ap.add_argument('--evict', nargs='*', default=[])
    ap.add_argument('--evict-prefix', nargs='*', default=[], help='evict every file starting with these prefixes')
    ap.add_argument('--stdout', default='-', help='file for the command stdout')
    ap.add_argument('--label', default='')
    ap.add_argument('--cgroup-limit', type=float, default=0, help='memory cgroup limit in bytes (0 = none)')
    ap.add_argument('--warmup', type=int, default=0, help='unmeasured runs before the measured ones')
    ap.add_argument('cmd', nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == '--' else a.cmd
    warm = []
    for i in range(a.warmup):
        r = subprocess.run([sys.executable, '-I', '-c', CHILD, json.dumps(cmd), a.stdout, ''], capture_output=True)
        if r.returncode != 0:
            print(r.stderr.decode(), file=sys.stderr)
            sys.exit(1)
        res = json.loads(r.stdout.decode())
        if res['exit'] != 0:
            print('warm-up failed (%d): %s\n%s' % (res['exit'], ' '.join(cmd), res['stderr_tail']), file=sys.stderr)
            sys.exit(res['exit'])
        warm.append(res['wall_s'])
    runs = []
    for i in range(a.repeat):
        targets = list(a.evict)
        for pre in a.evict_prefix:
            d = os.path.dirname(pre) or '.'
            b = os.path.basename(pre)
            targets += [os.path.join(d, f) for f in os.listdir(d) if f.startswith(b)]
        if targets:
            evict(targets)
        cg = make_cgroup(a.cgroup_limit) if a.cgroup_limit > 0 else ''
        r = subprocess.run([sys.executable, '-I', '-c', CHILD, json.dumps(cmd), a.stdout, cg], capture_output=True)
        if r.returncode != 0:
            print(r.stderr.decode(), file=sys.stderr)
            sys.exit(1)
        res = json.loads(r.stdout.decode())
        res['rep'] = i
        if cg:
            res['cgroup'] = {'limit_bytes': int(a.cgroup_limit),
                             'max_usage_bytes': read_int(os.path.join(cg, 'memory.max_usage_in_bytes')),
                             'failcnt': read_int(os.path.join(cg, 'memory.failcnt'))}
            try:
                os.rmdir(cg)
            except OSError:
                pass
        runs.append(res)
        if res['exit'] != 0:
            print('command failed (%d): %s\n%s' % (res['exit'], ' '.join(cmd), res['stderr_tail']), file=sys.stderr)
            break
    walls = [r['wall_s'] for r in runs]
    rss = [r['peak_rss_bytes'] for r in runs]
    out = {
        'label': a.label,
        'cmd': cmd,
        'repeat': len(runs),
        'warmup_wall_s': warm,
        'wall_s': {'min': min(walls), 'median': statistics.median(walls), 'max': max(walls),
                   'stdev': statistics.stdev(walls) if len(walls) > 1 else 0.0},
        'peak_rss_bytes': {'min': min(rss), 'median': statistics.median(rss), 'max': max(rss)},
        'runs': runs,
    }
    with open(a.json, 'w') as f:
        json.dump(out, f, indent=1)
    print('%s: wall median %.2fs (min %.2f max %.2f, n=%d) peak RSS %.1f MiB, exit %s' % (
        a.label or ' '.join(cmd[:3]), out['wall_s']['median'], out['wall_s']['min'], out['wall_s']['max'],
        len(runs), out['peak_rss_bytes']['max'] / 1048576, runs[-1]['exit']))
    if runs[-1]['exit'] != 0:
        sys.exit(runs[-1]['exit'])


if __name__ == '__main__':
    main()
