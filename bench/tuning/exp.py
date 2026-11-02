#!/usr/bin/env python3
"""Quick experiment runner: one sabertooth search + evaluation, one summary line.

(Development harness kept for reference: the paths below are those of the
machine the v0.4 experiments ran on; `bench/results/v0.4/experiments.jsonl` is
its output.)

exp.py NAME --db scale|scope40|<prefix> [--q tune|eval] [--bin BIN] [--threads 2] -- <extra sabertooth args>
Appends a JSON record to $EXP_OUT/results.jsonl (default scratchpad/exp).
"""
import json, os, subprocess, sys, time, argparse

SCR = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get('EXP_OUT', os.path.join(SCR, 'exp'))
EVAL = '/home/claude/work/sabertooth/bench/eval_scop.py'
DS = '/home/claude/bench_data/ds_scale'
DBS = {
    'scale': ('/home/claude/bench_results/scale/dbs/st_scale', DS + '/scale.extra.ids'),
    'scope40': ('/home/claude/bench_results/eval/dbs/st_scope40', None),
}

ap = argparse.ArgumentParser()
ap.add_argument('name')
ap.add_argument('--db', default='scale')
ap.add_argument('--q', default='tune')
ap.add_argument('--bin', default='/home/claude/work/sabertooth/target/release/sabertooth')
ap.add_argument('--threads', type=int, default=2)
ap.add_argument('--keep', action='store_true')
argv = sys.argv[1:]
if '--' in argv:
    k = argv.index('--')
    argv, extra = argv[:k], argv[k + 1:]
else:
    extra = []
a = ap.parse_args(argv)
os.makedirs(OUT, exist_ok=True)
if a.db in DBS:
    target, extra_ids = DBS[a.db]
else:
    target, extra_ids = a.db, None
q = '%s/q_%s.fa' % (DS, a.q)
qlab = '%s/q_%s.labels.tsv' % (DS, a.q)
tag = '%s__%s_%s' % (a.name, a.db if a.db in DBS else os.path.basename(a.db), a.q)
m8 = os.path.join(OUT, tag + '.m8')
sj = os.path.join(OUT, tag + '.stats.json')
cmd = [a.bin, 'search', q, target, '--evalue', '10', '--max-hits', '4000', '--threads', str(a.threads),
       '--out', m8, '--stats-json', sj] + extra
t0 = time.time()
r = subprocess.run(cmd, capture_output=True, text=True)
wall = time.time() - t0
if r.returncode != 0:
    sys.exit('FAILED: %s\n%s' % (' '.join(cmd), r.stderr[-2000:]))
ej = os.path.join(OUT, tag + '.eval.json')
ecmd = [sys.executable, EVAL, '--labels', DS + '/scope40.labels.tsv', '--queries', qlab, '--hits', m8,
        '--name', tag, '--json', ej, '--recall-levels', '0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.10']
if extra_ids and a.db == 'scale':
    ecmd += ['--db-extra', extra_ids]
subprocess.run(ecmd, check=True, capture_output=True)
e = json.load(open(ej))
s = json.load(open(sj))
c = s['counters']
nq = c['queries']
t = {x['evalue']: x for x in e['thresholds']}
mr = {round(x['recall'], 3): x['FP'] for x in e['matched_recall']}
rec = {
    'name': a.name, 'db': a.db, 'q': a.q, 'args': extra, 'wall': wall,
    'search_s': s['wall_s']['search'], 'rss_mib': s['memory']['peak_rss_bytes'] / 1048576,
    'user_cpu': s['memory']['user_cpu_s'],
    'cpu': {k[6:]: v for k, v in c.items() if k.startswith('cpu_s_')},
    'per_q': {'kmers': c['kmers_generated'] / nq, 'postings': c['postings_visited'] / nq,
              'targets_hit': c['targets_hit'] / nq, 'evaluated': c['targets_evaluated'] / nq,
              'xdrop': c['xdrop_extensions'] / nq, 'cand': c['candidates'] / nq, 'hits': c['hits'] / nq},
    'auc1': e['mean_auc1'], 'tp3': t[1e-3]['TP'], 'fp3': t[1e-3]['FP'], 'tp10': t[10.0]['TP'], 'fp10': t[10.0]['FP'],
    'fp_at': mr,
}
with open(os.path.join(OUT, 'results.jsonl'), 'a') as f:
    f.write(json.dumps(rec) + '\n')
if not a.keep:
    for p in (m8,):
        try:
            os.remove(p)
        except OSError:
            pass
cp = rec['cpu']
print('%-28s %-7s wall %6.2f cpu %5.1f [seed %4.1f hit %4.1f ev %4.1f gap %4.1f tb %3.1f] '
      'kmer %6.0f post %7.0f tgt %6.0f eval %6.0f cand %6.0f | AUC1 %.4f TP/FP@1e-3 %5d/%-4d TP@10 %5d FP@rec.07 %s .08 %s rss %.0f' % (
          a.name, a.db, wall, rec['user_cpu'], cp.get('seed_generation', 0), cp.get('hit_collection_sort', 0),
          cp.get('diagonal_evidence', 0), cp.get('gapped_scoring', 0), cp.get('traceback', 0),
          rec['per_q']['kmers'], rec['per_q']['postings'], rec['per_q']['targets_hit'], rec['per_q']['evaluated'],
          rec['per_q']['cand'], rec['auc1'], rec['tp3'], rec['fp3'], rec['tp10'], mr.get(0.07), mr.get(0.08), rec['rss_mib']))
