#!/usr/bin/env python3
"""Build every benchmark dataset from an ASTRAL SCOPe FASTA (deterministic).

Inputs:  --scope scope40.fa   (headers ">domain sccs ...", one sequence line or
                              wrapped; residues any case)
Outputs (in --out):
  scope40.fa / scope40.labels.tsv      the curated database (truth = SCOP labels)
  split.tsv                            superfamily -> tune|eval (sha256 parity)
  q_eval.fa / .labels.tsv              held-out queries (eval superfamilies, >=1 TP)
  q_tune.fa / .labels.tsv              tuning queries (tune superfamilies, >=1 TP)
  q_frag.fa / .labels.tsv              fragments of eval queries (30-60% length)
  q_multi.fa / .labels.tsv             two-domain fusions from different folds
  q_lc.fa / .labels.tsv                eval queries with a low-complexity insert
  q_null.fa / .labels.tsv              shuffled eval queries (no true homologs)
  decoy_shuffled.fa                    composition-preserving shuffles of all domains
  decoy_winshuf.fa                     windowed (10 aa) shuffles: local composition kept
  decoy_lc.fa                          shuffled backgrounds + low-complexity segments
  decoy_biased.fa                      skewed-composition sequences (TM-, acidic-, basic-,
                                       Pro-, Ser/Thr-, Cys-rich)
  db_decoy.fa / db_decoy.extra.ids     scope40 + all decoys (decoy hits are FPs)
  db_null.fa                           shuffled scope40 (null database)
  markov_<N>.fa (with --markov N)      order-2 Markov synthetic proteins (scaling filler)
  MANIFEST.tsv                         sha256 + record/residue counts of every file
"""
import argparse
import hashlib
import os
import random
import sys
from collections import defaultdict

AA = 'ACDEFGHIKLMNPQRSTVWY'


def read_fasta(path):
    recs = []
    h = None
    s = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if h is not None:
                    recs.append((h, ''.join(s)))
                h = line[1:]
                s = []
            else:
                s.append(line.strip())
    if h is not None:
        recs.append((h, ''.join(s)))
    return recs


def write_fasta(path, recs):
    with open(path, 'w') as f:
        for h, s in recs:
            f.write('>%s\n%s\n' % (h, s))


def write_labels(path, rows):
    with open(path, 'w') as f:
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')


def sf(sccs):
    return '.'.join(sccs.split('.')[:3])


def fold(sccs):
    return '.'.join(sccs.split('.')[:2])


def shuffle_seq(rng, s):
    l = list(s)
    rng.shuffle(l)
    return ''.join(l)


def window_shuffle(rng, s, w=10):
    out = []
    for i in range(0, len(s), w):
        out.append(shuffle_seq(rng, s[i:i + w]))
    return ''.join(out)


LC_KINDS = [
    lambda rng, n: rng.choice('QNSEPGAD') * n,                                   # homopolymer
    lambda rng, n: (rng.choice(['SR', 'GA', 'EK', 'PS', 'GS', 'QP', 'NQ', 'AE']) * n)[:n],  # dipeptide
    lambda rng, n: (rng.choice(['GPP', 'PEV', 'SGG', 'EEK', 'PAP']) * n)[:n],     # tripeptide (collagen-like)
    lambda rng, n: ''.join(('L' if i % 7 in (0, 3) else rng.choice('EKQRA')) for i in range(n)),  # heptad
    lambda rng, n: ''.join(rng.choice('ST') if rng.random() < 0.7 else rng.choice(AA) for _ in range(n)),  # S/T-rich
    lambda rng, n: ''.join(rng.choice('QN') if rng.random() < 0.8 else rng.choice(AA) for _ in range(n)),  # Q/N-rich prion-like
]


def lc_segment(rng, n):
    return rng.choice(LC_KINDS)(rng, n)


BIASES = {
    'tm': 'LLLLIIIVVVFFAAGGMWC',
    'acidic': 'DDDDEEEEESSAGTNQ',
    'basic': 'KKKKRRRRGSPAQT',
    'pro': 'PPPPPAASSGTEQ',
    'st': 'SSSSTTTTGNAPQ',
    'cys': 'CCCCGGKRYSTN',
}


def biased_seq(rng, kind, n):
    pool = BIASES[kind]
    return ''.join(rng.choice(pool) if rng.random() < 0.75 else rng.choice(AA) for _ in range(n))


def markov_model(seqs):
    counts = defaultdict(lambda: defaultdict(int))
    first = defaultdict(int)
    for s in seqs:
        s = ''.join(c for c in s if c in AA)
        if len(s) < 3:
            continue
        first[s[:2]] += 1
        for i in range(2, len(s)):
            counts[s[i - 2:i]][s[i]] += 1
    model = {}
    for k, d in counts.items():
        tot = sum(d.values())
        acc = 0
        cum = []
        for c, v in sorted(d.items()):
            acc += v
            cum.append((acc / tot, c))
        model[k] = cum
    ftot = sum(first.values())
    acc = 0
    fcum = []
    for k, v in sorted(first.items()):
        acc += v
        fcum.append((acc / ftot, k))
    return model, fcum


def pick(cum, u):
    lo, hi = 0, len(cum) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if cum[mid][0] < u:
            lo = mid + 1
        else:
            hi = mid
    return cum[lo][1]


def markov_seq(rng, model, fcum, n):
    s = list(pick(fcum, rng.random()))
    while len(s) < n:
        k = s[-2] + s[-1]
        cum = model.get(k)
        if cum is None:
            s.append(rng.choice(AA))
        else:
            s.append(pick(cum, rng.random()))
    return ''.join(s[:n])


def manifest_row(path):
    h = hashlib.sha256()
    nrec = 0
    nres = 0
    with open(path, 'rb') as f:
        for line in f:
            h.update(line)
            if line.startswith(b'>'):
                nrec += 1
            else:
                nres += len(line.strip())
    return (os.path.basename(path), h.hexdigest(), nrec, nres)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scope', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=20261018)
    ap.add_argument('--max-eval-queries', type=int, default=4000)
    ap.add_argument('--max-tune-queries', type=int, default=1000)
    ap.add_argument('--variant-queries', type=int, default=1000)
    ap.add_argument('--markov', type=int, action='append', default=[])
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)
    out = lambda name: os.path.join(args.out, name)

    raw = read_fasta(args.scope)
    recs = []
    lab = {}
    for h, s in raw:
        p = h.split()
        dom, sccs = p[0], p[1]
        s = s.upper()
        recs.append((dom, sccs, s))
        lab[dom] = sccs
    write_fasta(out('scope40.fa'), [('%s %s' % (d, c), s) for d, c, s in recs])
    write_labels(out('scope40.labels.tsv'), [(d, c, len(s)) for d, c, s in recs])

    # superfamily split
    sf_members = defaultdict(list)
    for d, c, s in recs:
        sf_members[sf(c)].append(d)
    split = {}
    for x in sorted(sf_members):
        split[x] = 'eval' if hashlib.sha256(x.encode()).digest()[0] % 2 == 0 else 'tune'
    write_labels(out('split.tsv'), sorted(split.items()))
    seq_of = {d: s for d, c, s in recs}

    def queries_for(half, limit):
        ids = [d for d, c, s in recs if split[sf(c)] == half and len(sf_members[sf(c)]) >= 2]
        r = random.Random(args.seed + (1 if half == 'eval' else 2))
        r.shuffle(ids)
        ids = sorted(ids[:limit])
        return ids

    q_eval = queries_for('eval', args.max_eval_queries)
    q_tune = queries_for('tune', args.max_tune_queries)
    for name, ids in (('q_eval', q_eval), ('q_tune', q_tune)):
        write_fasta(out(name + '.fa'), [('%s %s' % (d, lab[d]), seq_of[d]) for d in ids])
        write_labels(out(name + '.labels.tsv'), [(d, d, len(seq_of[d])) for d in ids])

    def variants(qids, suffix, seed_off):
        vr = random.Random(args.seed + seed_off)
        var_ids = sorted(vr.sample(qids, min(args.variant_queries, len(qids))))
        # fragments
        frag, frag_lab = [], []
        for d in var_ids:
            s = seq_of[d]
            if len(s) < 60:
                continue
            L = max(30, int(len(s) * vr.uniform(0.3, 0.6)))
            st = vr.randint(0, len(s) - L)
            qid = 'frag_%s_%d_%d' % (d, st + 1, st + L)
            frag.append((qid, s[st:st + L]))
            frag_lab.append((qid, d, L))
        write_fasta(out('q_frag%s.fa' % suffix), frag)
        write_labels(out('q_frag%s.labels.tsv' % suffix), frag_lab)
        # multidomain fusions (different folds)
        multi, multi_lab = [], []
        for d in var_ids:
            for _ in range(20):
                e = vr.choice(qids)
                if fold(lab[e]) != fold(lab[d]):
                    break
            qid = 'multi_%s_%s' % (d, e)
            s = seq_of[d] + seq_of[e]
            multi.append((qid, s))
            multi_lab.append((qid, '%s,%s' % (d, e), len(s)))
        write_fasta(out('q_multi%s.fa' % suffix), multi)
        write_labels(out('q_multi%s.labels.tsv' % suffix), multi_lab)
        # low-complexity inserts
        lcq, lcq_lab = [], []
        for d in var_ids:
            s = seq_of[d]
            seg = lc_segment(vr, vr.randint(30, 60))
            p = vr.randint(0, len(s))
            qid = 'lc_%s' % d
            q = s[:p] + seg + s[p:]
            lcq.append((qid, q))
            lcq_lab.append((qid, d, len(q)))
        write_fasta(out('q_lc%s.fa' % suffix), lcq)
        write_labels(out('q_lc%s.labels.tsv' % suffix), lcq_lab)
        # null queries: shuffled queries (independent shuffle seed)
        nr = random.Random(args.seed + seed_off + 1)
        nullq = [('null_%s' % d, shuffle_seq(nr, seq_of[d])) for d in var_ids]
        write_fasta(out('q_null%s.fa' % suffix), nullq)
        write_labels(out('q_null%s.labels.tsv' % suffix), [(q, 'NULL', len(s)) for q, s in nullq])
        return var_ids

    var_ids = variants(q_eval, '', 3)
    variants(q_tune, '_tune', 13)

    # decoys
    dr = random.Random(args.seed + 5)
    shuf = [('shuf_%s' % d, shuffle_seq(dr, s)) for d, c, s in recs]
    write_fasta(out('decoy_shuffled.fa'), shuf)
    win = [('winshuf_%s' % d, window_shuffle(dr, s)) for d, c, s in recs]
    write_fasta(out('decoy_winshuf.fa'), win)
    lcd = []
    for i in range(2000):
        d, c, s = recs[dr.randrange(len(recs))]
        bg = shuffle_seq(dr, s)
        for _ in range(dr.randint(1, 3)):
            seg = lc_segment(dr, dr.randint(25, 70))
            p = dr.randint(0, len(bg))
            bg = bg[:p] + seg + bg[p:]
        lcd.append(('lcdecoy_%d' % i, bg))
    write_fasta(out('decoy_lc.fa'), lcd)
    bd = []
    kinds = sorted(BIASES)
    for i in range(2000):
        k = kinds[i % len(kinds)]
        bd.append(('biased_%s_%d' % (k, i), biased_seq(dr, k, dr.randint(80, 400))))
    write_fasta(out('decoy_biased.fa'), bd)
    db_decoy = [('%s %s' % (d, c), s) for d, c, s in recs] + win + lcd + bd
    write_fasta(out('db_decoy.fa'), db_decoy)
    with open(out('db_decoy.extra.ids'), 'w') as f:
        for h, _ in win + lcd + bd:
            f.write(h.split()[0] + '\n')
    nd = random.Random(args.seed + 6)
    write_fasta(out('db_null.fa'), [('nulldb_%s' % d, shuffle_seq(nd, s)) for d, c, s in recs])

    # scaling filler
    if args.markov:
        model, fcum = markov_model([s for _, _, s in recs])
        lens = [len(s) for _, _, s in recs]
        for n in args.markov:
            mr = random.Random(args.seed + 100 + n)
            path = out('markov_%d.fa' % n)
            with open(path, 'w') as f:
                for i in range(n):
                    L = lens[mr.randrange(len(lens))]
                    f.write('>mk%d\n%s\n' % (i, markov_seq(mr, model, fcum, L)))

    rows = []
    for name in sorted(os.listdir(args.out)):
        if name.endswith('.fa'):
            rows.append(manifest_row(out(name)))
    with open(out('MANIFEST.tsv'), 'w') as f:
        f.write('file\tsha256\trecords\tresidues\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')
    print('eval queries %d, tune queries %d, variants %d' % (len(q_eval), len(q_tune), len(var_ids)))
    for r in rows:
        print('%-24s %8d records %10d residues' % (r[0], r[2], r[3]))


if __name__ == '__main__':
    main()
