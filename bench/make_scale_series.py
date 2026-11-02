#!/usr/bin/env python3
"""Database-size series for the scaling comparison.

Every database = SCOPe40 (all 14,323 domains, the only true homologs) + a
prefix of an order-2 Markov filler (SCOPe composition and length
distribution; no homology to anything, so every filler hit is a false
positive). Larger databases are supersets of smaller ones.

The filler is markov_360000.fa from make_datasets.py (seed 20261018 + 100 +
360000), followed — for databases above 68 M residues — by further blocks of
360,000 sequences from the same model with seeds 20261018 + 1000 + block.

  make_scale_series.py --data DS_SCALE --out DIR --filler 80000,172000,360000,720000

writes DIR/series_<residues>M.fa (named by total residues, rounded),
DIR/series_<residues>M.extra.ids (filler ids) and DIR/MANIFEST.tsv (sha256,
records, residues). With 360,000 filler sequences the file is byte-identical
to the v0.3 scale database (scope40.fa + markov_360000.fa).
"""
import argparse
import hashlib
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_datasets as md  # noqa: E402

BLOCK = 360000


def read_fasta(path):
    out = []
    h = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if h is not None:
                    out.append((h, ''.join(seq)))
                h = line[1:]
                seq = []
            elif line:
                seq.append(line.strip())
    if h is not None:
        out.append((h, ''.join(seq)))
    return out


def filler_blocks(data, seed, model, fcum, lens):
    """Yield (id, seq) of the filler in order: markov_360000.fa, then new blocks."""
    for h, s in read_fasta(os.path.join(data, 'markov_360000.fa')):
        yield h.split()[0], s
    b = 1
    while True:
        mr = random.Random(seed + 1000 + b)
        for i in range(BLOCK):
            L = lens[mr.randrange(len(lens))]
            yield 'mkb%d_%d' % (b, i), md.markov_seq(mr, model, fcum, L)
        b += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True, help='make_datasets.py output with markov_360000.fa')
    ap.add_argument('--out', required=True)
    ap.add_argument('--filler', default='80000,172000,360000,720000',
                    help='filler sequences per database (comma-separated)')
    ap.add_argument('--seed', type=int, default=20261018)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    scope = read_fasta(os.path.join(a.data, 'scope40.fa'))
    seqs = [s for _, s in scope]
    model, fcum = md.markov_model(seqs)
    lens = [len(s) for s in seqs]
    scope_res = sum(lens)
    counts = sorted(int(x) for x in a.filler.split(','))
    rows = []
    gen = filler_blocks(a.data, a.seed, model, fcum, lens)
    filler = []  # (id, seq) generated so far
    for n in counts:
        while len(filler) < n:
            filler.append(next(gen))
        total = scope_res + sum(len(s) for _, s in filler)
        name = 'series_%dM' % round(total / 1e6)
        path = os.path.join(a.out, name + '.fa')
        with open(path, 'w') as f:
            for h, s in scope:
                f.write('>%s\n%s\n' % (h, s))
            for fid, s in filler:
                f.write('>%s\n%s\n' % (fid, s))
        with open(os.path.join(a.out, name + '.extra.ids'), 'w') as f:
            for fid, _ in filler:
                f.write(fid + '\n')
        rows.append(md.manifest_row(path))
        print('%s: %d sequences, %d residues' % (name, rows[-1][2], rows[-1][3]), flush=True)
    with open(os.path.join(a.out, 'MANIFEST.tsv'), 'w') as f:
        f.write('file\tsha256\trecords\tresidues\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')


if __name__ == '__main__':
    main()
