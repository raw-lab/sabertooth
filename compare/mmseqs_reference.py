#!/usr/bin/env python3
"""
Independent reference implementation of MMseqs2's PSSM algorithm.

This is a *separate code path* from the Rust (`Sabertooth`) implementation,
written directly from the formulas in MMseqs2's PSSMCalculator.cpp / BaseMatrix.cpp:

  1. probMatrix P(a,b) = exp(lambda * S_ab) * pBack[a] * pBack[b]
  2. R[a][b] = P(a,b) / pBack[b]                       (= P(a|b))
  3. seqWeight (Henikoff 1994, position-based)
  4. matchWeight[pos][aa] = sum_k w[k]*[msa==aa], column-normalised
  5. Neff_M per column (entropy-based)
  6. pseudocounts: g = R.f ; tau = min(1, pca/(1+Neff/pcb)) ; p = (1-tau)f + tau g
  7. PSSM = round(8*log2(p/pBack)), clamped to [-128,127]

Agreement between this and Sabertooth (two independent implementations of the
same documented algorithm) is evidence the Rust port is faithful.

Usage:  python3 mmseqs_reference.py <matrix.out> <msa.fasta> [--pssm]
"""
import sys
import math

AA = "ACDEFGHIKLMNPQRSTVWY"  # 20 real amino acids, MMseqs2 order (X is index 20)
IDX = {c: i for i, c in enumerate(AA)}
N = 20
ANY_BACK = 1e-5


def parse_matrix(path):
    lines = open(path).read().splitlines()
    pback = None
    lam = None
    col_letters = None
    raw = {}
    started = False
    for line in lines:
        if line.startswith("# Background (precomputed optional):"):
            vals = line.split(":", 1)[1].split()
            pback = [float(v) for v in vals]  # 21 values incl X
            continue
        if line.startswith("# Lambda     (precomputed optional):"):
            lam = float(line.split(":", 1)[1].strip())
            continue
        if line.startswith("#") or not line.strip():
            continue
        toks = line.split()
        if not started:
            if all(len(t) == 1 and t.isalpha() for t in toks):
                col_letters = toks
                started = True
            continue
        row = toks[0]
        for j, v in enumerate(toks[1:]):
            if j < len(col_letters):
                raw[(row, col_letters[j])] = float(v)

    # xIsPositive check over full 21-letter alphabet ('X' present)
    full = AA + "X"
    x_pos = False
    for c in full:
        if raw.get(("X", c), 0.0) > 0.0 or raw.get((c, "X"), 0.0) > 0.0:
            x_pos = True
    # background rescale for non-positive X (matches MMseqs2)
    pb = {}
    for i, c in enumerate(full):
        pb[c] = pback[i]
    if not x_pos:
        pb["X"] = ANY_BACK
        for c in AA:
            pb[c] *= (1.0 - pb["X"])

    # reconstruct joint prob and R over the 20 real AAs
    P = [[0.0] * N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            s = raw[(AA[a], AA[b])]
            P[a][b] = math.exp(lam * s) * pb[AA[a]] * pb[AA[b]]
    R = [[P[a][b] / pb[AA[b]] for b in range(N)] for a in range(N)]
    pbg = [pb[AA[a]] for a in range(N)]
    return pbg, lam, P, R, raw, pb


def int_score_matrix(raw, pb, bit_factor=2.0, bias=-0.2):
    """Integer substitution scores (as Sabertooth builds them)."""
    full = AA + "X"
    out = {}
    for a in full:
        for b in full:
            P = math.exp(_lam * raw[(a, b)]) * pb[a] * pb[b] if (a, b) in raw else 0.0
            if P <= 0:
                out[(a, b)] = -1
                continue
            s = P / (pb[a] * pb[b])
            val = bit_factor * math.log2(s) + bias
            r = val - 0.5 if val < 0 else val + 0.5
            out[(a, b)] = max(-128, min(127, int(r)))
    return out


def parse_msa(path):
    seqs = []
    cur = ""
    started = False
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if started:
                seqs.append(cur)
            cur = ""
            started = True
        else:
            cur += line.strip()
    if started:
        seqs.append(cur)

    is_a3m = any(any(c.islower() for c in s) for s in seqs)
    msa = []
    if is_a3m:
        for s in seqs:
            row = []
            for c in s:
                if c.islower():
                    continue
                if c in "-.":
                    row.append(-1)  # gap
                elif c.isupper():
                    row.append(IDX.get(c, 20))
            msa.append(row)
    else:
        width = max(len(s) for s in seqs)
        first = seqs[0]
        cols = [i for i in range(width) if i < len(first) and first[i] not in "-."]
        for s in seqs:
            row = []
            for i in cols:
                c = s[i] if i < len(s) else "-"
                row.append(-1 if c in "-." else IDX.get(c.upper(), 20))
            msa.append(row)
    return msa


def seq_weights(msa, L, K):
    w = [1e-6] * K
    nres = [sum(1 for p in range(L) if msa[k][p] != -1) for k in range(K)]
    for pos in range(L):
        nl = [0] * N
        for k in range(K):
            c = msa[k][pos]
            if 0 <= c < N:
                nl[c] += 1
        distinct = sum(1 for x in nl if x > 0)
        if distinct == 0:
            continue
        for k in range(K):
            c = msa[k][pos]
            if 0 <= c < N:
                w[k] += 1.0 / (nl[c] * distinct * (nres[k] + 30.0))
    return w


def normalize(v, bg=None):
    s = sum(v)
    if s != 0:
        return [x / s for x in v]
    return list(bg) if bg else v


def match_weights(msa, w, L, K, pbg):
    freq = [[0.0] * N for _ in range(L)]
    for pos in range(L):
        for k in range(K):
            c = msa[k][pos]
            if 0 <= c < N:
                freq[pos][c] += w[k]
        freq[pos] = normalize(freq[pos], pbg)
    return freq


def neff(freq, w, msa, L, K):
    neff_hmm = 0.0
    for pos in range(L):
        ssum = 0.0
        for aa in range(N):
            f = freq[pos][aa]
            if f > 1e-10:
                ssum -= f * math.log2(f)
        neff_hmm += 2 ** ssum
    neff_hmm /= max(1, L)
    nlim = max(10.0, neff_hmm + 1.0)
    scale = math.log2((nlim - neff_hmm) / (nlim - 1.0))
    out = []
    for pos in range(L):
        wm = -1.0 / K
        for k in range(K):
            if msa[k][pos] != -1:
                wm += w[k]
        out.append(1.0 if wm < 0 else nlim - (nlim - 1.0) * (2 ** (scale * wm)))
    return out


def pssm(msa, pbg, R, pca=1.0, pcb=1.5, bias=0.0):
    K = len(msa)
    L = len(msa[0])
    w = normalize(seq_weights(msa, L, K))
    freq = match_weights(msa, w, L, K, pbg)
    nf = neff(freq, w, msa, L, K)
    prob = [[0.0] * N for _ in range(L)]
    for pos in range(L):
        tau = min(1.0, pca / (1.0 + nf[pos] / pcb))
        for a in range(N):
            g = sum(R[a][b] * freq[pos][b] for b in range(N))
            prob[pos][a] = (1 - tau) * freq[pos][a] + tau * g
    M = [[0] * N for _ in range(L)]
    for pos in range(L):
        for a in range(N):
            lp = math.log2(prob[pos][a] / pbg[a])
            val = 8.0 * lp + 8.0 * bias
            r = val - 0.5 if val < 0 else val + 0.5
            M[pos][a] = max(-128, min(127, int(r)))
    return M, nf


if __name__ == "__main__":
    mat_path = sys.argv[1]
    pbg, _lam, P, R, raw, pb = parse_matrix(mat_path)
    if len(sys.argv) > 2 and sys.argv[2] != "--intmat":
        msa = parse_msa(sys.argv[2])
        M, nf = pssm(msa, pbg, R)
        # print PSSM as "pos aa0 aa1 ... aa19"
        for pos in range(len(M)):
            print(str(pos) + " " + " ".join(str(x) for x in M[pos]))
    else:
        # print integer substitution matrix rows for A..Y
        sm = int_score_matrix(raw, pb)
        for a in AA:
            print(a + " " + " ".join(str(sm[(a, b)]) for b in AA))
