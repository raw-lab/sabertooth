#!/usr/bin/env python3
"""
Differential comparison: Sabertooth (Rust) vs the independent MMseqs2-formula
reference (mmseqs_reference.py). Diffs (a) the integer substitution matrix and
(b) the full PSSM for a given MSA, cell by cell.
"""
import subprocess
import sys
import re

BIN = sys.argv[1] if len(sys.argv) > 1 else "target/release/sabertooth"
MATRIX = "src/data/blosum62.out"
MSA = sys.argv[2] if len(sys.argv) > 2 else "demo/family.a3m"
AA = "ACDEFGHIKLMNPQRSTVWY"

import importlib.util
spec = importlib.util.spec_from_file_location("ref", "compare/mmseqs_reference.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)


def sabertooth_intmatrix():
    """Parse the integer substitution matrix from `sabertooth info`."""
    out = subprocess.run([BIN, "info"], capture_output=True, text=True).stdout
    rows = {}
    lines = out.splitlines()
    # find the score table: a header line of residue letters, then A.. rows
    start = None
    for i, ln in enumerate(lines):
        if re.match(r"^\s+A\s+C\s+D\s+E\s+F\s+G\s+H\s+I\s+K\s+L", ln):
            start = i + 1
            break
    if start is None:
        raise RuntimeError("could not find score table in info output")
    for ln in lines[start:]:
        toks = ln.split()
        if not toks or toks[0] not in AA + "X":
            continue
        rows[toks[0]] = [int(x) for x in toks[1:1 + 21]]
    return rows


def sabertooth_pssm():
    """Parse the PSSM from `sabertooth msa2profile`."""
    res = subprocess.run([BIN, "msa2profile", MSA], capture_output=True, text=True)
    out = res.stdout
    M = []
    for ln in out.splitlines():
        toks = ln.split()
        # rows look like: <pos> <res> v0..v19 <neff>
        if len(toks) >= 22 and toks[0].isdigit():
            vals = [int(x) for x in toks[2:22]]
            M.append(vals)
    return M


def diff_matrix():
    sab = sabertooth_intmatrix()
    pbg, lam, P, R, raw, pb = ref.parse_matrix(MATRIX)
    ref._lam = lam
    refm = ref.int_score_matrix(raw, pb)
    maxd = 0
    total = 0
    exact = 0
    for a in AA:
        for j, b in enumerate(AA):
            sv = sab[a][j]
            rv = refm[(a, b)]
            d = abs(sv - rv)
            maxd = max(maxd, d)
            total += 1
            if d == 0:
                exact += 1
    print(f"  integer substitution matrix (20x20):")
    print(f"    cells compared : {total}")
    print(f"    exact matches  : {exact}/{total} ({100*exact/total:.1f}%)")
    print(f"    max |diff|     : {maxd}")
    return maxd


def diff_pssm():
    sabM = sabertooth_pssm()
    pbg, lam, P, R, raw, pb = ref.parse_matrix(MATRIX)
    msa = ref.parse_msa(MSA)
    refM, nf = ref.pssm(msa, pbg, R)
    if len(sabM) != len(refM):
        print(f"    LENGTH MISMATCH: sabertooth {len(sabM)} vs ref {len(refM)}")
        return 999
    maxd = 0
    total = 0
    exact = 0
    within1 = 0
    sumd = 0
    for pos in range(len(refM)):
        for a in range(20):
            sv = sabM[pos][a]
            rv = refM[pos][a]
            d = abs(sv - rv)
            maxd = max(maxd, d)
            sumd += d
            total += 1
            if d == 0:
                exact += 1
            if d <= 1:
                within1 += 1
    print(f"  PSSM ({len(refM)} columns x 20):")
    print(f"    cells compared : {total}")
    print(f"    exact matches  : {exact}/{total} ({100*exact/total:.2f}%)")
    print(f"    within +-1     : {within1}/{total} ({100*within1/total:.2f}%)")
    print(f"    mean |diff|    : {sumd/total:.4f}")
    print(f"    max |diff|     : {maxd}")
    return maxd


if __name__ == "__main__":
    print("=" * 60)
    print("DIFFERENTIAL COMPARISON: Sabertooth (Rust) vs MMseqs2-formula reference (Python)")
    print("=" * 60)
    md1 = diff_matrix()
    print()
    md2 = diff_pssm()
    print()
    if md1 == 0 and md2 <= 1:
        print("  VERDICT: implementations agree (exact matrix; PSSM within rounding).")
        sys.exit(0)
    else:
        print(f"  VERDICT: differences found (matrix maxd={md1}, pssm maxd={md2}).")
        sys.exit(1)
