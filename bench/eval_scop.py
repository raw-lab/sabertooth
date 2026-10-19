#!/usr/bin/env python3
"""Evaluate homology-search results against SCOP(e) truth.

Truth (MMseqs2 / foldseek SCOP benchmark convention):
  * true positive  (TP): target in the query's superfamily (self-hits excluded)
  * false positive (FP): target in a different fold, or any synthetic decoy
  * ignored            : same fold, different superfamily (ambiguous)

Query records may be composite (multidomain / fragment / low-complexity
variants): their labels come from a query-label file mapping query id ->
one or more source domains; a target is TP if it shares a superfamily with ANY
source domain and FP only if its fold differs from ALL source folds.

Reports, with explicit denominators:
  * global precision / recall / FPR / FDP over a sweep of E-value thresholds
    (a precision-recall curve), recall = TP / all TP pairs,
    FPR = FP / all FP pairs, FDP = FP / (TP + FP)
  * false positives at matched recall levels (interpolated E-value threshold)
  * per-query sensitivity up to the first FP (AUC1), mean over queries
  * recall split into family-level and superfamily-only ("distant") TPs
  * stratification by query length

Usage:
  eval_scop.py --labels scope40.labels.tsv --queries queries.labels.tsv \
               --db-extra decoys.ids --hits run.m8 [--name X] [--json out.json]
"""
import argparse
import json
import math
import sys
from collections import defaultdict


def read_labels(path):
    lab = {}
    length = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 2:
                continue
            lab[p[0]] = p[1]
            if len(p) > 2:
                length[p[0]] = int(p[2])
    return lab, length


def sf(sccs):
    return '.'.join(sccs.split('.')[:3])


def fold(sccs):
    return '.'.join(sccs.split('.')[:2])


def family(sccs):
    return sccs


def read_query_labels(path, lab):
    """query id -> (list of source sccs, query length)."""
    q = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            srcs = p[1].split(',')
            q[p[0]] = ([lab[s] if s in lab else 'NULL.x.x.x' for s in srcs], int(p[2]) if len(p) > 2 else 0, srcs)
    return q


def read_hits(path, evalue_col=10, bits_col=11):
    """query -> list of (evalue, bits, target), best first (as the tool ranked
    them, re-sorted by E-value then bits for tie robustness)."""
    hits = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('query\t') or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            try:
                e = float(p[evalue_col])
                b = float(p[bits_col])
            except (IndexError, ValueError):
                continue
            hits[p[0]].append((e, -b, p[1]))
    for k in hits:
        hits[k].sort()
        # keep the best hit per (query, target)
        seen = set()
        out = []
        for e, nb, t in hits[k]:
            if t in seen:
                continue
            seen.add(t)
            out.append((e, -nb, t))
        hits[k] = out
    return hits


WINSHUF_LAB = {}


def classify(qsrc, tlabel, tid=None):
    """'TP', 'FP' or 'IGN' for a target label against query source labels."""
    if tlabel is None:
        # windowed (10 aa) shuffles keep residual similarity to their source:
        # ignore those derived from the query's own superfamily
        if tid is not None and tid.startswith('winshuf_'):
            src = WINSHUF_LAB.get(tid[8:])
            if src is not None and any(sf(s) == sf(src) for s in qsrc):
                return 'IGN'
        return 'FP'  # synthetic decoy
    if any(sf(s) == sf(tlabel) for s in qsrc):
        return 'TP'
    if all(fold(s) != fold(tlabel) for s in qsrc):
        return 'FP'
    return 'IGN'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--labels', required=True, help='domain<TAB>sccs<TAB>len for the SCOP part of the DB')
    ap.add_argument('--queries', required=True, help='query<TAB>src1,src2<TAB>len')
    ap.add_argument('--db-extra', default=None, help='ids of synthetic decoys in the DB (count toward FP denominators)')
    ap.add_argument('--hits', required=True)
    ap.add_argument('--name', default='run')
    ap.add_argument('--json', default=None)
    ap.add_argument('--recall-levels', default='0.10,0.20,0.30,0.40,0.50,0.60')
    args = ap.parse_args()

    lab, _ = read_labels(args.labels)
    WINSHUF_LAB.update(lab)
    queries = read_query_labels(args.queries, lab)
    extra = []
    if args.db_extra:
        with open(args.db_extra) as f:
            extra = [l.strip() for l in f if l.strip()]
    db_ids = list(lab.keys()) + extra
    n_db = len(db_ids)

    # per-superfamily / per-fold / per-family counts in the DB
    sf_count = defaultdict(int)
    fold_count = defaultdict(int)
    fam_count = defaultdict(int)
    for d, s in lab.items():
        sf_count[sf(s)] += 1
        fold_count[fold(s)] += 1
        fam_count[family(s)] += 1

    # denominators
    tp_total = 0
    fp_total = 0
    tp_fam_total = 0
    q_tp = {}
    for q, (srcs, qlen, src_ids) in queries.items():
        sfs = {sf(s) for s in srcs}
        n_tp = sum(sf_count[x] for x in sfs)
        # remove self-hits for real-domain sources present in the DB
        n_tp -= sum(1 for d in src_ids if d in lab)
        folds = {fold(s) for s in srcs}
        # targets in the same fold but another superfamily are ignored
        n_same_fold = sum(fold_count[f] for f in folds)
        n_fp = n_db - n_same_fold
        # family-level TPs (same family as some source), excluding self
        fams = {family(s) for s in srcs}
        n_fam = sum(fam_count[x] for x in fams) - sum(1 for d in src_ids if d in lab)
        q_tp[q] = n_tp
        tp_total += n_tp
        fp_total += n_fp
        tp_fam_total += n_fam

    hits = read_hits(args.hits)
    labelled = []  # (evalue, is_tp, is_fp, is_family, query)
    auc1 = []
    unknown_queries = 0
    for q, (srcs, qlen, src_ids) in queries.items():
        hl = hits.get(q, [])
        fams = {family(s) for s in srcs}
        tp_seen = 0
        first_fp = False
        for e, b, t in hl:
            if t in src_ids:
                continue  # self hit
            c = classify(srcs, lab.get(t), t)
            is_fam = c == 'TP' and lab.get(t) in fams
            labelled.append((e, c == 'TP', c == 'FP', is_fam, q))
            if c == 'TP' and not first_fp:
                tp_seen += 1
            if c == 'FP':
                first_fp = True
        if q_tp[q] > 0:
            auc1.append(tp_seen / q_tp[q])
    for q in hits:
        if q not in queries:
            unknown_queries += 1

    labelled.sort(key=lambda x: x[0])
    # PR sweep at every distinct E-value
    curve = []
    tp = fp = tpf = 0
    i = 0
    n = len(labelled)
    while i < n:
        e = labelled[i][0]
        while i < n and labelled[i][0] == e:
            _, a, b, f, _ = labelled[i]
            tp += a
            fp += b
            tpf += f
            i += 1
        curve.append((e, tp, fp, tpf))

    def at_threshold(x):
        best = (0, 0, 0)
        for e, a, b, f in curve:
            if e <= x:
                best = (a, b, f)
            else:
                break
        return best

    thresholds = [1e-10, 1e-5, 1e-3, 1e-2, 0.1, 1.0, 10.0]
    rows = []
    for x in thresholds:
        a, b, f = at_threshold(x)
        rows.append({
            'evalue': x,
            'TP': a, 'FP': b,
            'recall': a / tp_total if tp_total else 0.0,
            'recall_family': f / tp_fam_total if tp_fam_total else 0.0,
            'recall_superfamily_only': (a - f) / (tp_total - tp_fam_total) if tp_total > tp_fam_total else 0.0,
            'precision': a / (a + b) if a + b else 1.0,
            'FDP': b / (a + b) if a + b else 0.0,
            'FPR': b / fp_total if fp_total else 0.0,
            'FP_per_query': b / max(1, len(queries)),
        })

    # FP at matched recall: smallest FP count reaching each recall level
    matched = []
    for lvl in [float(x) for x in args.recall_levels.split(',')]:
        need = lvl * tp_total
        hit = None
        for e, a, b, f in curve:
            if a >= need:
                hit = (e, a, b)
                break
        matched.append({
            'recall': lvl,
            'reached': hit is not None,
            'evalue_threshold': hit[0] if hit else None,
            'FP': hit[2] if hit else None,
            'precision': (hit[1] / (hit[1] + hit[2])) if hit else None,
        })

    # downsampled PR curve for plotting (<= 200 points)
    step = max(1, len(curve) // 200)
    pr = [{'evalue': e, 'recall': a / tp_total if tp_total else 0, 'precision': a / (a + b) if a + b else 1.0,
           'TP': a, 'FP': b} for k, (e, a, b, f) in enumerate(curve) if k % step == 0 or k == len(curve) - 1]

    # by query length
    bins = [(0, 50), (50, 100), (100, 200), (200, 400), (400, 100000)]
    by_len = []
    for lo, hi in bins:
        qs = [q for q, (s, l, _) in queries.items() if lo <= l < hi and q_tp[q] > 0]
        if not qs:
            continue
        a = [x for x in auc1_by(queries, hits, lab, q_tp, qs)]
        by_len.append({'len_lo': lo, 'len_hi': hi, 'queries': len(qs), 'mean_auc1': sum(a) / len(a)})

    out = {
        'name': args.name,
        'queries': len(queries),
        'queries_with_tp': sum(1 for q in q_tp if q_tp[q] > 0),
        'db_size': n_db,
        'denominators': {'TP_pairs': tp_total, 'FP_pairs': fp_total, 'TP_family_pairs': tp_fam_total},
        'mean_auc1': sum(auc1) / len(auc1) if auc1 else 0.0,
        'thresholds': rows,
        'matched_recall': matched,
        'by_query_length': by_len,
        'pr_curve': pr,
        'hits_total': n,
        'hit_queries_not_in_label_file': unknown_queries,
    }
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(out, f, indent=1)
    # human summary
    print('%s: %d queries (%d with TPs), DB %d; TP pairs %d, FP pairs %d' % (
        args.name, out['queries'], out['queries_with_tp'], n_db, tp_total, fp_total))
    print('  mean AUC1 (sensitivity to 1st FP): %.4f' % out['mean_auc1'])
    print('  %-8s %8s %8s %8s %8s %9s %10s' % ('E<=', 'TP', 'FP', 'recall', 'prec', 'FDP', 'FPR'))
    for r in rows:
        print('  %-8g %8d %8d %8.4f %8.4f %9.5f %10.3e' % (r['evalue'], r['TP'], r['FP'], r['recall'], r['precision'], r['FDP'], r['FPR']))
    print('  FP at matched recall: ' + ', '.join(
        '%.2f:%s' % (m['recall'], m['FP'] if m['reached'] else 'n/a') for m in matched))


def auc1_by(queries, hits, lab, q_tp, qs):
    for q in qs:
        srcs, qlen, src_ids = queries[q]
        tp_seen = 0
        for e, b, t in hits.get(q, []):
            if t in src_ids:
                continue
            c = classify(srcs, lab.get(t), t)
            if c == 'FP':
                break
            if c == 'TP':
                tp_seen += 1
        yield tp_seen / q_tp[q]


if __name__ == '__main__':
    main()
