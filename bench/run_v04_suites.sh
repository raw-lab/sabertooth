#!/bin/bash
# The v0.4 SCOPe40 suites and the 68 M memory-limited runs, for v0.4 and the
# v0.3 "before" binary (each in its own output tree with its own databases).
#
#   run_v04_suites.sh DS_SCALE OUT ST04 ST03 MMSEQS DIAMOND
#
# DS_SCALE = make_datasets.py output that also holds scale_db.fa.
# Timed workloads (held-out `eval`, the 68 M memory limits) are medians of 3
# measured runs (after one warm-up where the cache is warm); the other
# workloads are single runs (their accuracy is deterministic).
set -euo pipefail
D=$1; OUT=$2; ST4=$3; ST3=$4; MM=$5; DM=$6
B=$(cd "$(dirname "$0")" && pwd)
O4=$OUT/v0.4; O3=$OUT/v0.3
R="python3 $B/run_suite.py --data $D --skip-existing"

# held-out queries, timed
$R --out $O4 --suite eval --sabertooth $ST4 --mmseqs $MM --diamond $DM --repeat 3 --warmup 1 \
   --only st_s2__eval st_s5.7__eval st_s7.5__eval mm_s5.7__eval mm_s7.5__eval dm_sensitive__eval dm_ultra__eval
$R --out $O3 --suite eval --sabertooth $ST3 --repeat 3 --warmup 1 --only st_s2__eval st_s5.7__eval st_s7.5__eval
# low-complexity / fragment / two-domain workloads and exhaustive Smith-Waterman: single runs
$R --out $O4 --suite eval --sabertooth $ST4 --mmseqs $MM --diamond $DM --repeat 1
$R --out $O3 --suite eval --sabertooth $ST3 --repeat 1 --only st_s2 st_s5.7 st_s7.5
# E-value calibration on null data
$R --out $O4 --suite null --sabertooth $ST4 --mmseqs $MM --diamond $DM --repeat 1
$R --out $O3 --suite null --sabertooth $ST3 --repeat 1 --only st_s5.7 st_s7.5 st_s5.7_plain
# every -s preset (tuning queries; SCOPe40 and 68 M residues)
$R --out $O4 --suite presets --sabertooth $ST4 --repeat 1
$R --out $O3 --suite presets --sabertooth $ST3 --repeat 1
# 68 M residues: out-of-core search on 512 MiB / 256 MiB machines, cold cache
$R --out $O4 --suite scale --sabertooth $ST4 --repeat 3 --only st_s5.7_mem512m st_s5.7_mem256m st_s5.7_cold
$R --out $O3 --suite scale --sabertooth $ST3 --repeat 3 --only st_s5.7_mem512m st_s5.7_mem256m st_s5.7_cold
echo SUITES_DONE
