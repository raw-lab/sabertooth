#!/bin/bash
# Collect the v0.4 run records into bench/results/v0.4, regenerate the tables
# (bench/results/v0.4/*.md and the README blocks) and build the tarball.
#
#   package_v04.sh SERIES_OUT SUITES_OUT SERIES_DATA OUT_DIR [NAME=DIR ...]
#
# NAME=DIR copies the records (*.json, and README.md / *.sh / *.py) of DIR into
# bench/results/v0.4/NAME (layout=run_layouts.py output, ab=run_ab.py output,
# dev=development records).
set -euo pipefail
SERIES=$1; SUITES=$2; SDATA=$3; OUT=$4; shift 4
SRC=$(cd "$(dirname "$0")/.." && pwd)
R=$SRC/bench/results/v0.4
mkdir -p "$R" "$OUT"
# database-size series: every run and build record (no hit files, no databases)
rm -rf "$R/series"
mkdir -p "$R/series"
cp "$SERIES/provenance.json" "$SERIES/summary.json" "$R/series/"
cp "$SDATA/MANIFEST.tsv" "$R/series/MANIFEST.tsv"
for d in "$SERIES"/*M; do
  n=$(basename "$d")
  mkdir -p "$R/series/$n"
  cp -p "$d"/*.json "$R/series/$n/"
done
# suites (v0.4 and the v0.3 binary): run records + database build records
rm -rf "$R/suites"
for v in v0.4 v0.3; do
  for suite in eval null presets scale; do
    [ -d "$SUITES/$v/$suite" ] || continue
    mkdir -p "$R/suites/$v/$suite"
    find "$SUITES/$v/$suite" -maxdepth 1 -name '*__*.json' ! -name '*.stats.json' -exec cp -p {} "$R/suites/$v/$suite/" \;
    [ -f "$SUITES/$v/$suite/summary.md" ] && cp "$SUITES/$v/$suite/summary.md" "$R/suites/$v/$suite/"
    if [ -d "$SUITES/$v/$suite/dbs" ]; then
      mkdir -p "$R/suites/$v/$suite/dbs"
      cp -p "$SUITES/$v/$suite"/dbs/*.build.json "$R/suites/$v/$suite/dbs/" 2>/dev/null || true
    fi
  done
done
# extra record sets (index layouts, A/B runs, development records)
for kv in "$@"; do
  name=${kv%%=*}; dir=${kv#*=}
  rm -rf "${R:?}/$name"
  mkdir -p "$R/$name"
  (cd "$dir" && find . \( -name '*.json' -o -name 'README.md' -o -name '*.sh' -o -name '*.py' \) -print0 \
     | xargs -0 -I{} cp -p --parents {} "$R/$name/")
done
# which host produced each record (the container moved once during the runs;
# bench/ENVIRONMENT.md): records written before the restart are host A
python3 - "$R" <<'PY'
import os, sys, datetime
R = sys.argv[1]
cut = datetime.datetime(2026, 10, 20, 2, 42, tzinfo=datetime.timezone.utc).timestamp()
rows = {'A': [], 'B': []}
for root, _, files in os.walk(R):
    for f in sorted(files):
        rel = os.path.relpath(os.path.join(root, f), R)
        series_run = rel.startswith('series/') and rel.count('/') == 2
        suite_run = rel.startswith('suites/') and ('__' in f or f.endswith('.build.json'))
        extra_run = rel.split('/')[0] in ('layout', 'ab', 'dev')
        if f.endswith('.json') and (series_run or suite_run or extra_run):
            p = os.path.join(root, f)
            rows['A' if os.path.getmtime(p) < cut else 'B'].append(os.path.relpath(p, R))
with open(os.path.join(R, 'HOSTS.md'), 'w') as out:
    out.write('# Host of every v0.4 run and build record\n\nHost A / host B: `bench/ENVIRONMENT.md`. Records written before the container restart at 2026-10-20T02:42Z are host A, later ones host B (file modification times, preserved when the records were copied here). Listed: `series/`, `suites/`, `layout/`, `ab/` and `dev/`. Not listed: `tuning/`, `experiments.jsonl` and `variant_cache/`, which are all host A (they were written before the restart; their file times in this tree are copy times).\n')
    for h in ('A', 'B'):
        out.write('\n## Host %s (%d records)\n\n' % (h, len(rows[h])))
        for p in sorted(rows[h]):
            out.write('- `%s`\n' % p)
PY
cd "$SRC"
python3 bench/series_tables.py "$R/series" --md "$R/series.md" > /dev/null
python3 bench/readme_tables_v04.py --v04 "$R/suites/v0.4" --v03 "$R/suites/v0.3" --series "$R/series" > "$R/tables.md"
python3 bench/readme_tables_v04.py --v04 "$R/suites/v0.4" --v03 "$R/suites/v0.3" --series "$R/series" --fill README.md
rm -rf bench/__pycache__ bench/tuning/__pycache__
cd "$SRC/.."
rm -f "$OUT/sabertooth-0.4.0.tar.gz"
tar --exclude='sabertooth/target' --exclude='sabertooth/vendor/hydra-mpp-core/target' \
    --exclude='__pycache__' -czf "$OUT/sabertooth-0.4.0.tar.gz" sabertooth
ls -la "$OUT/sabertooth-0.4.0.tar.gz"
tar -tzf "$OUT/sabertooth-0.4.0.tar.gz" | wc -l
