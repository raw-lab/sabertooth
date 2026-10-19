# Benchmark suite `build`

Version: `Sabertooth v0.3.0`; database: `scale`.

| run | --memory-limit | machine (cgroup) | wall s | peak RSS MiB | cgroup peak MiB | partitions | index MiB |
|---|---|---|---|---|---|---|---|
| st_createindex | default (1 GiB partitions) | — | 13.4 | 1457 | — | 1 | 833 |
| st_createindex_mem1g | 768M | 1024 MiB | 15.5 | 685 | 1024 | 3 | 833 |
| st_createindex_mem512m | 384M | 512 MiB | 16.7 | 452 | 512 | 6 | 833 |
| st_createindex_mem256m | 192M | 256 MiB | 23.2 | 259 | 256 | 11 | 833 |
