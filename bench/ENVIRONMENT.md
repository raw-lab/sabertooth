# Benchmark / test environment

## v0.4 (every number in `bench/results/v0.4/` and in the v0.4 tables)

The v0.4 runs used **two hosts**, because the session's container moved
during the runs. **Host A**: the tuning, the whole database-size series
(including its 1 GiB-machine runs) and the SCOPe40 suites (`eval`, the v0.4
`null` suite and v0.3's three `null_q` runs). **Host B**, after a container
restart at 2026-10-20T02:42Z (the same machine type as the v0.3 release,
2.1 GHz): v0.3's three `null_db` runs (accuracy only), the `presets` suite
and the 68 M `scale` runs (512 MiB / 256 MiB machines, cold cache) for both
versions, the index-layout comparison (`layout/`) and the A/B runs (`ab/`,
including MMseqs2's memory options on a 1 GiB machine).
Every timed row of a table comes from the host its table names
(`bench/results/v0.4/HOSTS.md` lists the records per host); timings from
different hosts are never compared with each other, except once on purpose
(the same development binaries ran 1.3–1.4× faster on host B, README
Benchmarks). The v0.3-release
records in `bench/results/*.json` are kept as history only.

### Host A (tuning, series, SCOPe40 timed suites)

Recorded 2026-10-20T00:25Z. Every v0.3, MMseqs2 and DIAMOND number of a
host-A table was re-measured on it with the same binaries.

```
$ rustc -vV
rustc 1.97.0 (2d8144b78 2026-07-07)
host: x86_64-unknown-linux-gnu
LLVM version: 22.1.6

$ cargo --version
cargo 1.97.0 (c980f4866 2026-06-30)

$ uname -srm
Linux 6.18.44-fc-v139 x86_64

$ lscpu (excerpt)
CPU(s):                                  2
Model name:                              Intel(R) Xeon(R) Processor @ 2.80GHz   (family 6 model 85 stepping 7: Cascade Lake)
Thread(s) per core:                      1
Core(s) per socket:                      2
L1d cache:                               64 KiB (2 instances)
L2 cache:                                2 MiB (2 instances)
L3 cache:                                33 MiB (1 instance)
Flags: sse2 popcnt bmi2 avx2 avx512f avx512bw (among others)

$ free -m
               total        used        free      shared  buff/cache   available
Mem:            8029         710        5095          14        2492        7318
Swap:              0           0           0

$ cat /sys/kernel/mm/transparent_hugepage/enabled
always [madvise] never
```

### Host B (presets, 68 M memory limits, index layouts, A/B runs)

Recorded 2026-10-20T02:42Z (after the restart): same kernel, compiler and
container limits as host A; CPU:

```
Model name:                              Intel(R) Xeon(R) Processor @ 2.10GHz   (family 6 model 207 stepping 2)
Thread(s) per core:                      1
Core(s) per socket:                      2
L1d cache:                               96 KiB (2 instances)
L2 cache:                                4 MiB (2 instances)
L3 cache:                                260 MiB (1 instance)
Flags: sse2 avx2 avx512f avx512bw (among others)
Mem: 8029 MiB, no swap
```

### Both hosts

Build profile unchanged (opt-level=3, thin LTO, codegen-units=1,
panic=abort, no `target-cpu`); SIMD kernels are chosen at run time (AVX2 by
default; the AVX-512 kernels added in v0.4 are opt-in with
`SABERTOOTH_SIMD=avx512` because they were not faster on either host's CPU —
see `CHANGES.md`). Container: 2 vCPU KVM guest, no swap, cgroup-v1 memory
controller (used for the 1 GiB runs).

| binary | SHA-256 | `version` |
|---|---|---|
| Sabertooth 0.4.0 (this tree, `cargo build --release --features distributed` — the last build of the test cycle; the distributed backend is not used by any benchmark run, see note) | `48c9d9d45bd5a5269d7eb5ec41e164bf07a6d7a0c194d78654b56bdb4a86ecd9` | `Sabertooth v0.4.0` |
| Sabertooth 0.3.0 (delivered v0.3 build, the "before") | `549806e86771203e26f00d7c79db5502e137e2b460ec055c2200db504622f6a4` | `Sabertooth v0.3.0` |
| MMseqs2 18-8cc5c (official `mmseqs-linux-avx2.tar.gz`) | `b7ef6e0e33df5dd4fa9cf988cbd8b4988c11a3a1255d2c377f13e1bb40c157fb` | `8cc5ce367b5638c4306c2d7cfc652dd099a4643f` |
| DIAMOND 2.2.8 (official `diamond-linux64.tar.gz`) | `d3d72578fa97e87cc0846a06d2f9d7dc481a2106265f9546bc91512972fcd1f1` | `diamond version 2.2.8` |

**Note on the v0.4 build.** The benchmarked binary is the `--features
distributed` build (the last build of the test cycle, copied before any
benchmark; no source file changed afterwards). A plain `cargo build
--release` of the same tree (SHA-256 `e74800dc02c60142b3643d0b559308c2ac66a20a1051772e846da3766fc4d86b`
on host B) gives byte-identical output on the tuning queries against SCOPe40
(default and `-s 7.5`) and against the 68 M-residue database (12,104 / 16,343 /
11,115 output lines, header included; re-checked by the independent review).
The crates.io packaging changes made after the review (license, crate
metadata, a version requirement on `hydra-mpp-core`, README) do not change
either binary: both rebuild to the hashes above (2026-11-02). Rust release
builds embed no timestamps, but a rebuild elsewhere may still differ in bytes
(paths, toolchain); the outputs are the reproducible quantity.

## v0.3 (history: the records in `bench/results/*.json`)


Recorded 2026-10-18T16:23Z on the machine used for every v0.3 number in bench/results.

```
$ rustc -vV
rustc 1.97.0 (2d8144b78 2026-07-07)
binary: rustc
commit-hash: 2d8144b7880597b6e6d3dfd63a9a9efae3f533d3
commit-date: 2026-07-07
host: x86_64-unknown-linux-gnu
release: 1.97.0
LLVM version: 22.1.6

$ cargo --version
cargo 1.97.0 (c980f4866 2026-06-30)

$ uname -srm
Linux 6.18.44-fc-v139 x86_64

$ lscpu (excerpt)
CPU(s):                                  2
Model name:                              Intel(R) Xeon(R) Processor @ 2.10GHz
Thread(s) per core:                      1
Core(s) per socket:                      2
Socket(s):                               1
Flags: (x86-64 with SSE2, AVX2, AVX-512 present; baseline uses SSE2 only)
L1d cache:                               96 KiB (2 instances)
L2 cache:                                4 MiB (2 instances)
L3 cache:                                260 MiB (1 instance)

$ free -m
               total        used        free      shared  buff/cache   available
Mem:            8029         544        6256          14        1469        7485
Swap:              0           0           0

$ df -h (work dir)
/dev/vda        252G   15G   29G  34% /
```

Build profile (Cargo.toml [profile.release]): opt-level=3, lto=thin, codegen-units=1, panic=abort;
no target-cpu=native (portable x86-64 baseline; SSE2 kernels need no runtime detection).

Features tested: default (no features) and `--features distributed` (vendored HydraMPP).

Container: 2 vCPU KVM guest, no swap, cgroup-v1 memory controller writable (used for memory-limit runs).

Reference tools (official prebuilt release binaries, downloaded 2026-10-18):

| tool | archive | binary SHA-256 | `version` |
|---|---|---|---|
| MMseqs2 | `mmseqs-linux-avx2.tar.gz` (GitHub release "18-8cc5c") | `b7ef6e0e33df5dd4fa9cf988cbd8b4988c11a3a1255d2c377f13e1bb40c157fb` | `8cc5ce367b5638c4306c2d7cfc652dd099a4643f` |
| DIAMOND | `diamond-linux64.tar.gz` (GitHub release v2.2.8) | `d3d72578fa97e87cc0846a06d2f9d7dc481a2106265f9546bc91512972fcd1f1` | `diamond version 2.2.8` |

The v0.2.0 baseline is the untouched 0.2.0 source tree built with the same
compiler and release profile (`sabertooth version` → `Sabertooth v0.2.0`).

The delivered source differs from the second-pass build only outside the
default search path: help text, `--help` per command, a stderr note for
temporary indexes, `.stidx` files accepted as general sequence input, and an
index no longer being adopted when its seed or masking differs from the
command's (no benchmark run uses such an index).

Benchmarked Sabertooth binaries (default features, `cargo build --release`):

| build | SHA-256 |
|---|---|
| v0.2.0 baseline | `a9fe5dbfea035bc5231e405e76e22c36454dd1ed8db28e0f7d90f619d5f67f13` |
| v0.3.0 delivered (this tree): reproduced eleven benchmark outputs of the two builds below byte for byte (eval, lc, frag, multi, tune, lc_tune `nomask`, tune `plain`, null_db, scale, scale out-of-core, scale `-s 7.5`) | `549806e86771203e26f00d7c79db5502e137e2b460ec055c2200db504622f6a4` |
| v0.3.0 interim (help text, temporary-index note): produced the re-run `tune` `plain` ablation | `5f20fd020b55b075709c8bff2bad55be022cb5ef3da77dbf02d30ebbc6a03ae9` |
| v0.3.0 second pass: every `scale` and `build` record, and the other Sabertooth runs of `tune` / `eval` / `null` | `8b3f4e9e00917bdbe4fa78d9603f4e551da3512b7a8fbf43b4bb14e4fca1af81` |
| v0.3.0 first pass (the `tune` / `eval` / `null` records of MMseqs2, DIAMOND and v0.2 were measured alongside it; its own outputs were byte-identical to the second pass) | `ab2c9381ab65f9709a6a65ee519e165edf460754d9905238fdd3446b3fc884a9` |

Note on the build profile line above: v0.3 adds AVX2 kernels selected by
**runtime** detection (`is_x86_feature_detected!`), still without
`target-cpu=native`; SSE2 remains the portable fallback.
