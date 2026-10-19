# Benchmark / test environment

Recorded 2026-10-18T16:23Z on the machine used for every number in bench/results.

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
