# Parallel founder scaling diagnosis

## Decision

The committed 0.12.0 founder implementation is healthy (decision class A). On
this controlled workload, packed materialization sped up 1.60x, 2.50x, and
3.43x with 2, 4, and 8 threads. Workers were active and well balanced, setup
and join costs were negligible, output writes were timed separately, and static
inspection found neither synchronization in the hot loop nor overlapping output
writes. No production code or package version changed.

The flattening beyond four threads is most consistent with the cost of shared,
irregular donor-panel reads and marker-major output read/modify/write traffic
meeting cache/memory and heterogeneous-core frequency limits. Reliable hardware
cache-miss, bandwidth, and core-frequency counters were not available on this
Windows installation, so the measurements do not distinguish those effects
individually. This is deliberately a cautious locality/platform diagnosis, not
a measured DRAM-bandwidth claim.

## Entry and environment

The measurement was made on 2026-09-05 from branch `main`, commit
`24c8921c40ba5db250c3db6dffec8a78ea2f1ad0`, package version 0.12.0. The
worktree and index were clean at entry. The four public workflow functions
`gsim_import_vcf()`, `gsim_reference()`, `gsim_simulate_founders()`, and
`gsim_simulate_pedigree()` were present. A small committed public founder smoke
workflow passed; its 4-reference, 8-marker, 8-founder HAP had MD5
`26b1877ac5def6ed8e4f3dadcea9f030`.

The test machine was Windows 11 x64 build 26200 with a 13th Gen Intel Core
i7-1365U (10 physical, 12 logical, hybrid cores), 16,794,288,128 bytes physical
memory, R 4.4.1, and GCC/G++ 13.2.0 using C++17. The processor's registry
nominal frequency was 2,688 MHz. BLAS, OpenMP, and MKL thread environment
settings were fixed at one; only gsim's explicit founder workers varied.

The previously prepared 1000 Genomes HAP panel was not present at the available
external paths. Following the qualification fallback, the run used a
deterministic external HAP fixture with 500 reference individuals and 49,000
markers on chromosome `22`. Its phases use varying allele-frequency thresholds
and deliberately asymmetric deterministic patterns. Genetic position runs
linearly from 0 to 5 cM. The reference HAP was 6,272,112 bytes, with MD5
`c8b31f900813da89f70b9e77a9d31fd4`; its two packed biological planes occupy
6,272,000 bytes. This fixture is realistic in dimensions and allele variation,
but it is not the 1000 Genomes LD structure. All fixture and output files lived
outside the repository and were removed after qualification.

## Frozen workload and method

The primary workload generated 10,000 founders across 49,000 markers, with one
ancestry of weight 1, `N=500`, `Ne=10000`, `rho=0.02`, mutation age `1e9`, seed
`20260905`, batch size 4,096, and phased HAP output. One warm-up of 512 founders
at eight threads preceded measurement. VCF import was excluded.

Each 1/2/4/8-thread invocation opened and loaded the same prepared HAP, planned
the same events, materialized the complete workload three times, and wrote one
output. The table reports the median materialization time and the independently
accounted single-pass workflow time. The materialization ranges expose run
noise. R process user and system CPU times were recorded with `proc.time()`.
An external PowerShell sampler measured peak process working set every 100 ms.
Filesystem cache effects were not controlled; HAP load and output timings must
therefore be read as bounded workflow observations, not storage benchmarks.

Every primary output was byte-identical: HAP MD5
`6ed6cc03d1524334d7556d0b366aa835`, BIM MD5
`4e4206cb4b90985b5b5591ce96e3cf3c`, FAM MD5
`97d8672dfb0a4223e93fbcdef22461e2`, and HAP size 123,088,112 bytes.

## Thread scaling

| Threads | Materialization median, s (range) | Material CPU, % | Speed-up | Efficiency | Material throughput, founder-markers/s | Single-pass total, s | Total speed-up | Peak working set, bytes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 40.98 (37.76--42.16) | 97.5 | 1.00 | 100.0% | 11,957,052 | 48.01 | 1.00 | 173,428,736 |
| 2 | 25.69 (25.62--25.75) | 193.3 | 1.60 | 79.8% | 19,073,569 | 32.43 | 1.48 | 176,377,856 |
| 4 | 16.42 (16.09--16.63) | 375.6 | 2.50 | 62.4% | 29,841,657 | 23.40 | 2.05 | 186,638,336 |
| 8 | 11.96 (11.55--12.24) | 578.3 | 3.43 | 42.8% | 40,969,900 | 18.99 | 2.53 | 194,322,432 |

The 8-thread process consumed about 5.8 effective cores during materialization,
despite seven spawned workers plus the calling thread. This is consistent with
the hybrid processor and a shared-data locality ceiling. It is not consistent
with a single large serialized critical section: utilization grew through all
thread counts and materialization wall time continued to fall.

Single-pass stage times were:

| Stage, seconds | 1 thread | 2 threads | 4 threads | 8 threads |
|---|---:|---:|---:|---:|
| HAP open | 0.56 | 0.53 | 0.51 | 0.53 |
| HAP chromosome load | 3.80 | 2.94 | 3.09 | 3.28 |
| Founder event planning | 1.52 | 1.30 | 1.44 | 1.48 |
| Packed materialization | 40.98 | 25.69 | 16.42 | 11.96 |
| HAP initialization | 0.29 | 0.42 | 0.22 | 0.25 |
| Positional batch writes | 3.51 | 4.59 | 4.95 | 4.54 |
| HAP finalization | 0.46 | 0.67 | 0.65 | 0.65 |

These independently timed boundaries show that output I/O did not contaminate
the materialization number. Planning, initialization, writing, and finalization
are serial in this experiment and explain why total speed-up is lower than
materialization speed-up.

## Batch sensitivity and memory

One complete workflow was run at batch sizes 1,024 and 8,192 in addition to the
primary 4,096 case. All batch/thread outputs had the same HAP/BIM/FAM hashes as
the primary run.

| Batch | 1-thread materialization, s | 8-thread materialization, s | Speed-up | 8-thread total, s | Reference + batch payload, bytes | 8-thread peak working set, bytes |
|---:|---:|---:|---:|---:|---:|---:|
| 1,024 | 26.49 | 7.08 | 3.74 | 20.60 | 18,816,000 | 158,298,112 |
| 4,096 | 40.98 | 11.96 | 3.43 | 18.99 | 56,448,000 | 194,322,432 |
| 8,192 | 44.32 | 12.69 | 3.49 | 18.51 | 106,624,000 | 238,219,264 |

The smaller batch's much lower materialization time is strong locality evidence,
although its additional positional HAP writes make its total eight-thread time
slower. Batch 8,192 improved total time only 0.48 seconds (2.5%) relative to
4,096 while adding about 50 MB of packed biological payload. Batch 4,096 remains
a sound memory/latency default.

The primary maximum event-plan object was 3,682,632 bytes. No complete founder
population, dense haplotypes, genotype matrix, per-thread reference copy, or
per-thread output batch was allocated. With fixed batch size, working biological
payload retains the contract `O((R+B) M_c / 8 + S_B)` for two phases; total
founder count affects the on-disk output and number of sequential batches, not
the largest packed batch.

## Native diagnostics

Static inspection of `materialize_founders()` established this execution model:

- one native materialization call per founder batch;
- `worker_count - 1` new `std::thread` objects per call, with the caller doing
  worker zero, followed by one join for each spawned worker;
- deterministic static partitioning into contiguous, disjoint 64-sample packed
  output words;
- distinct H1 and H2 destination planes and no overlapping worker writes;
- a single immutable shared reference panel;
- no mutex, atomic counter, R API call, shared RNG, or heap allocation in the
  interval-copy hot loop.

The three primary batches consequently created 0, 3, 9, and 21 native threads
at 1, 2, 4, and 8 requested threads. Aggregated minimal-work probes estimated
setup/join cost per batch at approximately 0.97, 1.61, and 1.94 ms for 2, 4,
and 8 threads. Thus all eight-thread setup/join activity was under 6 ms against
11.96 seconds of materialization.

Analytical worker-balance coefficients of variation across the whole primary
workload were:

| Threads | Output words CV | Segment operations CV | Marker visits CV |
|---:|---:|---:|---:|
| 2 | 0.0090 | 0.0017 | 0.0023 |
| 4 | 0.0127 | 0.0084 | 0.0032 |
| 8 | 0.0264 | 0.0299 | 0.0260 |

The small residual variation is dominated by the final partial batch. A separate
10 ms operating-system thread sample on one 4,096-founder, eight-thread batch
observed three materialization repetitions. In each repetition, the seven
spawned workers lived for approximately 5.05--5.29 seconds; the last worker
finished only 0.10--0.15 seconds after the first (under 3% of call duration).
This bounds join-tail idle time and rules out severe worker imbalance.

False sharing at partition boundaries cannot be proven absent because allocator
base alignment was not instrumented. It is unlikely to be dominant: partitions
never share a 64-bit output word, and the 1,024 batch gives each of eight workers
only two adjacent words per marker--the most boundary-sensitive case--yet it had
the fastest absolute materialization and the best eight-thread speed-up.

The dominant operations remain phase-specific reads from segment-selected donor
rows and marker-major output updates. Donor changes make the shared reference
access irregular, while each output word is repeatedly updated across different
individual segments. The strong batch-size response and diminishing effective
CPU utilization therefore support a cache/memory-locality limit. All-core turbo
reduction and scheduling onto the i7-1365U's heterogeneous cores may contribute.
Windows Performance Recorder/xperf or equivalent reliable cache/bandwidth
counters were not available, and no profiler or administrator access was added.

## Correctness and recommendation

The controlled benchmark proved exact HAP, BIM, and FAM bytes across every
thread and batch setting. The package's focused separated-simulation parity test
additionally checks phase-asymmetric donors, segment audits, batch divisions,
1/2/4/8-thread equality, chromosome identity, R RNG isolation, pedigree output,
and Mendelian consistency. Those results are recorded with the final validation
for this milestone.

No targeted production correction was justified. Thread creation, joining,
synchronization, load balance, false sharing, partitioning, serial work in the
materialization timer, and output contamination were each either measured small
or contradicted by the observations. Version 0.12.0 is retained.

For a typical chromosome job, four threads provide the better efficiency and
eight provide the shortest measured latency. Use four when running concurrent
chromosomes or other jobs; use eight for a single latency-sensitive chromosome
when memory bandwidth and CPU capacity are otherwise idle. Prefer multiple
independent approximately four-thread chromosome jobs over adding more than
four threads to one chromosome, subject to aggregate memory and storage-I/O
limits. This is an operational recommendation, not a request for process-level
sharding inside gsim.

The next smallest experiment, only if more speed is required, is a hardware-
counter run on a platform with available cache-miss, memory-bandwidth, and
per-core frequency telemetry. Do not change packed copying or add SIMD before
that evidence exists.

## Final validation

Because production code did not change, validation used the measurement-only
path authorized for this milestone rather than a separate duplicate full-suite
run. `test-separated-simulation.R` passed all 26 focused expectations. It covers
the public VCF-to-reference-to-founders-to-pedigree workflow, exact HAP bytes
and phases across batch sizes and 1/2/4/8 threads, ordered segment audits, R RNG
isolation, different-seed divergence, exact BED dosage, and zero Mendelian
inconsistencies.

A fresh isolated `R CMD INSTALL --preclean` installed and loaded gsim 0.12.0.
A built-source-tarball `R CMD check --no-manual` then completed with status
`OK`; that check also ran the package test suite once. `git diff --check`
passed. Compiled objects, DLLs, check/install trees, benchmark outputs, and the
external synthetic reference were removed after validation.

## Reproduction

`tools/benchmark/benchmark_parallel_founders.R` takes an external HAP prefix,
an external output directory, and a thread count. Its environment variables are
`GSIM_PARALLEL_FOUNDERS`, `GSIM_PARALLEL_BATCH_SIZE`,
`GSIM_PARALLEL_MATERIALIZATION_REPS`, `GSIM_PARALLEL_PROBE_REPS`,
`GSIM_PARALLEL_SEED`, and `GSIM_PARALLEL_KEEP_OUTPUT`. It never downloads data
and rejects an output directory inside the repository. With output retention
disabled, it hashes and then removes each generated HAP/BIM/FAM triplet.
