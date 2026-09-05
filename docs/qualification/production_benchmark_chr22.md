# Production benchmark: 1000 Genomes chromosome 22

## Scope and entry state

This measurement was run on 2026-09-05 against committed standalone `gsim`
0.10.0 on branch `main`, revision
`312bdd8138d1413f3dc66ccef33c98315a9ecf14`. The entry index and worktree were
clean. The committed bounded public VCF workflow passed all 40 expectations
before measurement.

The requested local HAPNEST checkout was not present at
`C:/Users/au223366/Documents/GitHub/synthetic_data`, and no alternative local
checkout was found. Julia was also unavailable on `PATH`. Consequently no
HAPNEST timing, memory, or distributional comparison was run. In particular,
this report does not infer HAPNEST performance from source code and cannot
establish whether gsim is more memory-efficient or faster than HAPNEST.

The requested HAPNEST scientific oracle remains revision
`ba52da1a63cf609306ea92540b3d130fa1efd213`. A later comparison requires an
already provisioned checkout and Julia environment. A suitable setup and run
sequence is:

```powershell
git -C <hapnest-source> checkout --detach ba52da1a63cf609306ea92540b3d130fa1efd213
julia --project=<hapnest-source> -e "using Pkg; Pkg.instantiate(); Pkg.status()"
$env:JULIA_NUM_THREADS = "1"
julia --threads 1 --project=<hapnest-source> <benchmark-driver.jl> <prepared-input> <output>
```

These commands were not run because this milestone prohibited obtaining or
altering runtimes and checkouts. The future driver must retain the frozen
samples, variants, map, phase-specific donor rule, and model parameters below.

## Environment and inputs

- Windows build 26200, x86-64
- 13th Gen Intel Core i7-1365U, 10 physical and 12 logical cores
- 16,794,288,128 bytes physical memory (15.64 GiB)
- R 4.4.1 (ucrt), x86_64-w64-mingw32
- GCC/G++ 13.2.0
- `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `OPENBLAS_NUM_THREADS=1`
- Julia: unavailable

The VCF was the official 1000 Genomes Phase 3 GRCh37 call set:

`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`

- bytes: 205,612,353
- SHA-256: `A90C16C4FF2B3196476D506AE13CB3047FAE8670163C7C932C4B0239AEF3DAF5`
- MD5: `afdf383a32fcef939ca8a468bef23759`

The genetic map was the official 1000 Genomes/HapMap Phase II GRCh37 map from:

`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110106_recombination_hotspots/HapmapII_GRCh37_RecombinationHotspots.tar.gz`

- archive SHA-256: `80F22D9E6CB0E497074ED1BC277E765FA9D8E22F21B2F66C3B10286520F6B68F`
- extracted `genetic_map_GRCh37_chr22.txt` SHA-256:
  `FE688FB1722C0F2AAEB691D813AD5D759FD033E345369E50694F5CB04B12D8A8`
- coordinate build: GRCh37
- caller-side label conversion: map `chr22` to VCF label `22`

The accompanying README says the map was lifted from HapMap Phase II build 35
to GRCh37 and repaired around assembly rearrangements. No physical-distance
approximation was used.

## Frozen configuration

- chromosome and inclusive region: `22:20000000-22200000`
- reference samples: first 500 VCF samples, in header order
- selected-sample file MD5: `e0da4f7d01783fdcc5e49187705d9760`
- retained variants: 48,088, in VCF order
- first/last IDs: `22:20000018:G:A`, `22:22199982:T:C`
- donor populations: one (`P1`), ancestry weight 1
- founders: 1,000
- pedigree extension: 1,000 children (2,000 total individuals)
- mutation age: `1e9`
- `N=500`, `Ne=10000`, `rho=0.02`
- seed: `20260905`
- output: separate phased HAP and dosage BED runs
- execution: single-threaded

One calibration used 20 references, 50 founders, and
`22:20000000-20500000`, retaining 11,357 markers. It was used only to confirm
that the main run was bounded. Its results were not pooled with the main run.
A tiny pack/unpack warm-up preceded the measured main stages.

## Timing

Each primary stage was measured once using R `proc.time()` on the same machine.
Download time is excluded. Times are seconds.

| Stage | User | System | Wall |
|---|---:|---:|---:|
| VCF scan/filter/map interpolation and HAP/BIM/FAM creation | 74.17 | 2.03 | 80.41 |
| Open prepared HAP dataset | 0.88 | 0.67 | 1.57 |
| Load packed reference chromosome | 0.92 | 5.88 | 6.88 |
| Founder event plan only | 0.11 | 0.00 | 0.13 |
| Founder plan plus packed materialization | 17.64 | 0.53 | 18.86 |
| Founder HAP output | 0.20 | 0.08 | 0.34 |
| Founder BED output, identical-seed separate run | 0.34 | 0.06 | 0.40 |
| Pedigree construction | 0.04 | 0.02 | 0.06 |
| Packed pedigree meiosis | 9.50 | 0.02 | 9.79 |
| Repeated packed HAP chromosome load | 0.37 | 2.11 | 2.49 |
| Repeated founder simulation | 7.84 | 0.27 | 8.32 |

The approximate one-time path from VCF through founder HAP output was 108.06
seconds (import, open/load, simulation, and output). A warm prepared-reference
repeat was 12.72 seconds including open/load and HAP output, or 11.15 seconds
when reusing the open dataset. Thus conversion to HAP and reuse avoided the
80.41-second import in each later run and made this observed repeat about
8.5 times shorter than the one-time path.

The initial and repeated founder stages differed substantially (18.86 versus
8.32 seconds), as did cold and warm HAP loading (6.88 versus 2.49 seconds).
Windows filesystem cache and workstation background activity were not
controlled, so those differences are not algorithmic speedup estimates. The
HAPNEST run could not be alternated with gsim, which is a further fairness
limitation.

## Import profile

A separate non-invasive native-reader scan of the same file and selection took
88.97 seconds wall (85.47 user, 1.19 system). Map interpolation over the retained
variants took 0.05 seconds wall. The scan reported:

- total VCF records: 1,103,547
- outside selected region: 1,053,205 (95.43%)
- retained biallelic phased SNPs: 48,088
- skipped indels: 1,934
- skipped multiallelic records: 296
- skipped symbolic/breakend records: 24
- other, missing, unphased, non-diploid, or duplicate-ID records: 0
- maximum parsing buffer: 85,614 bytes

These categories reconcile exactly to the total record count. The benchmark did
not add invasive per-record hooks, so gzip decompression, line parsing,
outside-region scanning, sample extraction, and packed-bit construction cannot
be separated reliably. The evidence does show that sequential compressed-file
scanning dominates one-time import and that most scanned records are outside
the chosen region. The profile-only scan ran later yet was slower than the
public import, illustrating the timing variability; it should not be subtracted
from import time to estimate writing.

Within founder generation, event planning was only 0.13 seconds. The remaining
18.73 seconds of the first combined stage includes validation, allocation,
audit construction, and packed segment application. This bounds that combined
materialization path but does not attribute all of it to bit copying.

## Memory and output sizes

Peak working set was sampled every 200 ms from the benchmark `Rscript` process.
The maximum individual-process working set was 252,051,456 bytes (240.37 MiB).
No child workload process was expected, so this is not an aggregate process-tree
measurement.

| Payload | Packed bytes | Byte-per-allele equivalent | Reduction |
|---|---:|---:|---:|
| Reference, 500 individuals, H1+H2 | 6,155,264 | 48,088,000 | 7.8125x |
| Founders, 1,000 individuals, H1+H2 | 12,310,528 | 96,176,000 | 7.8125x |
| Pedigree, 2,000 individuals, H1+H2 | 24,621,056 | 192,352,000 | 7.8125x |

The packed calculation is `2 * 8 * markers * ceiling(individuals / 64)`.
The shortfall from exactly 8x is final-word sample padding. Event-plan and
founder-audit R objects were 943,032 and 1,048,088 bytes respectively. At the
pedigree stage the founder and pedigree packed payloads total 36,931,584 bytes;
the difference to the observed process peak is about 215 MB, but the sampler
does not identify which stage set the peak, so this is not a precise overhead
allocation.

| File | Bytes |
|---|---:|
| Reference HAP | 6,155,376 |
| Reference BIM | 2,435,919 |
| Reference FAM | 13,500 |
| Founder HAP | 12,310,640 |
| Founder BED | 12,022,003 |
| Founder BIM | 2,435,919 |
| Founder FAM | 24,893 |

The reference and founder HAP files each have 112 bytes of bounded format
overhead (64-byte header plus one 48-byte chromosome-table record). No dense
production haplotype or genotype matrix was requested or observed.

## Bounded scientific checks

Checks used 256 evenly spaced markers and did not constitute a new scientific
qualification. The output retained two phases. Founder H1/H2 allele frequencies
were 0.03512109 and 0.03572266. Dosage frequencies for 0, 1, and 2 were
0.95168359, 0.02578906, and 0.02252734.

There were 13,070 founder segment records across 2,000 simulated haplotypes.
Segments per haplotype had minimum 1, median 6, mean 6.535, and maximum 19.
Sampled segment lengths had median 0.50620 and mean 1.25396; copied genetic spans
had median 0.457380 and mean 0.816126. Donor-population proportion was exactly 1
for `P1`, and mutation retention was exactly 1 under mutation age `1e9`. The
pedigree audit contained 95 crossover records. Event audits from the standalone
plan and packed materializer were exactly equal.

No cross-tool allele-frequency, LD, segment, or dosage comparison is reported
because HAPNEST was not runnable.

The benchmark script completed its count reconciliation, sample/variant-order,
event-plan/audit, packed-layout, and no-dense-genotype assertions. After the
benchmark, the bounded public VCF-to-HAP-to-founder/pedigree-to-HAP/BED workflow
again passed all 40 expectations. No complete package suite or package check was
run, as production code was unchanged.

## Decision

The primary next performance milestone should be tabix/index-aware regional VCF
access. The one-time 80.41-second import dominated the VCF-to-founder workflow,
and 95.43% of VCF records were outside the requested interval. This matters most
for repeated regional imports; users who convert a whole reference once and
reuse HAP already avoid that cost.

The secondary later option is to consolidate or batch packed founder segment
application across the native boundary before considering threading or SIMD.
Only 0.13 seconds was event planning, while the combined materialization stage
took 18.86 seconds on its first run. More precise native-side attribution should
precede selecting threading versus packed-copy optimization.

Faster BGZF decompression cannot be distinguished from parsing with the current
coarse profile. HAP memory mapping is low priority because chromosome-local
packed memory remained bounded and a warm load took 2.49 seconds. There is no
measured evidence yet for SIMD. HAPNEST must be run before making any claim that
gsim is at least as memory-efficient or faster.

## Reproduction

Run `tools/benchmark/benchmark_1000G_chr22.R` with external data and output
directories. It never downloads unless explicitly enabled and refuses to write
large benchmark output inside the repository. Configuration, stage timings,
session information, checks, and compact summaries are stored outside the
repository.
