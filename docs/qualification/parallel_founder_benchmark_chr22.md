# Batched and parallel founder benchmark

## Contract and environment

This bounded measurement used standalone `gsim` 0.11.0 on Windows build 26200,
R 4.4.1, GCC/G++ 13.2.0, and an Intel Core i7-1365U with 10 physical/12 logical
cores and 16,794,288,128 bytes RAM. BLAS/OpenMP environment thread counts were
fixed at one. The gsim worker count was explicitly one or eight.

The prepared reference is the official 1000 Genomes Phase 3 GRCh37 chromosome
22 VCF, converted once to HAP using the official HapMap Phase II GRCh37 map.
The input checksums and map provenance are recorded in
`production_benchmark_chr22.md`. VCF scanning and conversion are excluded here.

The frozen founder configuration was:

- chromosome `22`, region 20,000,000--22,200,000;
- 500 reference samples and 48,088 retained variants;
- one donor population, `N=500`, `Ne=10000`, `rho=0.02`, ancestry weight 1;
- mutation age `1e9`, seed `20260905`;
- 10,000 founders, batch size 4,096, phased HAP output;
- identical inputs and static packed-word partitioning for one and eight threads.

The two outputs were byte-identical: each was 120,797,168 bytes with MD5
`47c8f12ad6abbfc291ec7f86d49f5946`.

## Separation and exactness gate

Before enabling batches or worker threads, the former combined public workflow
was frozen on a two-chromosome founder/full-sib/half-sib/grandchild fixture.
Running the separated single-thread founder call followed by the separated
pedigree call reproduced founder H1/H2, pedigree H1/H2, segment records,
crossover records, sample order, and variant order exactly. There were zero
Mendelian inconsistencies. The corresponding files retained these entry hashes:

- phased pedigree HAP: `6b71725b956f2526e5efe39f430762be`;
- BIM: `e246025e0c287ab4fd4cd11ba45fd1ae`;
- FAM: `1071969cdaf231ace1cd1878bc77b196`;
- dosage BED: `21de9746f3b36c9e307b8e99e3a0e8b8`.

The bounded post-threading fixture then compared a single batch with unaligned,
64-aligned, and final-partial batch divisions at one, two, four, and eight
workers. HAP bytes, reloaded phases, segment audits, IDs, and metadata were
exactly equal in every case; changing the seed changed inheritance as expected.

## Timing

Wall times are seconds from one measured run per thread count. A 200 ms external
PowerShell sampler measured process working set.

| Stage | 1 thread | 8 threads |
|---|---:|---:|
| HAP open | 0.56 | 0.53 |
| HAP chromosome load | 2.71 | 2.57 |
| Founder event planning | 0.99 | 1.01 |
| Packed founder materialization | 37.56 | 11.76 |
| Positional HAP initialization | 0.28 | 0.30 |
| Positional batch writes | 4.07 | 4.25 |
| HAP metadata finalization/publication | 0.63 | 0.62 |
| Total | 47.05 | 21.31 |

Native materialization accelerated 3.19x (39.9% eight-thread efficiency). Total
time accelerated 2.21x (27.6% efficiency) because HAP loading, event planning,
and output remain sequential. Total throughputs were 10.22 and 22.57 million
individual-marker updates/second; materialization throughputs were 12.80 and
40.89 million updates/second.

## Memory

| Quantity | Bytes |
|---|---:|
| Packed reference H1+H2 | 6,155,264 |
| Largest packed founder batch H1+H2 | 49,242,112 |
| Reference plus batch biological payload | 55,397,376 |
| Largest current-batch event plan | 3,840,128 |
| Peak process working set, 1 thread | 164,601,856 |
| Peak process working set, 8 threads | 164,950,016 |

The output was completed in three sequential batches (4,096, 4,096, and 1,808
founders). Completed handles were closed before the next batch. Thus peak packed
working payload was bounded by the 4,096-founder batch even as the on-disk
10,000-founder payload reached 120,797,056 bytes plus the 112-byte HAP
header/table. No dense haplotype or genotype matrix was allocated.

The working-memory contract is
`O((R + B) M_c / 8 + S_B)`, while final HAP size remains
`2 * 8 * M * ceiling(N / 64)` bytes. Increasing the number of sequential batches
increases output and work, not the largest packed batch allocation.

## Extrapolation

The following linear projections use measured total throughput, not additional
runs. They exclude whole-genome scheduling effects, filesystem capacity/cache
changes, and any departure from this segment distribution, so they are planning
figures rather than performance claims.

| Founders x markers | 1-thread projection | 8-thread projection | Approximate HAP payload |
|---|---:|---:|---:|
| 100,000 x 20,000 | 3.26 min | 1.48 min | 500 MB |
| 100,000 x 100,000 | 16.31 min | 7.39 min | 2.50 GB |
| 1,000,000 x 6.8 million | 7.70 days | 3.49 days | 1.70 TB |

The last case necessarily needs explicit chromosome-level scheduling and very
large output infrastructure. Linear scaling from a single chromosome is not a
substitute for that qualification.

## Published HAPNEST context

Wharrie et al. report HAPNEST times for 100,000 samples of 15.0 minutes with one
thread and 6.3 minutes with eight threads at approximately 20,000 SNPs, and
41.6/11.7 minutes at approximately 100,000 SNPs. Their CentOS benchmark used
Intel Xeon E5-2680 v3 2.50 GHz processors, `Ne=500`, `rho=2.185`, different
reference/preparation and output paths, and averages of five trials. See the
[HAPNEST paper](https://academic.oup.com/bioinformatics/article/39/9/btad535/7255913).

Those published figures are not machine-, parameter-, or output-matched to this
gsim measurement. The gsim extrapolations therefore must not be presented as a
direct speed claim. HAPNEST's published eight-thread speedups are approximately
2.38x and 3.56x for its two variant scales; gsim observed 2.21x total and 3.19x
inside packed materialization on this distinct workload.

## Conclusion

Founder batching removes total-founder-count scaling from packed working memory,
and deterministic native word partitioning gives useful parallel acceleration.
The remaining measured founder bottleneck is packed materialization itself
(55.2% of eight-thread total), followed by sequential HAP I/O. The recommended
next milestone is focused C++ simplification around the now-stable separated
event-plan/materializer/sink architecture, with profiling of word-level copy
costs before considering SIMD or a different storage schedule.

Reproduce the measurement with `tools/benchmark/benchmark_parallel_founders.R`.
It accepts only an external prepared-reference prefix and external output
directory and never downloads data.
