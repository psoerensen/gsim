# Bounded BED Glist phenotype accumulation

## Revision and scope

Measured locally on 2026-09-07: gsim 0.15.0 working-tree implementation on
`main`, based on `e2bac30ee0425f2eaeb5efa2f30d4e2b5f69743c` (0.14.0).
Entry worktree and index were clean. This is a focused integration qualification,
not a new model, biological-data benchmark, external HAPNEST comparison, or
qualification of HAP phenotype input or gsuite-generated Glist objects.
The [phenotype contract](../design/marker_specific_phenotype.md#bounded-bed-accumulation)
is authoritative for current semantics, supported metadata and memory bounds.

The old R path retained chromosome chunks, concatenated all causal columns,
imputed/standardized that matrix, and multiplied it twice when calibrating
effects. The replacement freezes statistics in capped blocks and uses one
native selected-marker scan for all traits. Full-marker effects and truth are
retained; `return_genotypes = TRUE` now errors for Glist input.

Read-only implementation references were gbits
`02298cbb30e5d88acc0759d608631a0559675840`, gmat
`e8ac1b24dc66accef83b8b6a7328edf5493769cb`, and qgg
`bfac8b2c388afb7ae1c88019bcfef8588f81aedb`. The reused implementation is gbits
BedReader checked record I/O and `xmat_scalar`, adapted privately without its
external adapters, transform/selection framework or scheduling machinery.
The existing gsim packed BED diagnostic reader now shares that reader.
[Notices](../../inst/COPYRIGHTS) retain gbits MIT, gmat and HAPNEST provenance.

## Parity and failure checks

Six results and post-call RNG states were frozen from the entry R sources before
replacement, using installed qgg 1.1.6. The independent fixture encoder creates
two 32-sample, 12-marker BED files with allele-asymmetric 0/1/2 dosage, missing
values, invariant and all-missing columns. Sixteen selected samples are reordered;
eight eligible markers are noncontiguous. Cases cover BayesR one/two traits,
BayesC non-unit weights without calibration, MAF-dependent effects, derived
weights using Glist MAF, and fixed effects without standardization. Baseline
comparisons included optional summary statistics.

Causal IDs, components, probability surfaces, variance weights and post-call
R RNG states matched exactly. The largest absolute difference across compared
B/B_causal, G, E, Y, variance/correlation quantities and residual covariance was
1.77635683940025e-15. All comparisons passed 1e-12 tolerance, appropriate for
small scalar versus BLAS sums and subsequent calibration. This is not a promise
of bitwise products or an error bound for arbitrary large/cancelling sums.
Returned B is calibrated when requested; raw effect identity was separately
verified with calibration disabled, not inferred from rounded rescaled effects.

[Focused tests](../../tests/testthat/test-glist-streaming.R) independently
encode small BEDs and compare to the retained in-memory W path. They cover:

- One/two traits, covariance, fixed effects, explicit and derived variance
  weights, exact raw B, causal/component identity and RNG, block sizes 1/3/64.
- Selected-sample imputation/SD, monomorphic/all-missing behavior, sample order,
  noncontiguous markers, filtered BIM mappings and different per-file FAM order.
- Native registration/use without loading qgg; rejected unsupported storage,
  bad allele metadata, absent/duplicate IDs, invalid native indices, >64-column
  native requests, non-SNP-major mode and truncated BED.
- Deterministic bounded custom callbacks and bounded optional summary statistics.
  A 130-marker/two-trait fixture checks exact legacy summary row identity across
  the private decode cap, alongside numerical agreement.

The native tests pass without qgg in loaded namespaces. The existing phenotype
regression includes the real installed qgg PLINK fixture (489 samples, 2,000
markers), so support is not based solely on a hand-built descriptor. The focused
end-to-end phenotype regression and BED sink tests pass too. BED sink verification
was included because its private diagnostic reader now shares checked I/O.
No full test suite, R CMD check, large simulation or external data download ran.

## One bounded resource probe

Reproduction: install into an isolated library using [repository instructions](../../AGENTS.md),
explicitly load that installed gsim, then run from the repository root:

```r
source("tools/qualification/glist_streaming.R")
result <- glist_streaming()
```

The [script](../../tools/qualification/glist_streaming.R) independently writes
one synthetic BED record at a time. Fixed dimensions: 1,024 samples, 4,096
eligible markers, three traits, block size 32. Causal probabilities select
64/512/4,096 markers exactly. Generation seed 9107, simulation seed 552.
No summary statistics in this resource probe; their cap and numerical behavior
are tested separately. The synthetic dosages do not represent real 1000G LD.

Environment: Windows NT 10.0.26200.0, Intel Core i7-1365U, 12 logical processors,
R 4.4.1 x64, Rtools44 GCC/G++ 13.2.0, default local R BLAS. Scalar accumulation
has no new threads or scheduling. Timings are single observations with Rprofmem
active, no repeated trials, cache control or confidence intervals; they are
integration resource evidence, not throughput estimates.

| Causal markers | Elapsed seconds | Avoided full causal payload | Returned B_causal | Total returned object |
| ---: | ---: | ---: | ---: | ---: |
| 64 | 0.14 | 524,288 B | 6,320 B | 2,572,240 B |
| 512 | 0.28 | 4,194,304 B | 45,744 B | 2,722,768 B |
| 4,096 | 2.41 | 33,554,432 B | 361,136 B | 3,926,992 B |

For every run, decoded-block capacity was 262,144 bytes (256 KiB), native packed
record storage 256 bytes, and the largest Rprofmem allocation 262,192 bytes
(including its allocation overhead). No full causal genotype allocation appeared.
Native code inspection confirms that multiplication allocates only an n ? traits
R result and one packed record; the registered decode boundary rejects more than
64 columns. Rprofmem does not observe C++ heap allocations, hence both checks
matter. Several bounded decoded/copy matrices may coexist; 256 KiB is not a
claim about total genotype-stage peak memory.

Full B occupied 361,136 bytes in each run; G/E/Y together 272,400 bytes; input
Glist metadata 329,392 bytes. These object sizes include names and R overhead,
can double-count shared strings, and exclude transient R copies and the internal
BIM/FAM plan. The observed whole-process peak working set was 130,052,096 bytes
(124.027 MiB), sampled every 100 ms across the single process containing fixture
generation and all three calls. Sampling can miss peaks; process memory also
includes R, package loading, allocator retention and profiling. It is not a
per-kernel allocation measurement. Generated fixtures and profiling/sampling artifacts were temporary; reusable
test fixtures and the bounded probe are retained.

## Standalone build and integration

One clean native build used `R CMD INSTALL --preclean` into an empty temporary
library outside the repository. To correct legacy summary row order after the
initial checks, R sources/help were refreshed in the same library; native make
reported nothing to rebuild. A `--no-libs` refresh attempt could not stage the
DLL and restored the prior installation; the normal refresh succeeded.
Explicit isolated loading and registered native symbols were verified.
DLL imports contain only Windows runtime libraries and R.dll: no gbits/gmat/qgg
shared-library dependency, external ABI discovery or configured sibling paths.
DESCRIPTION keeps qgg in Suggests and adds no dependency. Startup locale warnings
and the installed testthat build-version notice were environmental; test cases
reported no warnings, failures or skips in the final focused runs.

Public Rd was regenerated through roxygen2 from R sources. Website source
preparation continues to consume the authoritative Markdown/Rd; no deployment
or site publication is part of this milestone. Historical performance reports,
scientific contracts, single-threaded pedigree meiosis and deferred external
HAPNEST comparison retain their separate scope.

Implementation SHA-256 (working-tree bytes):

| File | SHA-256 |
| --- | --- |
| `R/gsim.R` | `83c2206144353b70657191b4a3e00ef6567afb4b7d008551f2516e61fcc7245c` |
| `R/gsim_internal.R` | `6921879e8229e843c4f9ed5dfcaaf2842a701623a771d6f4c2f49ddfebfbf19d` |
| `R/glist_bed_internal.R` | `0af31ca8f14746fcf066d177a6dde20e6903badcdba67e5013b250c02c87ec43` |
| `src/bed_reader.h` | `246d31f1aefbca2d0516c50a2018ed2758f5bc951181726b6b256f22d89c0475` |
| `src/bed_reader.cpp` | `0b2574f5bd45c533b58ecee2ae3e8f966e20782a4d21ecd8ea0026299660f24d` |
| `src/bed_phenotype_r.cpp` | `aed994ab8cef572bfb76eda6adcc227f477f1498fdb469297d7fabc6dca73bb8` |
