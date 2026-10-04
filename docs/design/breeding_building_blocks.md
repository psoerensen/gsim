# Composable breeding simulation

These building blocks reuse the simulation setup from the related-population,
selective-genotyping and paired-pool experiments. They expose the simulation
mechanisms independently of any particular genomic predictor or experiment
driver. Previous experiment outputs remain historical evidence; replacing
their handwritten driver with these APIs does not silently reproduce or
overwrite their seed streams.

## Public workflow

| Function | Contract |
| --- | --- |
| `gsim_pedigree_from_table()` | Validate supplied identities, parents, sex and generation; reject absent parents, self-parentage and cycles; rebuild topological order, mapping and checksums. Preserve extra metadata. |
| `gsim_select()` | Choose a fixed quota randomly or using externally supplied named scores, with eligibility and strata. Return selection masks and conditional inclusion probabilities. Label deliberately oracle-based scores with `oracle=TRUE`. |
| `gsim_mate()` | Generate offspring parentage from disjoint sire/dam pools. Return actual parents used and a freshly validated combined pedigree. Preserve existing labels; new extra labels are missing. |
| `gsim_simulate_cohort()` | Reuse phased parents and write only new offspring, in separate HAP and BED triplets. Do not regenerate ancestral genotypes. |
| `gsim_trait_state()` | Freeze marker effects and dosage centers in a named, checksummed state. |
| `gsim_genetic_values()` | Accumulate `(dosage - fixed center) %*% effects` directly from packed BED, with fixed-center missing-call imputation. |
| `gsim_records()` | Generate observable scalar records from genomic trait values or basis coefficients; return truth separately. Support residual trait covariance within an observation unit and permanent-environment covariance within an animal. |
| `gsim_sample()` | Random recording/genotyping selection plus explicit retention of animals such as actual parents. Quotas are additional draws among non-retained animals. |
| `gsim_pool()` | Pool BED allele frequencies and observable phenotypes using separately declared DNA/phenotype weights, optional read-count error and streamed marker chunks. |
| `gsim_pool_covariance()` | Query a bounded block of pooled measurement-noise covariance, accounting for shared records, visits, traits and animals. |

All public help is generated from `R/breeding_*.R`. The package remains
standalone: it does not fit predictions, infer marker effects, estimate variance
components or depend on gsolve/greml. The caller can feed phenotype scores or
predictions into `gsim_select()`; the simulator cannot verify how those scores
were obtained. It never reads hidden genetic truth to choose parents implicitly.

## Genetic base and phenotype records

Apply one trait state to every generation. Effects are not redrawn or rescaled,
and dosage centers are not re-estimated in selected cohorts. This preserves
genetic gain on one reporting base. State creation accepts a named effect
matrix or the `B` matrix from a `gsim(standardize_W=FALSE)` simulation; creating the architecture
and applying it to later cohorts are separate operations.
Effects supplied as matrices must already be in raw-dosage units. A complete
standardized `gsim` result is rejected because the existing result does not
retain the dosage scaling required to convert its effects safely.

For multiple traits, one observation unit identifies simultaneous measurements
of one animal. Distinct visits use distinct units. Records from the same unit
share a vector of residual innovations; records from the same animal share
permanent-environment innovations across visits. Covariance matrices must be
positive definite or entirely zero, with at most 64 observed traits. Partial
trait patterns are allowed; the covariance dimensions match the set of traits
observed anywhere in the supplied record table.

For longitudinal genetic trajectories, supply an animal-by-basis-coefficient
genomic-value matrix and a record-by-coefficient `genetic_design`. For example,
an intercept/slope model uses rows `(1, time)`. Named design rows and columns
are aligned to record and coefficient identities. Permanent environment is a
trait-specific animal intercept, not a permanent random-regression trajectory.
The caller supplies fixed contributions such as herd/year/season effects.

Truth is separate from the observable `records$value`. Sampling and pooled
summaries use identities and observable values. There is no automatic selection
correction: retaining parents and using representative samples makes such
strategies testable, rather than certifying unbiased genomic estimates.

## Randomness and selection probabilities

New kernels use independent identity-keyed SplitMix64 streams. Mating keys use
child identity and parent role; record innovations use observation-unit or
animal identity, trait identity, and separate residual/permanent domains.
Reordering inputs and appending visits preserve existing draws when the trait
set, covariance and other model inputs remain fixed. Adding a covariance trait
or changing a trait's covariance factor is a model change, not an append-only
operation. No new kernel consumes R's global RNG.

Random selection ranks are keyed by identity and seed, and quotas apply within
each stratum. For a quota
`k` among `N` eligible animals in a stratum, the declared inclusion probability
is `k/N` under the pseudo-random sampling design. Score selection is
deterministic, with UTF-8 identity ties, and reports conditional 0/1 inclusion.
Those zeros do not support inverse-probability recovery of unobserved animals.
Retained parents have inclusion one; that alone does not make them
representative. Sampling quotas describe an additional sample; recording and
genotyping policies can be called separately with independent seeds.

Normal innovations use a specified Box-Muller transform. Binomial assay draws
use the compiled C++ standard-library distribution with keyed streams; their
reproducibility is scoped to the implementation, not guaranteed between
different toolchains. Meiosis reuses the existing chromosome-wise Poisson
no-interference model and identity streams. Neither model represents all
biological recombination or measurement mechanisms.

## Paired pools and measurement uncertainty

Each contribution identifies a pool, animal and observable record. DNA and
phenotype weights are normalized separately within pools. The DNA frequency
is `sum(normalized DNA weight * dosage / 2)` for the BIM A1/bit1 allele. The
phenotype is the separately weighted mean of observed values. Repeated animals
across visits may contribute more than once; their DNA weights add. Missing
contributing BED calls are rejected rather than silently imputed in pools.

The matching-weight diagnostic compares the total normalized weights on each
animal. Differing DNA/phenotype weights deliberately model, for example,
milk-volume-dependent DNA yield versus a differently weighted tank phenotype.
Even matching animal weights do not establish a common SNP regression row
when traits, visits, genetic designs or fixed effects differ. The caller must
preserve those model details when fitting paired summaries.

For read depth `d > 0`, the assay samples a binomial read count with probability
equal to the true pooled allele frequency. Conditional frequency variance is
`p * (1-p) / d`, returned in `assay_variance` (an attribute on streamed chunks).
This variance uses simulated true frequencies, rather than noisy estimates.
Errors are independent between pool/marker keys; extraction
bias, correlated reads and uncertain DNA weights are not automatically modeled.

Phenotype noise includes residual and permanent environment. Covariance
between pools is computed from shared unit/animal loadings on their covariance
factors. It excludes genetic covariance, fixed-effect uncertainty, and genotype
assay error. Genomic dependence belongs in the external prediction model.
Every pool reports its noise variance; covariance queries allow at most
256 pools per side. An overlapped individual/pool analysis must use the same
record/noise identities to construct its joint covariance externally; this
API currently queries pool-to-pool blocks only.

## Storage and execution bounds

The new parentage ordering, random/score selection, mating, stochastic records,
weight normalization, pooled values and covariance products run in C++.
Checksum reduction also runs natively using the existing serialization/checksum
convention, avoiding an R loop over large marker-effect or pedigree byte arrays.
Existing private native packed meiosis, HAP/BED I/O and BED effect accumulation
are reused. R validates contracts, aligns labels and orchestrates chromosomes
and bounded chunks. Small trait covariance factorization uses base R Cholesky;
no new general numerical library is introduced.

An incremental cohort holds one packed source-parent chromosome plus compact
used parents, a parent/offspring execution panel and the offspring output
panel. Memory depends on the current source cohort and offspring, not every
historical generation. Execution is single threaded and still has R calls per
animal/parental chromosome; this is a reuse-oriented simulation primitive,
not a newly qualified HPC population simulator. All used parents currently
must be available in one phased dataset. Cross-cohort mating requires caller
consolidation; there is no automatic multi-file parent archive.

The destination must be a new directory. Each triplet uses the existing
transactional publisher, but publishing both HAP and BED is not one atomic
transaction. On a failure between publications, a partial destination can
remain for inspection, and rerunning to that directory is refused. Existing
input datasets are never overwritten by this function.

BED truth accumulation holds one packed record and animal-by-trait output.
Pool frequencies hold one packed record and at most 64 pool-by-marker columns
per native call. Without a callback, the requested pool-by-marker result is
retained; use `allele_consumer` to stream large outputs. Covariance loadings
use storage proportional to contributions times trait count, without a
complete pool-by-pool matrix. The metadata and requested record/output tables
still occupy memory.

## Reproduction and qualification

The download-free [example](../../inst/examples/breeding_building_blocks.R)
creates a tiny synthetic phased reference, repeatedly selects parents,
generates only the next cohort, keeps a fixed marker architecture, produces
multitrait records, samples individuals with actual-parent retention, and
compares matched and mismatched pool weights. It illustrates composition,
not a statistical validation of a breeding strategy.

```r
library(gsim)
source(system.file("examples", "breeding_building_blocks.R", package="gsim"))
demo <- breeding_building_blocks(tempfile("breeding-demo-"))
demo$cohorts[[1]]$pool$summaries
```

Focused checks are in
[`test-breeding-blocks.R`](../../tests/testthat/test-breeding-blocks.R): parentage
failure cases, selection probabilities and order invariance, actual parents,
incremental/full-pedigree transmission parity across two cohorts, fixed-base
genomic truth, record append invariance, longitudinal design, native pool
products and dense measurement-noise covariance oracles. No new NAV-scale
benchmark or complete source-package qualification is claimed.

### Focused qualification, 2026-10-04

gsim 0.17.0 was installed with `R CMD INSTALL --preclean` into a fresh isolated
library outside the gsim repository, explicitly loaded from that library, and
tested on Windows with R 4.4.1 and Rtools44 GCC 13.2.0/C++17. The final native
build completed without compiler warnings. Generated help and namespace were
regenerated with roxygen2 7.3.3.

| Focused test file | Passing assertions |
| --- | ---: |
| `test-breeding-blocks.R` | 60 |
| `test-plink-dataset.R` | 77 |
| `test-hap-dataset.R` | 55 |
| `test-public-vcf-workflow.R` | 40 |
| `test-glist-streaming.R` | 394 |
| **Total** | **626** |

There were zero test failures, errors, test warnings or skips. R emitted its
existing startup locale warnings, and testthat reported that its installation
was built under R 4.4.3. These are distinct from the test results.

The installed download-free example also completed three offspring cohorts,
with 112 pedigree animals. Checks verified actual-parent retention, matching
and deliberately mismatching pool weights, and offspring dataset identities.
The combined focused checks and example took approximately 11 seconds; this
is execution evidence on tiny fixtures, not a population-scale benchmark.

The isolated installation and generated outputs stay in the existing gsuite
developer cache, outside source control. The previous gsim 0.16.0 study library
and archived selection/pool experiment outputs were preserved. Full
`R CMD check`, cross-toolchain assay parity, cross-cohort parent consolidation,
and large-population performance qualification were not run or established.

### Approximate development time

| Date | Work | Approximate time |
| --- | --- | ---: |
| 2026-10-04 | Native building blocks, packed cohort integration, public simulation contracts, generated help, example, targeted verification and documentation | 45 minutes |

This is an approximate elapsed development entry, including build/check time.
It is additional to the earlier prediction and selection experiments, whose
own development/evidence records remain separate.
