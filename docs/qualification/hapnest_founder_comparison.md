# Deferred HAPNEST founder comparison

## Local status

The bounded distributional comparison was not run on 2026-09-04.  The audited
machine had no `julia` executable and the source-only HAPNEST checkout had no
reference-data directory.  No Julia installation, package environment,
container, Python `bed_reader` environment, or external biological dataset was
downloaded.  Exact native primitive tests remain the qualification evidence
available in this checkout.

## Pinned oracle and required fixture

Use HAPNEST commit `ba52da1a63cf609306ea92540b3d130fa1efd213` with Julia
1.6.4 (the version in HAPNEST's Dockerfile) and one Julia thread.  A tiny local
fixture must contain, in identical variant order:

- `H1` and `H2`: `Int8` donor-individual-by-variant matrices with only 0/1;
- donor `SampleID` and `Superpopulation` metadata;
- cumulative genetic positions in cM;
- mutation ages in generations;
- marker metadata needed only if the standard BED writer is exercised;
- ancestry weights plus per-population `N`, `Ne`, and `rho`;
- a fixed list of independent replicate seeds.

For the first comparison, use 128 variants at 0.01 cM spacing, two populations
with eight donor individuals each, deliberately asymmetric `H1` and `H2`
patterns, weights 0.35/0.65, and mutation ages spanning both sides of the
resulting coalescent-age distribution.  The asymmetry must be strong enough to
detect accidental pooled donor sampling: for example, include phase-diagnostic
variants fixed to 0 in every H1 row and 1 in every H2 row.  Use 20 replicates
of 500 synthetic individuals: bounded qualification, not a performance
benchmark.

The oracle driver must include the pinned, unmodified
`algorithms/genotype/genotype_algorithm.jl`, construct `GenomicMetadata`, call
`create_reference_table()` and `get_genostr()`, and export for each replicate:

- the genotype-count matrix;
- the reference table columns `H,I,S,E,P,Q,T`;
- copied genetic span computed from `S`, `E`, and the aligned map (the sampled
  exponential draw is not retained by the unmodified HAPNEST function);
- per-segment copied and retained alternative-allele counts.

The copied span and mutation counts are not retained by HAPNEST's standard
PLINK output; derive them from the saved reference table and the pinned
`add_mutations()` function without changing generation.
Record the HAPNEST commit, Julia version, `Distributions.jl` version, native
`gsim` commit, compiler, and all seeds beside the outputs.

## Commands once the prerequisites are already local

From a pinned HAPNEST checkout and an existing Julia 1.6.4 environment:

```text
git -C <hapnest-source> checkout --detach ba52da1a63cf609306ea92540b3d130fa1efd213
julia --project=<hapnest-source> -e "using Pkg; Pkg.instantiate(); Pkg.status()"
set JULIA_NUM_THREADS=1
julia --threads 1 --project=<hapnest-source> <oracle-driver.jl> <fixture-dir> <oracle-output-dir>
Rscript --vanilla <native-and-compare-driver.R> <fixture-dir> <oracle-output-dir> <comparison-output-dir>
```

On PowerShell, set the thread count with
`$env:JULIA_NUM_THREADS = '1'`.  `Pkg.instantiate()` and any Python package
setup require authorization if their caches are not already local.  The
production `gsim` path must never call these commands.

## Monte Carlo acceptance rules

Do not compare simulated haplotypes for identity.  For every summary, retain
the 20 replicate values from each engine and compare the difference of engine
means with its replicate-based Monte Carlo standard error.  Accept when the
absolute difference is at most four standard errors; if the estimated standard
error is zero, require exact equality.  Apply this rule to:

- per-marker allele-frequency error summaries (mean signed error and RMSE);
- genetic-distance-bin mean pairwise `r^2` (LD decay);
- ancestry segment proportions and donor-population usage;
- mean segment count and mean copied genetic span;
- mutation retention among copied alternative alleles;
- mean per-marker heterozygosity.

For ancestry and mutation proportions, also report the analytic pooled
binomial standard error.  For segment lengths, report medians and quartiles as
descriptive diagnostics and use a fixed-seed within-replicate bootstrap for
their differences; the 99% bootstrap interval should contain zero.  Inspect
plots and raw summaries, but do not replace these predeclared uncertainty-based
rules with arbitrary decimal tolerances after seeing results.
