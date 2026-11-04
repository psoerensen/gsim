# Breed-origin simulation

`gsim_boa_genotypes()` transmits caller-supplied, purebred founder haplotypes
through an explicitly supplied pedigree using the existing packed Mendelian
meiosis kernel. It propagates one-byte breed labels on the identical phase and
crossover trajectory. Founder breed labels are zero-based in the returned raw
origin matrices, with an explicit ordered breed vector. Chromosomes are processed
separately. Founder inputs and origins are currently materialized matrices;
alleles remain packed during inheritance and genetic-value accumulation.
`maximum_ancestry_bytes` admits the two origin buffers before allocation.

The first public contract requires at least two breeds, purebred founders, two
known parents per descendant, unique marker IDs and contiguous chromosomes with
nondecreasing Morgan positions. It inherits the no-interference, sex-independent
Poisson crossover model of the existing meiosis kernel. It does not add a
breed-origin inference algorithm: these labels are simulation truth.

`gsim_boa_trait_state()` freezes physical per-marker effect covariances, effects,
marker partitions, purebred-founder frequencies and a reference-breed LBP
parameterization. Breed-major channels can represent traits or longitudinal basis
coefficients. Block diagonal SNP covariance corresponds to UNCOR effects; a full
SPD covariance permits correlated breed effects. Independent disjoint partitions
may have different covariances. Up to 64 breed/channel axes are supported by the
reused identity-keyed normal generator. Entirely zero covariance is supported;
arbitrary singular nonzero effect covariance is not silently repaired.

Each marker's breed-b SNP contribution uses
`[(a1=b)(h1-p_b) + (a2=b)(h2-p_b)] beta_b`.
Non-reference LBP contributions use
`[I(a1=b)+I(a2=b)-2 f_b] theta_b`, where f_b is the pedigree breed fraction.
The reference effect is zero. `gsim_boa_values()` accumulates these directly from
packed alleles and returns separate breed, total SNP, total LBP and total genetic
values. It adds no pedigree approximation and never rescales a later cohort to
force a target heritability or genetic correlation. Reuse the frozen state and
base frequencies when comparing generations. Nominal marker-effect covariance
and realized individual genetic covariance have different meanings.

Use `gsim_records()` to add correlated residual and permanent environmental
effects, fixed contributions, missing-trait patterns or explicit longitudinal
incidence. Sampling and phenotype withholding occur after all genetic truth is
generated; hidden genotypes and phenotypes must remain separate from fitting
inputs. Meiosis, marker effects and record noise use separate identity/seed
streams without consuming R's RNG. Founder-panel generation is caller-owned.

## Validation and limits

The isolated Windows installation passed two focused tests (55 assertions):
phase/origin replay against continuous crossover events, packed versus independent
dense accumulation, exact purebred/F1 LBP zeros, equivalent reference changes,
repeatability, R RNG preservation and component-label invariance with equal
covariance inputs. The final isolated installation also passed 72 existing packed-simulation and 60 breeding/record assertions. That initial qualification did not run full R CMD check; the later
[full-package follow-up](../qualification/genetic_extensions.md#full-source-package-check-follow-up)
passed all 1,643 assertions and PDF manual checks, with zero errors/warnings and
one existing license-file note.

The gsuite developer recipes `tools/examples/boa_simulation_pilot.R` and
`tools/examples/boa_simulation_pilot.py` connect simulation to standalone fitting.
Their 102-animal, 1,000-marker, three-breed pilot includes F1, F2, backcross and
three-way offspring, two traits, two repeated visits, missing traits, 62 typed
animals and eight phenotype-withheld candidates. Founder haplotypes have synthetic
LD blocks and breed differentiation; they are not calibrated to a real breed.
Observed exports and truth are separate. The bounded export recipe uses private
HAP sink helpers and is a developer example, not a promised public exporter API.

This is a structured-meiosis robustness simulation, **not** an exact Gaussian
draw from the fitted partial-pedigree hybrid H covariance. The fitter's fixed
pedigree K_b and C0 are declared working-model inputs. Do not interpret a single
pilot as REML recovery, PEV coverage, scientific calibration or population-scale
performance qualification. An independent model-matched Gaussian calibration
study and repeated empirical-panel studies remain distinct follow-up work.
R BOA fitting remains deferred; gsim owns simulation only.

Both fixed-covariance standalone pilot modes passed total/breed/LBP selected PEV against matching native libraries. These numerical development runs used the existing test-only entitlement fixture; production trust was unchanged.
