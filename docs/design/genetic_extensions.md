# Dominance, environmental interactions and offspring mutation

## Explicit diploid genetic model

`gsim()` accepts `dominance` and `interaction` as marker-named effect vectors
(single trait) or marker-by-trait matrices. Matrices with trait names are aligned
to D1, D2, ...; omitted markers have zero effects. Environmental covariates are
finite scalars or vectors named by every selected sample. There is one supplied
covariate; multiple covariates are not inferred or fitted.

For raw diploid dosage x, the genetic value is

`G_it = sum_j X_ij B_jt + sum_j 1(x_ij=1) D_jt + e_i sum_j x_ij I_jt`.

X follows the existing additive standardization policy. Dominance uses uncentered
heterozygosity; interaction uses raw dosage, independently of additive
standardization. Neither term is automatically orthogonalized against the
additive component. Environmental zero and units are caller-owned; no main
environmental effect is silently added. Add one explicitly in record `fixed`
contributions if needed. These supplied effects consume no new RNG draws.

Nonzero extra-effect markers require hard calls 0, 1, 2 or NA. In `gsim()`, missing
calls use selected-sample observed mean dosage or observed mean heterozygosity,
respectively, and an entirely missing active extra-effect marker is rejected.
BED processing remains bounded to at most 64 decoded columns; native products
retain one physical BED record. Custom readers must preserve raw calls/missingness.

Returned `G_additive`, `G_dominance`, `G_slope`, and `G_interaction` separate truth;
G is additive + dominance + interaction, and Y is G + residual E. Unrequested
extra components are NULL, avoiding unused animal-by-trait buffers. Returned
`dominance`/`interaction` matrices contain the actual applied effects. When
`scale_effects=TRUE`, one multiplier per trait scales *all* effects and components
to the supplied total genetic variance vg. Consequently h2 is total genetic
heritability with extensions, not narrow-sense heritability. Covariances among
components are retained; their variances need not add to total variance.
`component` and B describe the additive architecture. The causal table covers
the union of additive, dominance and interaction loci, with separate effect
columns; component 1 can therefore occur at a dominance-only locus. Fixed models
may have zero additive effects. Random additive architectures retain their
existing no-additive-causal-draw rejection policy.

## Persistent cohorts and repeated environments

`gsim_trait_state(effects, center, dominance, interaction)` freezes all effects
and the additive dosage center. `gsim_genetic_values(..., environment=...)`
uses these without cohort-dependent recentering or rescaling. Missing interaction
calls use fixed center; missing dominance calls use HWE heterozygosity
`2 p (1-p)`, where p = center/2. This is an explicit imputation convention, not
an assertion that the selected cohort is in Hardy-Weinberg equilibrium.
Extended outputs carry a `genetic_components` attribute including slope values.

`gsim_records(..., interaction_values=slopes, environment=...)` applies a scalar
or record-named environment to supplied animal-by-trait/basis slopes through
the existing genetic design. It supports repeated visits and returns the separate
record interaction truth while preserving noise streams. Supply baseline genomic
values here, not values that already include the same environmental interaction.

## Offspring mutation

`gsim_simulate_pedigree()`, `gsim_simulate_cohort()` and `gsim_boa_genotypes()`
accept `mutation_rate`, default 0: the probability of independently flipping an
inherited 0/1 allele at each existing marker, in each transmitted copy. Mutation
occurs immediately after recombination and is inherited by later descendants.
Founders, crossover plans and inherited breed-origin labels remain unchanged.
Rate 1 flips every transmitted allele; rate 0 preserves the previous allele path.
This is recurrent binary mutation at supplied markers, not novel-site generation,
sequence mutation, a nucleotide-specific substitution process or founder mutation.

Mutation keys contain a separate domain, seed, child identity, chromosome,
parental side and within-chromosome marker index. Reordering unrelated animals,
cohort batching and changes in mutation rate do not alter crossover streams or
R's RNG. Marker catalog/order must remain fixed when comparing mutation streams.
The kernel modifies packed alleles without building a dense mutation matrix.
Mutation rate is retained in cohort and pedigree simulation metadata.

## Qualification scope

Focused checks cover independent genetic formulas, allele reversal of dominance,
common total-variance scaling, BED/dense parity with missing calls, persistent
fixed-base truth, record-specific environments with unchanged noise, zero/one
mutation endpoints, unchanged founders/crossovers, exact replay, cohort/full
pedigree agreement, and a binomial mutation-count check. Execution results are
recorded in [the focused qualification record](../qualification/genetic_extensions.md). These checks do not establish
population-specific biological calibration or fitted-model recovery. Polyploidy,
epistasis and general multi-environment covariance models are outside this change.
