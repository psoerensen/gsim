# HAPNEST founder-haplotype model

## Status and provenance

This document freezes the founder-only model implemented by `gsim`.  It is a
clean native reimplementation of the stochastic model in HAPNEST, inspected at
commit [`ba52da1a63cf609306ea92540b3d130fa1efd213`](https://github.com/intervene-EU-H2020/synthetic_data/tree/ba52da1a63cf609306ea92540b3d130fa1efd213),
and of the model described in Wharrie et al. (2023), Bioinformatics 39:btad535,
[doi:10.1093/bioinformatics/btad535](https://doi.org/10.1093/bioinformatics/btad535).
HAPNEST and `gsim` are both GPL-3 licensed.  No HAPNEST source is copied and no
Julia, Python, HAPNEST, container, downloader, or external process is used at
runtime.

The HAPNEST source is the executable scientific oracle when the paper and code
differ.  In particular, the paper's Figure 1 uses a non-strict mutation-age
inequality, whereas the implementation uses strict `T < mutation_age`.

## Inputs, meanings, and units

- The reference is a pair of variant-aligned `H1` and `H2` matrices with donor
  individuals on rows and 0/1 alleles on columns.  HAPNEST preprocessing takes
  the first and second phased VCF alleles for every sample into separate
  `Int8` matrices with the same row order
  ([`preprocessing/utils.jl` lines 100-128](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/preprocessing/utils.jl#L100-L128)).
  `metadata.haplotypes[p]` contains population-file `SampleID` values and
  `metadata.index_map` maps each `SampleID` to this shared donor-individual row
  axis
  ([`parameter_parsing.jl` lines 156-181](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/utils/parameter_parsing.jl#L156-L181)).
  The native interface requires the two matrices separately and one population
  label per donor-individual row.
- `N[p]` is the HAPNEST reference-panel population size used in the coalescent
  scale.  HAPNEST derives it as the number of donor sample IDs in population
  `p`, not the number of synthetic samples
  ([`parameter_parsing.jl` lines 156-181 and 253-275](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/utils/parameter_parsing.jl#L156-L181)).
  The native compatibility interface verifies that caller-supplied `N[p]`
  equals the number of donor-individual rows assigned to active population
  `p`.
- `Ne[p]` is the population-specific effective population size.  `rho[p]` is
  HAPNEST's population-specific recombination-rate calibration parameter
  ([default values](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/config.yaml#L74-L89)).
- Genetic positions are cumulative chromosome positions in centimorgans
  ([`GenomicMetadata` declaration](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/utils/parameter_parsing.jl#L41-L62)).
  `rho` must consequently be calibrated for this cM scale; the implementation
  does not silently convert Morgans, base pairs, or rates.
- Mutation ages and coalescent age `T` are in generations.  HAPNEST
  preprocessing converts source mutation ages from years using 25 years per
  generation and rounds to three decimals
  ([`preprocessing/utils.jl` lines 40-69](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/preprocessing/utils.jl#L40-L69)).

## Stochastic algorithm

For each requested individual, generate two haplotypes independently.  For
each chromosome block, start at its first variant and repeat until its last
variant has been filled:

1. Select donor population `p` from the configured ancestry weights.  HAPNEST
   uses `StatsBase.sample(..., Weights(...))` independently for every segment
   ([`genotype_algorithm.jl` lines 55-72](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L55-L72)).
   The native API accepts finite nonnegative relative weights with a positive
   sum and normalizes them once; this is probabilistically equivalent to
   HAPNEST's configuration percentages, which must sum to 100
   ([`parameter_parsing.jl` lines 206-240](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/utils/parameter_parsing.jl#L206-L240)).
2. Draw coalescent age
   `T ~ Gamma(shape = 2, scale = Ne[p] / N[p])`.  Julia
   `Distributions.Gamma(alpha, theta)` uses shape and scale, and the HAPNEST
   call is explicit
   ([`genotype_algorithm.jl` lines 7-13](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L7-L13)).
   Because the shape is exactly two, the native implementation samples the
   mathematically identical sum of two independent exponentials with scale
   `Ne[p] / N[p]`.
3. Draw genetic segment length
   `L ~ Exponential(scale = 1 / (2 * T * rho[p]))`
   ([`genotype_algorithm.jl` lines 16-22](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L16-L22)).
   Equivalently, the exponential rate is `2 * T * rho[p]`.  C++ samples it as
   `-scale * log(U)` for `U` uniform on the open interval `(0, 1)`.
4. Select uniformly one donor `SampleID` assigned to `p`; this is a reference
   individual, not an independently addressable haplotype.  Resolve its shared
   row using `index_map`.  When constructing synthetic haplotype 1, copy only
   from that row of `H1`; when constructing synthetic haplotype 2, copy only
   from the same-indexed row of `H2`.  The reference individual is sampled
   independently for every segment and output haplotype
   ([`genotype_algorithm.jl` lines 66-76](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L66-L76),
   [`write_output.jl` lines 56-93](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/write_output.jl#L56-L93)).
5. Let `objective = genetic_position[start] + L`.  Starting at `start`, advance
   while the current position is less than or equal to `objective` and the
   chromosome has another variant.  The inclusive endpoint is therefore the
   first variant whose position exceeds `objective`, or the chromosome's final
   variant if none does.  This deliberate overshoot, including both endpoints,
   follows [`genotype_algorithm.jl` lines 25-35 and 47-54](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L25-L35).
6. For every position in the inclusive segment, retain a copied alternative
   allele only when `T < mutation_age[position]`; otherwise write zero
   ([`write_output.jl` lines 39-51](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/write_output.jl#L39-L51)).
7. Continue at `endpoint + 1`.  A segment and its donor can never cross a
   chromosome boundary.  HAPNEST itself invokes generation separately for
   each chromosome
   ([`genotype_algorithm.jl` lines 146-163](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L146-L163));
   the native call accepts contiguous chromosome blocks and resets at every
   boundary.

Pair the two generated haplotypes without reconstructing phase.  Optional
diploid counts are the exact elementwise sum `h1 + h2`, as in
[`write_output.jl` lines 97-99](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/write_output.jl#L97-L99)
and the paper's Figure 1c.

The native argument `donor_phase` defaults to `"hapnest"`, and that is the only
mode implemented in this milestone.  A pooled union of H1/H2 haplotypes would
be a distinct extension, not HAPNEST compatibility; `donor_phase = "pooled"`
is rejected rather than silently substituted.

## RNG and batch policy

The native core is single-threaded.  It uses specified SplitMix64 streams and
does not use R's global RNG or implementation-defined C++ distribution
objects.  Each `(seed, global haplotype index, chromosome-block index)` tuple
derives an independent deterministic stream.  Uniform doubles use the leading
53 random bits plus a half-unit offset, giving values strictly inside `(0,1)`;
bounded donor indices use rejection sampling rather than biased remainder
mapping.  The algorithm and constants are part of the reproducibility
contract.

`individual_offset` supplies global individual numbering for an internal
batch.  Generating a range in one call or in consecutive batches with the same
seed and offsets produces identical haplotypes and segment records.  HAPNEST
batches only to construct and write output
([`genotype_algorithm.jl` lines 88-133](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/algorithms/genotype/genotype_algorithm.jl#L88-L133));
the native batch boundary is likewise not a biological event.

## Validation and ownership

The native interface rejects empty or dimensionally misaligned H1/H2 matrices,
misaligned donor/variant names, missing or nonbinary reference alleles, missing
donor labels, `N` inconsistent with donor-individual population counts,
duplicated/unnamed population parameters, unsupported donor-phase modes,
negative or all-zero ancestry weights, positive-weight populations without a
donor, nonpositive or nonfinite `N`, `Ne`, or `rho`, nonfinite or negative
mutation ages, and nonfinite genetic positions.  Chromosome labels must form
contiguous blocks and positions must be nondecreasing within each block.
Inputs are rejected rather than reordered because reordering only the genetic
distance vector, as HAPNEST currently does
([`parameter_parsing.jl` lines 185-197](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/utils/parameter_parsing.jl#L185-L197)),
could break alignment with haplotypes and mutation ages.  Missing values are
never imputed in the production generator.

`gsim` owns the simulation-facing R validation and this statistical copying
algorithm.  The focused interface returns R `raw` matrices (one byte per
allele/count), not dense double matrices and not a new persistent storage
format.  Reusable bit-packed binary/ternary representation and PLINK coding
belong to `gbits`; panel identity, map alignment, transformations, and LD
summaries belong to `gmat`.  No production BED/BIM/FAM writer currently exists
in those repositories, so this milestone adds none.  Mutation-age filtering
belongs to the simulation model in `gsim`; mutation-map ingestion,
interpolation, and variant alignment remain preprocessing concerns outside
this milestone.

## Historical copying is not pedigree meiosis

The segment switches above approximate historical coalescent/recombination
events to synthesize unrelated founder chromosomes.  They are not biological
meioses between known parents.  A later pedigree milestone must take each
founder's two already-phased chromosomes and perform new parental meioses on
chromosome maps while preserving pedigree IDs and parentage.  It must not call
this founder-copying process to create descendants.

## Proposed pedigree/meiosis milestone

The next milestone should add a separate, explicitly named meiosis layer; it
must not extend or overload the HAPNEST founder primitive.

1. Read the existing `gsim_pedigree` object in `canonical_order`, without
   altering `pedigree`, `canonical_order`, `external_order`, `mapping`, or any
   `gsim_pedigree_records()` contract.  Define pedigree founders as rows with
   both `sire` and `dam` missing, including later new founders.
2. Require an explicit one-to-one `founder_map` with columns `animal` and
   `founder_individual`.  Validate completeness, uniqueness, and exact ID
   matching before simulation.  Record, for every chromosome, the two source
   founder haplotypes used by each mapped animal.
3. For a nonfounder in topological order, form the paternal chromosome by a
   new biological meiosis between the sire's two phased homologues and form the
   maternal chromosome independently from the dam's two homologues.  Crossover
   locations use the chromosome genetic map and reset at every boundary.  RNG
   streams must be keyed by seed, child ID, parental role, and chromosome so
   full-sib, paternal-half-sib, and maternal-half-sib families share the
   appropriate parent but receive independent auditable gametes.
4. The existing pedigree permits one missing parent.  Such a contribution must
   never be silently fabricated.  Either the first bounded API rejects partial
   parentage, or a later extension supplies an explicit `ghost_founder_map`
   whose generated chromosomes and stable synthetic IDs are recorded while
   the public pedigree's missing parent remains missing.  This choice must be
   frozen before implementation.
5. Emit audit tables for founder mapping and every meiosis: child, parental
   role, recorded parent ID, chromosome, starting homologue, crossover
   coordinates, and packed output row.  Family type is derived from the
   unchanged sire/dam columns, so full- and half-sib relationships remain
   inspectable.
6. Add or use a `gbits`-owned packed builder/constructor and hand transient
   generated bit planes directly to `gbits::PackedMatrix`; do not materialize
   an animals-by-whole-genome R double matrix.  The natural representation is
   two binary packed matrices with identical animal/variant axes.  Any diploid
   ternary view is derived from those phased planes and never used to recreate
   phase.
7. PLINK output remains downstream.  A future `gbits` writer should encode BED,
   while `gmat`-owned aligned `MarkerMetadata`/`SampleMetadata` supplies BIM/FAM
   identity.  FAM paternal ID, maternal ID, and sex must come directly from the
   unchanged pedigree table (`0` only for missing parents/unknown sex under the
   chosen PLINK convention), with an exported audit mapping between FAM order,
   canonical pedigree order, external order, and packed rows.

Before implementing step 3, separately freeze and validate the biological
meiosis model (crossover-count process, interference assumption, map-unit
conversion, chromatid starting choice, zero-length intervals, and sex-specific
maps).  None of those choices can be inferred from the historical HAPNEST
copying process.
