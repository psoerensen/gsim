# Packed simulation workflow and native architecture

This guide describes the supported public genotype workflow. `gsim` owns the
private packed storage and file I/O necessary for simulation; it is not a
general-purpose matrix library. The complete VCF-to-genotype workflow compiles
into the package shared library, with C++17, zlib, and native thread support.
No separately installed gbits or gmat library, dynamic backend loader, external
ABI, library-path setting, or fallback backend is used.

## Reference, base, and pedigree populations

The public packed workflow has three explicit biological roles:

1. `gsim_import_vcf()` or `gsim_reference()` describes a real phased reference
   panel stored as HAP/BIM/FAM.
2. `gsim_simulate_founders()` uses historical HAPNEST-compatible segment copying
   to create an unrelated phased synthetic base population. Its seed, `N`,
   `Ne`, `rho`, ancestry weights, and mutation ages belong only to this step.
3. `gsim_simulate_pedigree()` loads an existing phased base population and uses
   biological meiosis to transmit paternal H1 and maternal H2 to descendants.
   Its independent seed belongs only to meiosis.

The base population remains HAP because phase is required for recombination and
because the same founders may seed multiple pedigrees. BED is an unphased
ALT-dosage output and is available only after pedigree simulation.

Founder generation validates final IDs before writing. Supplying `n` creates
`syn1`, `syn2`, and so on; alternatively `founder_ids` preserves an explicit
unique order. Pedigree founder IDs must exactly equal the base FAM IDs. There is
no positional matching or implicit remapping.

The founder model is a clean native C++ implementation informed by
[HAPNEST](https://github.com/intervene-EU-H2020/synthetic_data) at pinned
revision
[`ba52da1a63cf609306ea92540b3d130fa1efd213`](https://github.com/intervene-EU-H2020/synthetic_data/tree/ba52da1a63cf609306ea92540b3d130fa1efd213)
and by Wharrie et al. (2023),
[Bioinformatics 39:btad535](https://doi.org/10.1093/bioinformatics/btad535).
No HAPNEST source is included. HAPNEST is
[GPL-3 licensed](https://github.com/intervene-EU-H2020/synthetic_data/blob/ba52da1a63cf609306ea92540b3d130fa1efd213/LICENSE),
as is gsim.

Internalized packed-haplotype, HAP v1, and BED components originated in our
`gbits` project at revision
`089bf1e69dea356248a62bb2d3bded4e84c64f7f` and retain the applicable MIT
notice installed in `COPYRIGHTS`. Internalized variant/sample metadata,
BIM/FAM, and strict phased-VCF components originated in our own `gmat` project
at revision `33d6751abf00c41a15223459df7cae028d54b4b5`; they are distributed as
part of gsim under gsim's GPL-3 license.

## Production call flow

The production path is intentionally direct:

```text
phased VCF records -> packed chromosome H1/H2 -> positional HAP sink
HAP chromosome -> founder event planner -> batched packed materializer -> HAP sink
founder HAP chromosome -> pedigree meiosis -> HAP or BED sink
```

`vcf_reader` owns bounded plain/gzip/BGZF line parsing. `metadata_storage`
validates BIM/FAM identities and serialization. `packed_chromosome` is the only
production allele representation and also owns interval, mutation-filter, and
gamete materialization. `hap_storage` and `bed_storage` own their respective
binary formats. `hapnest_founders.cpp` and `pedigree_meiosis.cpp` own stochastic
policy. `packed_r.cpp` and `metadata_r.cpp` are the two cohesive R/native
boundaries; registration is isolated in `init.c`.

There is one production founder algorithm and one production pedigree
algorithm. The raw byte founder and meiosis implementations remain independent,
bounded scientific oracles used by tests; they are not production routes.

## Ownership and lifetime

Native objects are concrete private `gsim::native` types owned directly by an R
external pointer. `native_r.h` centralizes pointer tags, validation, explicit
release, and finalization. HAP-loaded chromosome handles own their packed words
independently of the reader. Dataset sinks use RAII and same-directory staged
publication. R closes each reference and generated chromosome deterministically
after it is appended; finalizers are a safety net, not the normal lifecycle.

The canonical layout is marker-major, with `ceil(samples / 64)` little-endian
logical 64-bit words per marker and individual `i` at bit `i %% 64`. H1 and H2
are separate planes, bit zero is REF, bit one is ALT, and final-word padding is
always zero. Checked unsigned dimensions and offsets guard allocation and file
layout calculations.

## Founder batches and threads

The R layer constructs one compact columnar event plan per batch. Without
audits, it creates neither per-segment R objects nor an R data frame. One native
call validates and materializes that complete batch. Static workers receive
disjoint destination-word ranges, share read-only reference planes, never call
R, and report exceptions only after all workers join. RNG streams are keyed by
seed, global founder identity, phase, and exact chromosome label, so batch size,
thread count, and scheduling cannot affect H1/H2, audits, or HAP bytes.

The default requests 8,192 samples per batch. Explicit positive sizes are
rounded up to 64-sample words and capped at the population size; the final
batch may be partial. The positional HAP sink writes each completed batch to
its predetermined marker-major range without retaining earlier batches or creating temporary
per-batch datasets. Founder working memory is
`O((R + B) M_c / 8 + S_B)`, for reference samples `R`, effective batch samples
`B`, current-chromosome markers `M_c`, and current-batch segment records `S_B`.
No worker receives a private copy of the reference chromosome.

Pedigree meiosis intentionally remains single-threaded. Canonical
parent-before-offspring traversal is unchanged, child H1 is paternal, and child
H2 is maternal. Arbitrary pedigree sizes can place a parent and child in the
same mutable packed word, so safe generation-level parallel writes would need a
different storage design.

## Test-only native helpers

Small pack/unpack, word inspection, genotype decoding, fixed-event gamete and
segment primitives, raw founder generation, raw pedigree meiosis, and the
bounded BED reader remain registered for exact tests and diagnostics. They are
not invoked by the public production workflow and do not create an alternative
production representation.

## Shared deterministic event semantics

Founder packed generation requests the founder core in event-plan
mode: haplotype and genotype allocation are disabled while the canonical RNG
loop emits one record per copied segment.  Each record contains phase, donor
population and individual, inclusive start/end, coalescent age, and sampled
length.  The packed path selects H1 reference storage for output H1 and H2 for
output H2, then applies the strict `T < mutation_age` filter directly to packed
bits. Event plans and materialization share the same scientific contract.

Meiosis event sampling uses one native routine shared by the
raw gamete function and the packed event-plan function.  It returns the exact
starting homologue and sorted crossover positions.  It also resolves each
position with native `lower_bound` to the first marker at or to its right, so a
crossover at a marker switches before that marker.  The packed call receives
only these zero-based boundaries.  Exact audit equality between raw and packed paths is required in tests.

One packed pedigree call processes canonical animals parent-before-offspring
for one chromosome.  Founder IDs are matched explicitly, child H1 is the sire
gamete, child H2 is the dam gamete, and packed descendants immediately become
parents.  The one-known-parent rejection remains in force.  Operational batch
boundaries do not change event order or streams.

## Memory and complexity

For `I` individuals and `M_c` markers on the current chromosome, each phase
uses `8 * M_c * ceil(I / 64)` payload bytes.  Two byte matrices use
`2 * I * M_c` bytes.  With `I` divisible by 64 the biological allele payload
reduction is exactly eightfold. Handles and vectors add small fixed headers;
the final word is the only alignment padding.  Optional oracle/diagnostic decoded genotype
counts use `I * M_c` additional bytes and are not allocated by default.

The final on-disk HAP payload remains
`2 I M / 8` for total markers `M`, apart from word padding and format headers.
Founder materialization is `O(I M_c + S)` allele visits for `S` copied segments. Pedigree storage is `O(I M_c / 8)`, event storage is
`O(I + C)`, and materialization is `O(I M_c + C)`, where `C` is crossover
count.  Packed simulation does not unpack biological matrices or retain another
chromosome.  A caller can consume and release a chromosome's handles before
requesting the next chromosome.

## Dataset alignment and representations

HAP v1 preserves separate H1/H2 planes in marker-major packed storage;
BIM supplies ordered variant IDs, exact chromosome labels, cumulative cM,
physical positions, and alleles, while FAM supplies ordered sample identities
and pedigree metadata. `gsim_reference()` checks HAP counts and chromosome
ranges against BIM/FAM. Founder population labels are named by reference FAM
IID, mutation ages by reference BIM ID, and population parameters by population
name. The public interface aligns these identities before event generation.

VCF REF is bit 0/BIM A2 and ALT is bit 1/BIM A1. Left/right phased GT order
becomes H1/H2. SNP-major PLINK 1 BED stores unphased ALT dosage H1 + H2;
it cannot replace a phased base population for further meiosis. Founder output
is always HAP/BIM/FAM; pedigree output may be HAP/BIM/FAM or BED/BIM/FAM.
BIM cumulative cM is divided by 100 for pedigree meiosis in Morgans.
See the [VCF import contract](vcf-compressed-import.md) for selection, errors,
map interpolation, and parser memory bounds.

## Current limitations and phenotype handoff

The packed genotype stages do not allocate whole-genome dense matrices. The
subsequent BED-backed [phenotype workflow](marker_specific_phenotype.md#bounded-bed-accumulation)
uses bounded selected-column statistics and a native packed-record accumulation
for all traits. Optional summary statistics use bounded blocks too. qgg supplies
Glist construction, but is not required to accumulate from an existing supported
BED-backed Glist. HAP phenotype input and gsuite-generated Glist compatibility
are not established. Full-marker effects/truth and phenotype outputs still cost
memory; the phenotype contract specifies these bounds.

Plain, gzip, and BGZF VCF are supported through sequential scanning. BCF,
indexed regional access, missing or multiallelic retained alleles, HAP compression
or memory mapping, pedigree-parallel storage, and SIMD are unsupported here.
One full packed pedigree chromosome remains resident so descendants can become
parents. Detailed event audits add memory when requested.

The [founder scientific contract](hapnest_founder_model.md) specifies sampling,
mutation filtering, phase, and exact chromosome/seed identity. The
[pedigree scientific contract](pedigree_marker_meiosis.md) specifies inheritance
and crossover boundaries and distinguishes the byte oracle from public packed
simulation. Public [qualification evidence](../README.md#qualification-and-reproducibility)
retains measured environments and limitations; an external HAPNEST comparison
remains deferred.
