# Experimental chromosome-wise packed simulation

## Boundary and ownership

The committed byte-matrix founder and pedigree implementations remain the
independent exact oracles.  `gsim` owns validation, RNG streams, HAPNEST and
meiosis event generation, pedigree traversal, chromosome orchestration, IDs,
maps, and audit records.  `gbits` 0.18 provides the mutable one-bit storage and
policy-free materializers through stable C ABI 4.  `gmat` is not used here; it
remains the future owner of persistent panel, chromosome, variant, allele, and
map metadata.

The interface is deliberately internal and chromosome-local:

```
.gsim_gbits_backend(library)
.gsim_hapnest_founders_packed_chromosome(...)
.gsim_pedigree_genotypes_packed_chromosome(...)
.gsim_gbits_unpack(haplotypes)       # bounded validation only
.gsim_gbits_decode_genotypes(h1, h2) # optional validation output only
```

The library path is explicit or comes from `GSIM_GBITS_LIBRARY`.  R loads that
exact shared object, resolves every required phased-haplotype symbol, and the
native adapter requires ABI 4 before constructing a handle.  Missing files,
missing symbols, and ABI mismatches fail clearly.  The R external pointer owns
the `gbits` handle; its protected backend object keeps function addresses and
the loaded-library record alive until the handle finalizer calls the matching
close function.  The ordinary package and raw oracles remain usable without
`gbits`.

## Shared deterministic event semantics

Founder packed generation requests the committed founder core in event-plan
mode: haplotype and genotype allocation are disabled while the unchanged RNG
loop emits one record per copied segment.  Each record contains phase, donor
population and individual, inclusive start/end, coalescent age, and sampled
length.  The packed path selects H1 reference storage for output H1 and H2 for
output H2, then applies the strict `T < mutation_age` filter directly to packed
bits.  The raw oracle's default behavior and seeded outputs are unchanged.

Meiosis event sampling was factored into one native routine used by both the
raw gamete function and the packed event-plan function.  It returns the exact
starting homologue and sorted crossover positions.  It also resolves each
position with native `lower_bound` to the first marker at or to its right, so a
crossover at a marker switches before that marker.  The packed call receives
only these zero-based boundaries.  The draw order and SplitMix64 stream keys
are unchanged, and exact audit equality is required in tests.

One packed pedigree call processes canonical animals parent-before-offspring
for one chromosome.  Founder IDs are matched explicitly, child H1 is the sire
gamete, child H2 is the dam gamete, and packed descendants immediately become
parents.  The one-known-parent rejection remains in force.  Operational batch
boundaries do not change event order or streams.

## Memory and complexity

For `I` individuals and `M_c` markers on the current chromosome, each phase
uses `8 * M_c * ceil(I / 64)` payload bytes.  Two byte matrices use
`2 * I * M_c` bytes.  With `I` divisible by 64 the biological allele payload
reduction is exactly eightfold.  Handles and vectors add small fixed headers;
the final word is the only alignment padding.  Optional decoded genotype
counts use `I * M_c` additional bytes and are not allocated by default.

Founder reference storage is `O(R M_c / 8)`, generated storage is
`O(I M_c / 8)`, and event storage is `O(S)`, where `R` is reference-individual
count and `S` is copied-segment count.  Scalar materialization is `O(I M_c +
S)` allele visits.  Pedigree storage is `O(I M_c / 8)`, event storage is
`O(I + C)`, and materialization is `O(I M_c + C)`, where `C` is crossover
count.  No normal operation unpacks, allocates per marker, or retains another
chromosome.  A caller can consume and release a chromosome's handles before
requesting the next chromosome.

## Deferred work

VCF ingestion, persistent panel integration, multithreading,
SIMD, phenotype integration, and whole-genome retained objects remain outside
this milestone. `gbits` 0.19 now supplies the separately qualified direct
SNP-major BED sink described in `direct_packed_bed.md`; the subsequent
`plink_dataset_production.md` layer supplies BIM/FAM identity, ordering, and
allele-orientation metadata through `gmat`.
