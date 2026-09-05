# Experimental chromosome-wise packed simulation

## Boundary and ownership

The committed byte-matrix founder and pedigree implementations remain the
independent exact oracles.  `gsim` owns validation, RNG streams, HAPNEST and
meiosis event generation, pedigree traversal, chromosome orchestration, IDs,
maps, and audit records. The private `gsim::native` implementation provides
mutable one-bit storage and policy-free materializers. Packed and file-format
code originated from our `gbits` project at revision
`089bf1e69dea356248a62bb2d3bded4e84c64f7f` and retains its applicable MIT
notice. Metadata and VCF components originated in our own `gmat` project at
revision `33d6751abf00c41a15223459df7cae028d54b4b5` and are distributed as part
of gsim under gsim's GPL-3 license. No external library is used.

The interface is deliberately internal and chromosome-local:

```
.gsim_packed_backend()
.gsim_hapnest_founders_packed_chromosome(...)
.gsim_pedigree_genotypes_packed_chromosome(...)
.gsim_packed_unpack(haplotypes)       # bounded validation only
.gsim_packed_decode_genotypes(h1, h2) # optional validation output only
```

HAP-loaded packed reference handles now enter the same materializer directly;
the additive integration and explicit ownership contract are frozen in
`hap_reference_founders.md`.

R external pointers own private packed handles in the package DLL. Finalizers
close the matching object. There is no dynamic loader, external ABI, library
path, environment-variable lookup, or fallback implementation.

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

BGZF/BCF, missing or multiallelic alleles, compression, memory mapping,
threading, SIMD, and phenotype integration remain outside this path.
