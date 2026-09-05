# Direct packed HAP reference founders

## Ownership and integration boundary

This internal milestone composes the existing components without introducing a
second founder algorithm or packed representation. `gsim` 0.7.0 owns the
HAPNEST event model, validation, chromosome orchestration, IDs, audit records,
and deterministic cleanup. `gbits` 0.20.0 ABI/SOVERSION 4 already owns the
marker-major one-bit H1/H2 handles, HAP v1 reader, strict filtered segment copy,
and downstream HAP/BED sinks. `gmat` 0.3.0 remains the BIM/FAM metadata owner.
No `gbits` or `gmat` change is required.

The integration gap was in the former
`.gsim_hapnest_founders_packed_chromosome()` wrapper. It requested events by
passing complete byte reference matrices to the raw oracle and then packed
those matrices before copying. The native founder RNG loop is now factored into
one shared event-plan routine. The raw oracle supplies a byte-copy consumer;
the packed path supplies the existing gbits filtered-copy consumer. The event
draw order, SplitMix64 constants, phase-specific donor row, inclusive endpoint,
and strict `T < mutation_age` rule are unchanged.

The experimental chromosome-level interfaces are:

```
.gsim_hapnest_founders_packed_reference_chromosome(...)
.gsim_hapnest_founders_from_hap_chromosome(...)
```

The first accepts two caller-owned ordinary gbits phase handles. The second
accepts an open validated HAP/BIM/FAM reader, loads one chromosome, invokes the
first interface, and deterministically closes the two loaded reference handles
on success or failure. Generated handles are independent, caller-owned gbits
objects and can directly feed packed pedigree meiosis, HAP, or BED. The packed
pedigree oracle now accepts either its original byte founder matrices or a pair
of caller-owned packed founder handles; this does not alter its event model.

## Exact alignment contract

HAP sample count and chromosome ranges are validated against FAM/BIM when the
dataset is opened. Donor populations must be a named vector whose unique names
exactly equal FAM individual IDs; it is reordered only by explicit exact ID
matching. Genetic positions and mutation ages must likewise be uniquely named
by the selected BIM variant IDs. Missing, duplicate, or extra identifiers are
rejected.

The caller-provided genetic map must equal the BIM cM values under gmat's frozen
15-decimal fixed BIM serialization. This comparison is exact at the metadata
format boundary; it accommodates only information that BIM itself cannot retain
beyond that representation. Mutation ages are not stored in HAP/BIM/FAM and
remain explicit simulation metadata. `N` must equal the number of packed donor
rows assigned to each active population. Bit zero remains REF, bit one remains
ALT, and generated dosage downstream remains H1 + H2.

## Event and ownership contract

One native generator is the source of every founder event. Its event record is
linear in segment count and contains global synthetic individual/haplotype,
phase, stable chromosome label, inclusive marker interval, donor population and
individual, coalescent age, and sampled length. Packed copied/retained ALT
counts are obtained directly from the packed source during filtered copying;
no biological allele array is decoded.

HAP-loaded phase handles own their words independently of the reader. Closing a
reader therefore does not invalidate already loaded phases. Explicit packed
handle close clears the R external pointer; later use and double close fail
clearly. The HAP-reference helper closes its temporary reference handles, while
returned founders remain live. Sequential callers can consume and close both
reference and generated handles before loading the next chromosome.

## Memory and complexity

For `R` reference individuals, `I` generated founders, `M_c` current-chromosome
markers, and `S` segment records, the retained biological payload is

`O((R + I) * M_c / 8 + S)`.

Specifically the two reference planes use
`16 * M_c * ceil(R / 64)` bytes and the two generated planes use
`16 * M_c * ceil(I / 64)` bytes. No byte reference matrix or genotype matrix is
allocated. Scalar generation remains linear in event generation plus markers
visited by packed filtered copies. Audit-off mode uses the same event plan but
does not count ALT bits or retain the audit data frame in the result.

HAP v1 has no ID checksum. Dataset open validates counts and contiguous ranges;
identity alignment is established by coordinated BIM/FAM creation and exact IDs
held by the open dataset, not cryptographically proven against independently
substituted same-shape metadata.

## Deferred work

VCF/BCF import, compression, memory mapping, missing/multiallelic alleles,
threading, SIMD specialization, phenotype integration, and public APIs remain
outside this milestone. The next bounded milestone is direct phased biallelic
VCF import into chromosome-wise packed HAP/BIM/FAM.
