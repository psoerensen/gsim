# Marker-level pedigree meiosis oracle

## Scope and scientific separation

This document freezes the provisional, internal marker-level inheritance model
implemented by `gsim`.  HAPNEST-compatible founder generation and pedigree
meiosis are separate processes.  The former models historical copying from a
reference panel, including `N`, `Ne`, `rho`, coalescent ages, segment lengths,
and mutation-age filtering.  The latter transmits one recombinant copy of each
known parent's two chromosomes to a child.  None of the HAPNEST demographic or
mutation parameters participates in pedigree meiosis.

Offspring mutation, gene conversion, segregation distortion, crossover
interference, sex-specific maps, and pedigree inbreeding adjustments are not
part of this milestone.

## Pedigree and founder contract

- The input must be an existing `gsim_pedigree` object.  Its
  `canonical_order` is the authoritative parent-before-offspring traversal;
  the externally ordered `pedigree` table is matched to it by animal ID.
- An animal with both `sire` and `dam` missing is a founder, including a later
  founder.  An animal with exactly one recorded parent is rejected, naming the
  animal and missing parental side.
- Every founder must occur exactly once, by ID, in each supplied founder `H1`
  and `H2` raw 0/1 matrix.  The two matrices must have identical row and column
  identifiers and must contain exactly the pedigree founder set.  Row position
  alone is never used for matching.
- For a nonfounder, output `H1` is always the gamete transmitted by its recorded
  sire (paternal phase) and output `H2` is always the gamete transmitted by its
  recorded dam (maternal phase).  A founder's supplied H1/H2 labels are retained
  as given.  Diploid counts, when requested, are formed only as `H1 + H2`.

## Chromosomes and map

Variants occupy columns in caller order and require unique, nonempty IDs.
Chromosome labels must form unique contiguous blocks.  `genetic_position` is a
finite, nondecreasing cumulative position in **Morgans within each chromosome**;
positions reset only at chromosome boundaries.  For a chromosome whose first
and last marker positions are `a` and `b`, its simulated genetic length is
`L = b - a`.  A one-marker or zero-span chromosome has length zero.

No crossover can cross a chromosome boundary.  Each parental gamete starts a
fresh chromosome-specific process.

## No-interference meiosis

For every nonfounder, parental side, and chromosome independently:

1. Draw `C ~ Poisson(L)`.  The native deterministic sampler counts unit-rate
   exponential waiting times until their cumulative value exceeds `L`; this is
   an exact Poisson-process construction and uses open-interval uniforms.
2. Conditional on `C`, draw `C` locations independently and uniformly on
   `(a, b)`, then sort them.  This is distributionally the specified conditional
   uniform-location model.
3. Draw the starting parental homologue with probability one half after the
   count and location draws.  `U < 0.5` selects parental H1; otherwise H2.
4. Traverse markers from left to right and alternate parental homologues at
   every crossover.

The exact boundary convention is: a crossover at position `x` switches the
source **before** assigning a marker at `x`.  Equivalently, the allele at marker
position `p` uses the starting homologue toggled once for every crossover with
`x <= p`.  Fixed-event materialization accepts boundary-equal locations so this
rule can be tested directly; continuously sampled production locations are
strictly inside `(a, b)`.  Coincident fixed crossovers toggle once each.

The sire's gamete is generated independently and becomes child H1.  The dam's
gamete becomes child H2.  Descendant H1/H2 rows are immediately eligible as the
parental homologues of later descendants.

## Deterministic streams and traversal

The implementation is single-threaded and does not touch R's global RNG.  Each
gamete has a SplitMix64 stream derived from the numeric seed plus stable UTF-8
hashes of the child animal ID and chromosome label and a paternal/maternal side
tag.  As in the founder oracle, string identity is the exact, unnormalized,
case-sensitive UTF-8 byte sequence from `Rf_translateCharUTF8()`.  FNV-1a uses
offset `0xcbf29ce484222325` and prime `0x00000100000001b3`; all components are
combined with the meiosis-specific SplitMix64 constants in native code.
Consequently streams do
not depend on external row order, output options, chromosome collection order,
or operational batch boundaries.  Adding unrelated animals whose IDs and
parentage do not replace existing animals does not change existing streams.
As with any finite hash, a theoretical 64-bit collision is possible.

Canonical animals are processed in order.  Optional internal `batch_size`
changes only loop divisions, never stream keys or traversal order.  Simulating
one chromosome alone gives the same transmissions as selecting that chromosome
from a multi-chromosome call because the chromosome label, not block index,
keys its stream.

## Internal API and returned object

`.gsim_pedigree_genotypes(pedigree, founder_haplotypes, chromosome,
genetic_position, seed, return_haplotypes = TRUE, return_genotypes = TRUE,
return_crossovers = TRUE)` is deliberately unexported.  `founder_haplotypes`
is a named list containing `h1 = <raw matrix>` and `h2 = <raw matrix>` with
founder IDs as row names and variant IDs as column names.  Additional fields
from the qualified founder-core result are ignored, so that result can be
passed directly after its founder IDs have been explicitly aligned to the
pedigree.

The returned `gsim_pedigree_genotypes` list contains optional canonical-order
`h1`, `h2`, and `genotypes` raw matrices; optional `crossover_audit` tables for
every meiosis and crossover; `sample_ids`; `variant_ids`; an ordered variant
map; chromosome-block metadata; canonical/external pedigree alignment; and
model/RNG settings.  Disabled large outputs are `NULL` and do not change any
random event.  Detailed crossover rows are not allocated when their audit is
disabled.

`.gsim_meiosis_materialize()` is a separate fixed-event primitive.  It accepts
caller-supplied starting homologue and sorted crossover positions and performs
no random draws, allowing the boundary and alternating-copy rules to be tested
against hand-calculated gametes.

## Ownership and future production path

This deterministic simulation primitive and its internal R interface belong in
`gsim`.  Its raw byte matrices are a small-data oracle and are not bit-packed.
A production implementation should put persistent bit-packed haplotypes,
materialization, and BED coding in `gbits`, while `gmat` owns persistent aligned
variant, genetic-map, and panel metadata where appropriate.  No storage or file
writer is added here.
