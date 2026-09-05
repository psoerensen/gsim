# Standalone native architecture

`gsim` 0.12 compiles its complete supported VCF-to-simulation workflow into the
package shared library. It never searches for or loads `gbits` or `gmat`, and
the former `GSIM_GBITS_LIBRARY` and `GSIM_GMAT_LIBRARY` variables have no
meaning.

The minimal dependency inventory is marker-major one-bit H1/H2 storage and
interval operations; founder-segment and gamete materialization; HAP v1 and
SNP-major BED I/O; deterministic BIM/FAM validation and serialization; and the
strict streaming phased-VCF parser. Statistical event generation, chromosome
orchestration, pedigree policy, validation, and provenance remain gsim code.

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

One phase uses `8 * markers * ceiling(samples / 64)` bytes. Founder generation
is batch-local: only the current reference chromosome and one word-aligned
founder batch are resident, and batches are written positionally into the
marker-major HAP staging file. Work is otherwise chromosome-local. No
production path allocates a dense haplotype or genotype matrix.

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

The R layer constructs one compact columnar event plan per batch. One native
call validates and materializes that complete batch. Static workers receive
disjoint destination-word ranges, share read-only reference planes, never call
R, and report exceptions only after all workers join. RNG streams are keyed by
seed, global founder identity, phase, and exact chromosome label, so batch size,
thread count, and scheduling cannot affect H1/H2, audits, or HAP bytes.

Requested batches are aligned to 64-sample words except for the final partial
batch. The positional HAP sink writes each completed batch to its predetermined
marker-major range without retaining earlier batches or creating temporary
per-batch datasets. Founder working memory is
`O((R + B) M_c / 8 + S_B)`, for reference samples `R`, effective batch samples
`B`, current-chromosome markers `M_c`, and current-batch segment records `S_B`.
No worker receives a private copy of the reference chromosome.

Pedigree meiosis intentionally remains single-threaded. Canonical
parent-before-offspring traversal is unchanged, child H1 is paternal, and child
H2 is maternal. Arbitrary pedigree sizes can place a parent and child in the
same mutable packed word, so safe generation-level parallel writes would need a
larger storage design rather than a structural cleanup.

## Test-only native helpers

Small pack/unpack, word inspection, genotype decoding, fixed-event gamete and
segment primitives, raw founder generation, raw pedigree meiosis, and the
bounded BED reader remain registered for exact tests and diagnostics. They are
not invoked by the public production workflow and do not create an alternative
production representation.

## Consolidation in 0.12

The consolidation removed the internal packed and metadata C-ABI adapter pairs,
their function-pointer tables, status/version/ABI proxy calls, repeated opaque
handle layers, and R pseudo-backend contexts. It also removed two superseded
per-segment R/native entry points. Storage files were renamed by responsibility:
`packed_chromosome`, `hap_storage`, `bed_storage`, `metadata_storage`, and
`vcf_reader`.

Before consolidation the native tree contained 22 C/C++ headers or sources,
6,567 native lines, 51 registered routines, and two redundant adapter dispatch
layers. The final measurements are 18 files, 5,016 native lines, 47 registered
routines, and no adapter dispatch layer. The two scientific production paths
(founder generation and pedigree meiosis) are unchanged.
