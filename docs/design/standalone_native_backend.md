# Standalone native backend

`gsim` 0.9 compiles its complete supported VCF-to-simulation workflow into the
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

One phase uses `8 * markers * ceiling(samples / 64)` bytes. Work is
chromosome-local: import, founder generation, meiosis, HAP writing, and BED
writing release each chromosome before the next. No production path allocates
a dense haplotype or genotype matrix.
