# Streaming phased VCF import contract

This document freezes the bounded import contract implemented in `gsim` 0.10.0.
The production reader is private native C++ compiled into `gsim`; it does not
invoke bcftools or link to a separately installed htslib, gbits, or gmat.

## Input and selection

The reader recognizes plain VCF and gzip by file content. A gzip stream whose
first member contains the BGZF `BC` extra subfield is reported as BGZF. zlib
validates each ordinary gzip or concatenated BGZF member while the reader
consumes bounded lines; no complete decompressed file or dense biological
matrix is created.

VCF sample IDs must be unique and nonempty. `samples = NULL` retains header
order. An explicit unique vector is matched exactly and establishes output
FAM/bit order, including noncontiguous and reordered selection. Only those
sample fields are parsed. `chromosome` uses exact label identity. `region` is
an inclusive pair of positive integer base-pair coordinates and requires one
selected chromosome. The initial reader scans sequentially and is not an
indexed-access implementation.

Records outside the selected chromosome or region are counted before genotype
parsing where possible. Records retained in VCF order must be biallelic SNPs
with distinct, uppercase A/C/G/T REF and ALT alleles, and every selected call
must have complete phased diploid GT equal to `0|0`, `0|1`, `1|0`, or `1|1`.
GT is found by name in FORMAT. With `unsupported = "skip"`, ordinary indels,
multiallelic, symbolic/breakend and other unsupported alleles, missing calls,
unphased calls, and non-diploid calls are categorized and skipped. With
`unsupported = "error"`, the first such record reports chromosome, position,
ID, and reason. Malformed headers, fields, POS, FORMAT/GT syntax, sample-field
counts, duplicate selected header IDs, corrupt/truncated compression, and
duplicate final IDs are always fatal.

Missing VCF IDs are deterministically replaced by `CHROM:POS:REF:ALT`.
Identity collisions after replacement are errors.

## Map and allele contract

A sparse map contains `chromosome`, positive `base_pair_position`, and finite
nonnegative cumulative `genetic_position_cm`. Positions are strictly increasing
and cM is nondecreasing within each contiguous exact-label chromosome block.
Every retained variant must lie from the first through final knot. Exact knots
copy the supplied cM value. Between knots, R double arithmetic evaluates, in
this order:

```text
cm_left + (cm_right - cm_left) *
  (bp - bp_left) / (bp_right - bp_left)
```

There is no extrapolation, clamping, physical-distance substitution, or silent
`chr` normalization. BIM uses the repository's deterministic decimal formatter
and records cumulative cM. Retained variants tied at the same physical position
receive the same interpolated value; map-knot positions themselves cannot tie.
The older exact-by-`variant_id` map representation
remains accepted because it shares the same unambiguous alignment path.

VCF left/right alleles become H1/H2. REF is packed bit 0 and BIM A2; ALT is
packed bit 1 and BIM A1. Downstream BED dosage is H1 + H2. Filtering never
changes phase or allele orientation.

## Memory and reporting

One chromosome's marker-major H1/H2 planes are retained at a time, using
`ceil(N / 64)` 64-bit words per marker and zero padding. The parser retains one
logical VCF line, a fixed 64 KiB decompression buffer, selected sample fields
for that record, and a two-byte-per-selected-sample allele buffer. Completed
chromosome handles are appended to HAP/BIM/FAM and released.

The import descriptor reports input type and size, original and selected sample
counts, deterministic record/category counts, output chromosome counts and
file sizes, peak packed bytes, maximum parsing-buffer bytes, generated-ID
count, map policy, and explicit declarations that no dense haplotype/genotype
matrix or uncompressed temporary VCF was created.
