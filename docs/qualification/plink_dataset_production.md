# Packed PLINK dataset qualification

Qualification was run on Windows with R 4.4.1 and Rtools GCC 13.2.0 on
2026-09-05. Builds and installations used isolated temporary directories.
`gbits` was the unchanged committed 0.19.0 ABI/SOVERSION 4 library. `gmat`
0.2.0 retained experimental ABI/SOVERSION 0, and `gsim` was 0.5.0.

## Exact validation

The new `gmat.metadata_writer` test passed as part of all 8 `gmat` CTests.
It checks immutable metadata ownership, exact BIM/FAM bytes, ALT/A1 and
REF/A2 ordering, decimal formatting, UTF-8/space-containing paths, overwrite
refusal, C ABI lifecycle, and representative invalid metadata. The unchanged
complete `gbits` suite passed 13/13 tests.

The new `gsim` dataset file passed 76 exact expectations. Hand fixtures checked
the full byte sequences:

```
BED: 6c 1b 01 2b
BIM: chrZ\tv1\t0\t101\tG\tA\n
FAM: F\tS\t0\t0\t1\t-9\n ...
```

The BED byte represents phased 0|0, 0|1, 1|0, and 1|1 genotypes in sample
order and decodes to ALT dosages 0, 1, 1, and 2. Reversing the explicit
`bit1_allele`/ALT and `bit0_allele`/REF assertion is rejected.

Three chromosomes (`chrZ`, `01`, `CHR1`) with 2, 1, and 3 markers preserved
literal append and within-chromosome order. Expected and observed BED size was
15 bytes for 5 individuals and 6 markers. Founder and multigenerational
pedigree BED decoding had zero differences from both packed decoding and raw
oracle `H1 + H2`; missing-genotype count, sample-order mismatches,
variant-order mismatches, allele-count mismatches, and independently checked
Mendelian inconsistencies were all zero. The pedigree included full siblings,
paternal and maternal half siblings, a grandchild, and a later founder.

Invalid duplicate/missing identifiers, map values, BP positions, disjoint
chromosome blocks, unsupported or equal alleles, reversed orientation,
sample/variant order, self-parenting, one-known-parent records, absent parents,
parent-after-child order, contradictory sex, and nonmissing phenotype metadata
were rejected.

## Publication failures

Exact tests covered cancellation, failure before BED append, injected failure
after BED completion, injected failure after BIM completion, overwrite refusal,
successful overwrite, incomplete existing triplet rejection, and injected
failure after the first final rename. All staging files were removed after
cancellation or prepublication failure. The partial-publication test restored
all three pre-existing byte sequences exactly, left no owned temporary or
backup file, and returned no success manifest. Paths containing spaces and
non-ASCII characters succeeded on the qualification platform.

## Resource bound

For 257 individuals and chromosome marker counts 127, 65, and 1, with a
7-record BED conversion buffer:

- peak current-chromosome packed H1/H2 payload: 10,160 bytes;
- final BED: 12,548 bytes (exactly the expected size);
- BED conversion buffer: 455 bytes;
- BIM: 5,850 bytes, maximum formatted record 36 bytes;
- FAM: 4,775 bytes, maximum formatted record 19 bytes;
- bounded end-to-end elapsed time: 0.09 seconds.

The writer retained no decoded genotype matrix. A raw one-byte dosage matrix
for all 257 by 193 calls would require 49,601 biological-call bytes (and an R
integer matrix would require four times that), while a dense double matrix
would require 396,808 payload bytes. Each packed chromosome was released after
append; the 10,160-byte peak is the largest chromosome, not a whole-genome
haplotype total. Time is `O(NM + M + N)`. Conversion temporary memory is
`O(buffer_variants * ceil(N/4) + max_record_bytes)` plus validated identity
metadata; biological payload is `O(M_current * ceil(N/64))` per phase.

## Regression gates

A fresh isolated `gsim` install succeeded. All seven bounded test files passed:
existing `gsim`, pedigree/record, raw founder, raw pedigree meiosis, packed
simulation, BED sink, and new PLINK dataset tests. `git diff --check` passed in
each modified repository. No broad benchmark or large pedigree workload was
run.
