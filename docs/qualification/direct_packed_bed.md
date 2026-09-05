# Direct packed BED qualification

The qualification uses exact byte equality and the committed `gbits` BED
reader. No BIM/FAM files or whole-genome intermediate are involved.

## Exact fixtures

Independent reference encoding confirmed the complete `0x6c 0x1b 0x01`
header and every record byte for sample counts 1, 2, 3, 4, 5, 7, 8, 63, 64,
and 65. It covered dosage 0/1/2, all-zero and all-one phases, both
heterozygous orientations, alternating patterns, one and multiple markers,
packed-word padding, and final BED-byte padding. The four-sample dosage vector
`0,1,2,0` produced the hand-calculated record byte `0xcb`.

The phase-asymmetric founder, general founder, and multigenerational pedigree
fixtures had zero mismatches against raw-oracle genotype counts, packed bounded
decoding, and unpacked H1 + H2. The committed BED reader reported zero missing
calls. Sample, variant, chromosome, and allele-count mismatch counts were zero.
The pedigree included children, a grandchild, full siblings, paternal half
siblings, maternal half siblings, and a later founder.

Three chromosomes (`chrZ`, `01`, `CHR1`) with 1, 5, and 3 markers were generated,
appended, and released sequentially. Capacities of one and four records produced
identical bytes. Forward output followed forward order exactly; reverse output
contained the corresponding reversed chromosome blocks without changing any
within-chromosome calls. For seven individuals and nine variants, expected and
observed size were both

`3 + 9 * ceil(7 / 4) = 21 bytes`.

## Filesystem lifecycle

Tests confirmed refusal of an existing destination by default, preservation of
that file, explicit overwrite with the old file retained until finalization,
append-after-finalize and double-finalize rejection, invalid path rejection,
incompatible append failure, cancellation/destruction without publication, and
removal of the exact owned temporary file. No partial destination was published
after failure.

## Memory and bounded sanity run

For 513 individuals and 2,048 markers, the two packed inputs used 294,912 bytes,
the BED file used 264,195 bytes, a 64-record conversion buffer used 8,256 bytes,
and the fixed native lifecycle object occupied 160 bytes excluding dynamic path
strings. A decoded genotype matrix would contain 1,050,624 elements and was not
allocated by the writer.

A non-benchmark sanity run wrote 4,202,496 genotypes (513 by 8,192) three times
in 0.010, 0.030, and 0.010 seconds on the qualification host. Its packed input,
BED, and conversion-buffer sizes were 1,179,648, 1,056,771, and 8,256 bytes.
The coarse timings are recorded only to rule out accidental quadratic or per-
genotype-allocation behavior, not as a production throughput claim.

The writer is `O(N M)` in time and `O(B ceil(N/4))` temporary bytes for fixed
record capacity `B`; temporary memory is independent of chromosome and genome
marker counts.

Final focused totals were 13/13 `gbits` CTest programs, 49/49 direct-BED R
expectations, 71/71 packed-simulation expectations, 140/140 raw-founder
expectations, 80/80 raw-meiosis expectations, 51/51 existing `gsim()`
expectations, and 102/102 existing pedigree/record expectations. Fresh shared-
library and isolated R-package installations succeeded.
