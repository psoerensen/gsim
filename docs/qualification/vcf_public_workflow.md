# Phased VCF import and public workflow qualification

## Scope and exact results

Qualification used `gbits` 0.20.0 ABI/SOVERSION 4, `gmat` 0.4.0
experimental ABI/SOVERSION 0, and `gsim` 0.8.0. The plain-text parser and C ABI
passed the complete nine-test gmat CTest suite; the unchanged gbits library
passed all 14 CTests. An isolated gsim installation passed all 11 testthat
files, including the strict import and public workflow fixtures.

The import fixtures cover all four GT classes, both heterozygous phase
orientations, GT in different FORMAT positions, multiple chromosomes,
generated IDs, exact ID- and BP-keyed maps, caller sample metadata, and sample
counts 1, 63, 64, and 65. Reloaded packed words and bounded unpacked H1/H2 were
exact. BIM/FAM bytes and HAP byte counts were exact. Repeated import bytes were
identical and the R global RNG state was unchanged.

Rejection tests cover unphased, missing, haploid, polyploid, multiallelic,
symbolic, breakend, indel, lowercase/non-ACGT, identical REF/ALT, absent or
duplicate GT, malformed sample arity/count, invalid POS, empty/duplicate final
IDs, duplicate samples, missing fileformat declaration, disjoint chromosome
blocks, `.vcf.gz`, unmatched maps, and negative genetic positions. Existing
transaction tests continue to cover staged HAP/BIM/FAM cancellation and
publication rollback; import map failures published zero files.

The end-to-end fixture contains two donor populations, five founders including
a later founder, full siblings, paternal and maternal half siblings, and a
grandchild on two chromosomes. Public HAP output matched the raw founder and
raw pedigree oracles exactly for H1 and H2. Independent BED decoding matched
raw `H1 + H2` at every genotype. Mendelian inconsistencies were 0. Fixed-seed
HAP bytes reproduced exactly, a changed seed diverged, HAP versus BED output
did not change simulated dosage, and the R global RNG was unchanged. Reversed
reference chromosome order and standalone-versus-surrounded simulation matched
exactly after chromosome-ID alignment. Existing founder tests also retain exact
batch invariance.

## Bounded memory and timing evidence

A Windows/R 4.4.1 sanity run imported 64 samples, 96 markers, and three
32-marker chromosomes from a 27,236-byte VCF, then simulated five founders and
ten pedigree animals. Import took 0.07 seconds and BED simulation took 0.09
seconds; these are sanity timings, not production benchmarks.

The HAP file was 1,744 bytes: 1,536 bytes packed biological payload and 208
bytes header/table. Byte-per-allele H1/H2 would require 12,288 bytes, exactly
eight times the packed payload. Peak loaded import-chromosome payload was 512
bytes and its combined native/R record allele buffers were 256 bytes (128 bytes
each). The avoided dense reference
haplotype and genotype payloads were 12,288 and 6,144 bytes. The reported peak
reference/founder/pedigree biological payload during simulation was 1,536
bytes. BED was 291 bytes with a 192-byte bounded conversion buffer.

Import time is `O(NM)` and retains `O(M)` identity/offset metadata plus
`O(N M_c / 8 + N)` allele storage/buffer for the current chromosome. Simulation
retains `O((R + F + P) M_c / 8 + S + C)` biological/event storage for reference
individuals, founders, pedigree animals, founder segments, and crossovers on
the current chromosome. Neither path allocates a whole-genome phased matrix or
a dosage matrix.

## Package check

Build-then-check completed with no errors or warnings when unavailable Suggests
were not forced and network repository checks were disabled. The remaining
NOTE is the repository's pre-existing explicit C++11 declaration, which is
essential to its native sources. All tests were also run explicitly with
testthat available.
