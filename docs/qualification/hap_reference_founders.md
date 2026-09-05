# Direct HAP-reference founder qualification

The focused qualification uses exact comparisons among the committed raw
reference oracle, the prior raw-to-packed wrapper, and HAP-loaded packed donor
handles. It covers phase-asymmetric donors, two populations, copied and removed
mutations, three chromosome labels whose lexical and supplied orders differ,
standalone/surrounded/reversed processing, explicit individual batches, output
options, HAP and BED round trips, reader and handle lifetimes, and direct packed
pedigree use across full siblings, paternal and maternal half siblings, and a
grandchild.

Acceptance is zero mismatches in H1, H2, H1+H2, complete segment records,
chromosome/variant/sample identity, HAP reloads, BED decoding, packed pedigree
output, and crossover audits, plus zero independently counted Mendelian
inconsistencies. Fixed seed and the R global RNG must remain invariant. Invalid
population, map, mutation, chromosome, dimension, and handle lifecycle inputs
must fail clearly.

## Results

The new focused suite passed 57/57 expectations. Complete segment records from
the raw oracle, the prior raw-to-packed path, and the direct HAP-loaded path had
zero differences, including donor population, donor individual, inclusive
intervals, coalescent ages, sampled lengths, and copied/retained ALT counts.
H1, H2, H1+H2, generated-founder HAP reloads, generated-founder BED decoding,
packed pedigree H1/H2/genotypes/crossover audits, pedigree HAP reloads, and
pedigree BED decoding each had zero mismatches. The phase-decisive fixture
produced only zeroes in generated H1 and only ones in generated H2. The
independent Mendelian inconsistency count was zero.

For three distinguishable chromosomes, forward versus reverse, standalone
versus surrounded, removal/addition of other chromosomes, one versus two
`individual_offset` batches, and output/audit toggles each had zero mismatches.
An informative different seed diverged. Closing the HAP reader did not affect
independently owned loaded handles; explicit handle close made use-after-close
and double-close fail as specified. The R global RNG state was unchanged.

The reproducible bounded accounting script reported:

| Chromosome | Markers | Reference packed H1/H2 | Generated packed H1/H2 | Event records | Peak biological payload | Raw reference avoided | Dense genotype avoided | Segments |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `Z-2` | 129 | 4,128 B | 4,128 B | 535,344 B | 8,256 B | 33,024 B | 16,512 B | 6,660 |
| `01` | 65 | 2,080 B | 2,080 B | 397,264 B | 4,160 B | 16,640 B | 8,320 B | 4,934 |
| `chrA` | 33 | 1,056 B | 1,056 B | 272,728 B | 2,112 B | 8,448 B | 4,224 B | 3,377 |

The 128-sample phase payloads are exactly one eighth of byte-per-allele H1/H2.
Only one chromosome's reference and generated payloads were live in each loop;
each generated pair was explicitly closed before loading the next chromosome.
The audit sizes include R data-frame overhead and scale with segments. No
decoded genotype matrix was requested.

The final dependency/package gates passed 14/14 gbits CTests, the committed
gmat metadata test, an isolated native `gsim` 0.7.0 installation, and every
`gsim` test file: BED sink, HAP dataset, PLINK dataset, packed simulation, raw
founder, raw pedigree meiosis, ordinary `gsim()`, pedigree, and the new direct
HAP-reference integration. No broad benchmark or external data was used.
