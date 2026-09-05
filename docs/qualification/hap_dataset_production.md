# HAP/BIM/FAM phased-dataset production

## Frozen ownership and representation

This internal experimental path composes existing components rather than
redefining them. `gbits` 0.20.0 (ABI/SOVERSION 4) owns HAP v1 packed storage,
validation, sinks, readers, and loaded phase handles. `gmat` 0.3.0
(experimental ABI/SOVERSION 0) owns validated variant/sample records and exact
BIM/FAM reading and writing. `gsim` 0.6.0 owns simulation-facing identity,
chromosome orchestration, alignment checks, provenance, and staged three-file
publication.

HAP bits have exactly the established simulation meaning: zero is REF and one
is ALT. H1 and H2 remain separate; genotype dosage is only H1 + H2. BIM A1 is
ALT and A2 is REF. FAM and BIM order are respectively the sample-bit order and
the global marker order in HAP. Nothing frequency-normalizes or swaps alleles.

The exact binary layout is frozen in the `gbits` design document
`docs/design/hap_v1.md`: a 64-byte little-endian header with `HAP\\x01`, packed
chromosome payloads starting at byte 64, and a terminal table of 48-byte
entries. Each payload stores marker-major H1 followed by H2 using
`ceil(N / 64)` words per marker, sample `i` at LSB-first bit `i % 64`, with
zero padding.

## Internal lifecycle

`.gsim_hap_dataset_create()` validates and fixes final sample order through a
`gmat` metadata handle, creates same-directory staged paths, and opens the HAP
sink. `.gsim_hap_dataset_append()` requires tagged packed H1/H2 handles whose
sample and variant IDs exactly match the declared metadata, validates the BIM
block before writing, then appends packed words. Each chromosome may be
released after return. `.gsim_hap_dataset_finalize()` completes HAP, streams
BIM/FAM, verifies counts and exact HAP size, and calls the existing generalized
three-file backup-and-rollback publisher.

Existing destinations are refused by default. Explicit overwrite requires a
complete old triplet. Old files remain in place until all replacements are
ready. Publication renames three files sequentially, so no filesystem-wide
atomic observation is claimed; on failure, newly published files are removed
and the old triplet is restored wherever the platform permits. A rollback
failure is reported, never silently successful. Cancellation and prepublication
failure remove only owned staged/backup files.

`.gsim_hap_dataset_open()` asks `gmat` to read and revalidate BIM/FAM, asks
`gbits` to fully validate HAP v1, then checks HAP sample/marker counts and
zero-based chromosome ranges against contiguous BIM blocks.
`.gsim_hap_dataset_load_chromosome()` returns ordinary owning packed H1/H2
handles tagged with FAM/BIM IDs. They remain usable after reader closure and
can feed the direct BED sink without decoding.

HAP v1 deliberately contains no IDs or checksum. Coordinated creation proves
alignment by construction, while reopening detects shape and chromosome-block
mismatches. It cannot detect an independently substituted BIM/FAM with the
same counts and block sizes; the manifest states this limitation explicitly.

## Memory and complexity

For chromosome `c`, biological payload is
`2 * M_c * ceil(N / 64) * 8` bytes. Writing is `O(NM_c / 64)` packed-word I/O
and retains only the caller's current chromosome plus bounded stream state.
Reading is the same complexity and allocates exactly two chromosome-local
packed planes. Metadata is retained once as validated BIM/FAM records; text is
streamed record by record. No dense byte allele matrix, dosage matrix, second
packed chromosome, or whole-genome genotype payload is allocated.

The next bounded milestone is direct phased VCF/BCF import into chromosome-wise
packed HAP/BIM/FAM, with allele normalization and metadata provenance specified
before implementation.
