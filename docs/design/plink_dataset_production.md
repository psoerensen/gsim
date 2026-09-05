# Experimental packed PLINK dataset production

This internal milestone composes three existing ownership layers. `gbits`
0.19 ABI 4 writes SNP-major BED bytes directly from chromosome-local packed
H1/H2 handles. `gmat` 0.2 experimental ABI 0 validates immutable ordered
variant/sample metadata and streams BIM/FAM text. `gsim` validates simulation
alignment, orchestrates chromosome appends, records provenance, and publishes
the coordinated triplet. There is no decoded genotype matrix and no duplicate
BED encoder.

## Identity and allele contract

The canonical sample order is supplied once when the dataset is created and
must exactly match both packed phase handles for every chromosome. The
canonical variant order is the literal sequence of chromosome appends and the
literal order within each appended chromosome. Variant IDs are globally
unique, and chromosome labels occupy one contiguous block.

Simulation bit `0` is REF; bit `1` is ALT; `H1 + H2` is ALT dosage. BIM A1 is
ALT and A2 is REF. The committed `gbits` BED decoder counts A1, so BED decoding
is exactly `H1 + H2`: dosage 0 is PLINK pair `11`, dosage 1 is `10`, and dosage
2 is `00`; `01` missing is never emitted. Heterozygous phase orientation has
no dosage effect. This interface never silently reverses metadata alleles.

BIM and FAM use the exact contracts in `gmat`'s
`docs/design/plink_metadata_writer.md`. In particular, BIM is six tab-separated
fields (chromosome, ID, cM, BP, ALT/A1, REF/A2), and FAM is family ID,
individual ID, paternal ID, maternal ID, sex, and missing phenotype `-9`.
Lines end with the single byte `\n` and there are no headers.

## Internal lifecycle

The unexported interface is conceptually:

```
.gsim_plink_dataset_create(gbits_backend, gmat_backend, prefix,
                           sample_metadata, overwrite,
                           buffer_variants, provenance)
.gsim_plink_dataset_append(dataset, chromosome, h1, h2,
                           variant_metadata)
.gsim_plink_dataset_finalize(dataset)
.gsim_plink_dataset_cancel(dataset)
```

Creation resolves and checks every required `gbits` and `gmat` symbol and ABI
before creating output. It validates and freezes FAM metadata, canonical paths,
and pre-existing-target policy. Each append validates its complete chromosome
metadata before synchronously appending packed BED records. Handles are not
retained, so callers can release a consumed chromosome immediately. Finalize
validates all counts, finishes BED, writes BIM/FAM, verifies record and byte
counts, and only then begins publication.

All preparation uses uniquely owned same-directory staging paths. With
overwrite disabled, any existing `.bed`, `.bim`, or `.fam` target is rejected
before staging. With overwrite enabled, existing targets remain untouched
until all three replacements are closed and verified. Publication first moves
existing targets to uniquely owned backups, then renames staged BED, BIM, and
FAM into place. On any partial failure, newly published files are removed and
backups are restored wherever the platform permits. Only transaction-owned
staging and backup paths are removed. A failed or cancelled object cannot
return a successful manifest.

Three independent renames cannot be filesystem-wide atomic. Readers may observe
a short missing or mixed triplet during publication. The contract is instead:
no final is touched before all replacements are ready; success exposes a fully
matched triplet; detected failure attempts complete rollback and reports any
rollback failure explicitly rather than silently claiming success.

The success manifest records canonical paths; sample and marker counts and
orders; chromosome ranges; observed/expected BED bytes; BIM/FAM record counts;
`bit 1 = ALT = BIM A1`; `gbits` and `gmat` versions/ABIs; provenance; buffer
bounds; and `publication_status = "published"`.

## Resource bound

BED conversion retains only `buffer_variants * ceil(N/4)` bytes beyond the two
packed chromosome handles. BIM and FAM writers retain one formatted record at
a time. Dataset state retains IDs and metadata but no biological allele matrix.
Peak biological payload is the current chromosome's two packed phases, which
may be released after append. Time is `O(NM + M + N)` and writer temporary
memory is `O(buffer_variants * ceil(N/4) + max_record_bytes)`, excluding the
validated identity metadata required for the manifest.
