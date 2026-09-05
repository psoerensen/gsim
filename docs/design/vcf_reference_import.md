# Phased VCF reference import and public simulation contract

## Ownership

`gmat` owns strict streaming parsing of ordinary text VCF and the resulting
sample/variant identity. `gbits` continues to own the marker-major one-bit H1
and H2 handles and HAP v1 storage. `gsim` aligns the caller's genetic map and
sample metadata, fills one chromosome's packed handles from bounded VCF record
buffers, coordinates HAP/BIM/FAM publication, and exposes the supported R
workflow. No component introduces another packed representation.

## Supported VCF subset

The importer accepts an uncompressed VCF text file with a valid `#CHROM`
header, at least one uniquely named diploid sample, contiguous chromosome
blocks, positive nondecreasing POS values within each chromosome, and
biallelic SNP records whose REF and ALT are distinct uppercase A/C/G/T bases.
FORMAT must contain exactly one GT field; its position is discovered for every
record. Every sample field must have the declared FORMAT arity and GT must be
exactly `0|0`, `0|1`, `1|0`, or `1|1`. Missing, unphased, haploid, polyploid,
multiallelic, symbolic, breakend, indel, non-ACGT, malformed, or inconsistent
records are errors. `.vcf.gz` and BCF are unsupported.

VCF ID is retained when it is nonempty and not `.`. Otherwise the importer
uses the UTF-8 string `CHROM:POS:REF:ALT`. Final IDs must be globally unique.
The descriptor records which IDs were generated.

## Alleles, phase, samples, and map

The left and right alleles of GT become H1 and H2, respectively. VCF REF is
bit 0 and BIM A2; VCF ALT is bit 1 and BIM A1. Consequently downstream BED
dosage is exactly `H1 + H2`, the ALT/A1 dosage. Alleles are never flipped,
frequency-normalized, imputed, phased, split, or discarded.

VCF sample order is the HAP bit order and FAM row order. With no caller sample
metadata, FID is `reference`, IID is the VCF sample ID, both parent IDs are 0,
sex is 0, and phenotype is -9. Supplied metadata are matched by exact IID and
must contain each VCF sample once with no extras.

The required map is a data frame with `chromosome`, `genetic_position_cm`, and
exactly one usable key: `variant_id` or `base_pair_position`. Each imported
variant must match exactly once, with no unused rows. cM values must be finite,
nonnegative, and nondecreasing within chromosome; BP must already be
nondecreasing in the VCF. There is no interpolation, extrapolation, or implicit
unit conversion during import. Founder generation consumes cumulative cM as
before. Public pedigree simulation explicitly converts the imported cumulative
cM to Morgans by division by 100 for its biological-meiosis model and records
that conversion in provenance.

## Streaming and public objects

`gmat` validates the complete record stream while retaining only sample and
variant metadata plus record offsets. During import it emits one `2 * N` byte
allele record at a time. `gsim` writes those bits directly into two
chromosome-local `gbits` handles, appends the handles to HAP v1, and releases
them before the next chromosome. No allele or dosage matrix is created.

`gsim_import_vcf(vcf, map, output, sample_metadata = NULL,
overwrite = FALSE)` publishes HAP/BIM/FAM and returns a lightweight
`gsim_reference` descriptor. `gsim_reference(prefix)` validates an existing
triplet and reconstructs that descriptor. `gsim_simulate()` loads one reference
chromosome at a time, applies the unchanged founder event plan and packed
materializer, applies the unchanged packed pedigree meiosis, writes HAP or BED,
and returns a compact manifest. Packed handles are caller-owned by `gsim` for
the duration of one chromosome and are deterministically closed on success or
failure.
