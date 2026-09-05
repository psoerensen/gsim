# Direct compressed-VCF example using the official IGSR/1000 Genomes Phase 3
# GRCh37 chromosome 22 call set. The VCF download is approximately 196 MB and
# requires internet access. No bcftools, htslib, or other genomic software is
# used.

library(gsim)

vcf_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
  "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.",
  "20130502.genotypes.vcf.gz"
)
map_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/",
  "20130507_omni_recombination_rates/",
  "HapmapII_GRCh37_RecombinationHotspots.tar.gz"
)

data_dir <- Sys.getenv(
  "GSIM_1000G_DIR", file.path(tempdir(), "gsim-1000g-chr22")
)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
vcf_path <- file.path(data_dir, basename(vcf_url))
if (!file.exists(vcf_path)) {
  download.file(vcf_url, vcf_path, mode = "wb")
}

# This is the official GRCh37 map distributed by the 1000 Genomes project. Its
# source chromosome labels include "chr"; conversion to the VCF's exact "22"
# label is deliberate and visible here, never performed silently by gsim.
map_member <- "genetic_map_GRCh37_chr22.txt"
map_path <- file.path(data_dir, map_member)
if (!file.exists(map_path)) {
  map_archive <- file.path(data_dir, basename(map_url))
  if (!file.exists(map_archive)) {
    download.file(map_url, map_archive, mode = "wb")
  }
  utils::untar(map_archive, files = map_member, exdir = data_dir)
}
map_source <- utils::read.delim(map_path, check.names = FALSE)
stopifnot(identical(unique(map_source[["Chromosome"]]), "chr22"))
chr22_map <- data.frame(
  chromosome = "22",
  base_pair_position = map_source[["Position(bp)"]],
  genetic_position_cm = map_source[["Map(cM)"]]
)

selected_samples <- c(
  "HG00096", "HG00097", "HG00099", "HG00100",
  "HG00101", "HG00102", "HG00103", "HG00105"
)
output_dir <- file.path(tempdir(), "gsim-1000g-chr22-output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reference <- gsim_import_vcf(
  vcf = vcf_path,
  map = chr22_map,
  output = file.path(output_dir, "reference"),
  samples = selected_samples,
  chromosome = "22",
  region = c(20000000, 20500000),
  unsupported = "skip",
  overwrite = TRUE
)

# Setting every missing-parent probability to zero makes every nonfounder have
# both parents. Later founders, when drawn, have both parents missing.
pedigree <- gsim_pedigree(
  n_generations = 3,
  animals_per_generation = 6,
  founder_generations = 1,
  unknown_sire_probability = 0,
  unknown_dam_probability = 0,
  new_founder_probability = 0.1,
  permute_external_order = TRUE,
  seed = 22
)
populations <- stats::setNames(
  rep(c("P1", "P2"), each = 4), selected_samples
)
mutation_age <- stats::setNames(
  rep(1e9, length(reference$variant_ids)), reference$variant_ids
)
simulation_arguments <- list(
  reference = reference,
  pedigree = pedigree,
  populations = populations,
  ancestry_weights = c(P1 = 0.5, P2 = 0.5),
  mutation_age = mutation_age,
  N = c(P1 = 4, P2 = 4),
  Ne = c(P1 = 10000, P2 = 10000),
  rho = c(P1 = 0.02, P2 = 0.02),
  seed = 123
)
hap_result <- do.call(gsim_simulate, c(
  simulation_arguments,
  list(output = file.path(output_dir, "simulation-hap"),
       format = "hap", overwrite = TRUE)
))

# Set GSIM_1000G_WRITE_BED=true to produce BED/BIM/FAM from the same packed
# simulation. This is optional because it repeats the simulation deterministically.
if (identical(tolower(Sys.getenv("GSIM_1000G_WRITE_BED")), "true")) {
  bed_result <- do.call(gsim_simulate, c(
    simulation_arguments,
    list(output = file.path(output_dir, "simulation-bed"),
         format = "bed", overwrite = TRUE)
  ))
  print(list(paths = bed_result$paths,
             individual_count = bed_result$individual_count,
             marker_count = length(bed_result$variant_ids)))
}

print(reference$import)
print(list(paths = hap_result$paths,
           individual_count = hap_result$individual_count,
           marker_count = length(hap_result$variant_ids),
           simulation = hap_result$simulation))

# If an illustrative two-knot map is substituted for the official map above,
# the example remains a useful software smoke test, but its recombination
# distances are not scientifically calibrated and must not be used for
# scientific inference.
