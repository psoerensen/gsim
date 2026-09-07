The complete scripts are available below. Website builds display code only;
they never run these scripts or download biological inputs.

## Local end-to-end workflow

[Download end_to_end_phenotype.R](files/end_to_end_phenotype.R) or
[read the complete script on GitHub](https://github.com/psoerensen/gsim/blob/main/inst/examples/end_to_end_phenotype.R).

With gsim and qgg already installed, run from the repository root:

```r
library(gsim)
source("inst/examples/end_to_end_phenotype.R")
```

The script creates a small synthetic phased VCF fixture locally, generates a
reusable founder population, transmits chromosomes through a pedigree, writes
HAP and BED, prepares a qgg Glist, and simulates phenotypes. It compares default,
causal-probability, variance-weight, and combined configurations. Its synthetic
LD scores illustrate alignment and parameterization, not real reference LD.

qgg is optional for gsim generally but required for this Glist workflow. The
phenotype stage uses bounded BED decoding and native accumulation; HAP
phenotype input is not supported. The phenotype contract details memory bounds. Compatibility with gsuite-generated
Glist objects is not established.

## 1000 Genomes chromosome 22

[Download 1000G_chr22.R](files/1000G_chr22.R) or
[read the complete script on GitHub](https://github.com/psoerensen/gsim/blob/main/inst/examples/1000G_chr22.R).

This external-data workflow requires internet access when the official IGSR
GRCh37 VCF and HapMap map are absent. The VCF download is approximately 196 MB.
The script documents `GSIM_1000G_DIR` for local input reuse. It imports a selected
region and named samples, creates unrelated founders, and simulates pedigree
descendants. This is a real 1000G reference workflow, unlike the synthetic local
fixture above. No data are bundled with the website.

The following command is an external-input workflow, not a website build step:

```r
source("inst/examples/1000G_chr22.R")
```

See [Getting started](getting-started.qmd) for the stage-by-stage workflow
sketches and [Validation and performance](validation.qmd) for the retained
scientific evidence, benchmark environments, and uncertainty. External HAPNEST
comparison remains deferred; pedigree meiosis remains single-threaded.

## Controlled gene sets

Run `source("inst/examples/gene_sets.R")` from the repository root with gsim
installed. The example obtains causal SNPs from a small existing simulation,
uses an explicitly synthetic mapping, and prints exact set composition and
shared genes. Sets are deliberately conditional on causal truth.

[Download gene_sets.R](files/gene_sets.R). See the
[authoritative contract](contracts/gene_sets.qmd) for input, sampling and truth
semantics and instructions for substituting a real mapping.
