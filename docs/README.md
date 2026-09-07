# gsim documentation

Start with the [package README](../README.md) for installation and examples.
These documents describe the supported workflow; qualification reports retain
the versions, fixtures, and limits of their original measurements.

## Authoritative contracts

- [Packed simulation workflow and native architecture](design/packed_chromosome_simulation.md):
  reference, founder, and pedigree roles; HAP/BED/BIM/FAM; identity alignment;
  ownership, lifetime, deterministic batching/threading, and memory bounds.
- [HAPNEST-informed founder model](design/hapnest_founder_model.md): donor and
  phase sampling, mutation filtering, chromosome/seed identity, and attribution.
- [Pedigree inheritance and byte oracle](design/pedigree_marker_meiosis.md):
  parental transmission, crossover boundaries, maps, and public packed semantics.
- [VCF import](design/vcf-compressed-import.md): supported records, selection,
  skip/error policy, allele orientation, and map interpolation.
- [Marker-specific phenotype architecture](design/marker_specific_phenotype.md):
  causal probabilities, conditional mixture proportions, variance weights,
  metadata alignment, and returned truth.

The installed R help is authoritative for public arguments and return values:
[reference](../man/gsim_reference.Rd), [VCF import](../man/gsim_import_vcf.Rd),
[founders](../man/gsim_simulate_founders.Rd),
[pedigree simulation](../man/gsim_simulate_pedigree.Rd),
[phenotypes](../man/gsim.Rd), [printing](../man/print.gsim.Rd),
[pedigree domains](../man/gsim_pedigree.Rd), and
[record workloads](../man/gsim_pedigree_records.Rd).

Pedigree record generators supply deterministic solver workloads. Their latent
parent-average values are not exact draws from a pedigree relationship
covariance. Genotype-based phenotype simulation is a separate workflow.

## Examples and development

- [1000 Genomes chromosome 22](../inst/examples/1000G_chr22.R) is a complete
  GRCh37 example that downloads official reference/map inputs when absent.
- [End-to-end phenotypes](../inst/examples/end_to_end_phenotype.R) constructs a
  small synthetic local fixture. Run after `library(gsim)` with qgg installed;
  its illustrative LD scores are not estimates from real 1000G data.
- [Repository instructions](../AGENTS.md) describe installation, explicitly
  loading an isolated installed package for focused tests, and generated Rd.
- [Copyrights and provenance](../inst/COPYRIGHTS) retain HAPNEST attribution
  and applicable gbits/gmat notices.

qgg is optional generally and used to construct Glist objects. Supported BED
Glist accumulation is standalone and bounded; see the phenotype contract for
selection, storage, memory and custom-reader limits. HAP phenotype input and
gsuite-generated Glist compatibility are not established.

## Qualification and reproducibility

Reports retain the scope and revision of their own qualification. Historical recommendations are scoped to their measured revisions.

| Evidence | Scope and reproduction |
| --- | --- |
| [Bounded Glist accumulation](qualification/glist_bounded_accumulation.md) | 0.15.0 focused parity, malformed-input and bounded resource checks; [script](../tools/qualification/glist_streaming.R). |
| [Genotypes to phenotypes](qualification/genotypes_to_phenotypes.md) | Small local fixture, qgg dosage/ID parity and phenotype controls; [example](../inst/examples/end_to_end_phenotype.R) and [test](../tests/testthat/test-end-to-end-phenotype.R). |
| [Founder chromosome identity](qualification/hapnest_founder_chromosome_identity.md) | Exact byte/audit comparisons and intentional replacement of block-index streams; [tests](../tests/testthat/test-hapnest-founders.R). |
| [Pedigree meiosis](qualification/pedigree_marker_meiosis.md) | Byte oracle, Mendelian consistency, and uncertainty-based relationship experiment; [tests](../tests/testthat/test-pedigree-genotypes.R). |
| [Deferred HAPNEST comparison](qualification/hapnest_founder_comparison.md) | Pinned external oracle, prerequisites, proposed fixture and acceptance rules; comparison has not passed or run. |
| [Production chromosome 22 benchmark](qualification/production_benchmark_chr22.md) | Real 1000G inputs, checksums, measured revision/environment and timing uncertainty; [script](../tools/benchmark/benchmark_1000G_chr22.R). |
| [Batched/parallel founder benchmark](qualification/parallel_founder_benchmark_chr22.md) | Real 1000G prepared panel, plus a separately labelled synthetic 0.12 guard; [script](../tools/benchmark/benchmark_parallel_founders.R). |
| [Parallel scaling diagnosis](qualification/parallel_founder_scaling_diagnosis.md) | Synthetic 49,000-marker panel, not 1000G LD; counter limitations and removed external fixture are explicit; same [benchmark script](../tools/benchmark/benchmark_parallel_founders.R). |

The [50,000-animal record-workload script](../tools/qualification/pedigree_solver_workload.R)
reports construction sizes, times, and checksums, not covariance recovery or
solver performance. Current [tests](../tests/testthat.R), both examples, and all
benchmark/qualification scripts remain available. Benchmarks and broad
qualification require authorization; nothing in this index claims they were
rerun during documentation cleanup.
