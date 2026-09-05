# Reference, base, and pedigree populations

The public packed workflow has three explicit biological roles:

1. `gsim_import_vcf()` or `gsim_reference()` describes a real phased reference
   panel stored as HAP/BIM/FAM.
2. `gsim_simulate_founders()` uses historical HAPNEST-compatible segment copying
   to create an unrelated phased synthetic base population. Its seed, `N`,
   `Ne`, `rho`, ancestry weights, and mutation ages belong only to this step.
3. `gsim_simulate_pedigree()` loads an existing phased base population and uses
   biological meiosis to transmit paternal H1 and maternal H2 to descendants.
   Its independent seed belongs only to meiosis.

The base population remains HAP because phase is required for recombination and
because the same founders may seed multiple pedigrees. BED is an unphased
ALT-dosage output and is available only after pedigree simulation.

Founder generation validates final IDs before writing. Supplying `n` creates
`syn1`, `syn2`, and so on; alternatively `founder_ids` preserves an explicit
unique order. Pedigree founder IDs must exactly equal the base FAM IDs. There is
no positional matching, implicit remapping, or compatibility alias for the
removed combined `gsim_simulate()` operation.

Founder batches default to 8,192 requested samples. Explicit positive sizes are
rounded up to a 64-sample word boundary and capped at the population size; the
last batch may be partial. The current reference chromosome and current batch
are the only biological inputs resident during base generation. The positional
HAP sink preallocates the final chromosome plane and writes each batch into its
nonoverlapping marker-major word range, preserving HAP v1 bytes.

Each batch has one compact columnar event plan and one native materialization
call. The normal no-audit path does not promote those vectors to an R data
frame or construct per-segment R objects. Workers use static contiguous
packed-word partitions, write disjoint words, read the reference only, and never
call R. RNG streams remain functions of seed, global founder index, phase, and
exact chromosome label, so batching and thread scheduling cannot change
results.

Pedigree simulation remains single-threaded. Parent-before-offspring generations
are dependency-safe conceptually, but parents and children may share a mutable
64-bit marker word under arbitrary valid pedigree sizes and later-founder
placement. Parallel reads and writes of that word would be a C++ data race.
Changing storage or adding generation snapshots solely to parallelize this step
would be a larger, memory-increasing redesign; founder generation remains the
priority.

Founder working memory is `O((R + B) M_c / 8 + S_B)`, where `R` is reference
individuals, `B` is the aligned batch size, `M_c` is current-chromosome markers,
and `S_B` is current-batch events. Pedigree output still holds one full packed
chromosome because descendants must immediately become parents. No production
path creates a dense allele or dosage matrix.
