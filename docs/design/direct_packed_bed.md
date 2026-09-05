# Experimental direct packed BED orchestration

`gsim` owns final chromosome, individual, and variant ordering, plus provenance.
The binary encoder and transactional sink are `gbits` 0.19 C ABI 4 facilities.
`gmat` is deliberately absent: BIM/FAM, persistent marker/sample metadata, and
allele-orientation policy remain its next-milestone responsibility.

The unexported orchestration surface is:

```
.gsim_bed_sink_create(backend, path, sample_ids, overwrite,
                      buffer_variants, provenance)
.gsim_bed_sink_append(sink, chromosome, h1, h2, variant_ids)
.gsim_bed_sink_finalize(sink)
.gsim_bed_sink_cancel(sink)
```

Creation validates one canonical sample order and creates no final BED file.
Each append requires both packed phases to carry exactly that sample order and
the supplied unique variant order. Chromosomes may be appended once and in any
caller-declared final order. The sink consumes the handles synchronously and
does not retain them, so a simulation driver may release each chromosome before
generating the next.

Finalization returns a compact manifest with path, sample and variant order,
chromosome blocks, counts, observed and expected byte size, conversion-buffer
bytes, `gbits` version/ABI, and caller provenance. It explicitly records that
BIM/FAM and the mapping between simulated alternative alleles and BIM allele 1
are deferred. No decoded genotypes are returned.

The dynamically loaded backend resolves every phased, BED-reader, and BED-sink
symbol and requires ABI 4. Missing or pre-0.19 libraries fail before output is
created. Ordinary `gsim`, raw founder/meiosis, and packed simulation remain free
of a required `gbits` runtime dependency until this experimental path is called.

The production call allocates only the already-existing two packed phases and
the `gbits` conversion buffer. `.gsim_gbits_bed_read_all()` is a bounded test-
only decoder using the committed `gbits` reader; it is not called by the writer.
Founder and meiosis event streams are not involved in BED coding and are
unchanged.
