#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "metadata_storage.h"
#include "vcf_reader.h"
#include "native_r.h"

namespace {
namespace md = gsim::native::metadata;

enum class MetadataKind { variant, sample };

struct Metadata {
    explicit Metadata(md::ValidatedVariantMetadata value)
        : kind(MetadataKind::variant), variants(
              std::make_unique<md::ValidatedVariantMetadata>(std::move(value))) {}
    explicit Metadata(md::ValidatedSampleMetadata value)
        : kind(MetadataKind::sample), samples(
              std::make_unique<md::ValidatedSampleMetadata>(std::move(value))) {}
    MetadataKind kind;
    std::unique_ptr<md::ValidatedVariantMetadata> variants;
    std::unique_ptr<md::ValidatedSampleMetadata> samples;
};

struct VcfReader {
    VcfReader(const std::string& path,
              const std::vector<std::string>& selected_samples,
              const std::string& selected_chromosome, bool has_region,
              std::uint64_t region_start, std::uint64_t region_end,
              bool skip_unsupported)
        : value(path, selected_samples, selected_chromosome, has_region,
                region_start, region_end, skip_unsupported) {}
    md::PhasedVcfReader value;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

Metadata* require_metadata(SEXP pointer, MetadataKind kind) {
    Metadata* value = gsim::native::r::require<Metadata>(
        pointer, "native metadata handle");
    if (value->kind != kind) fail("native metadata handle has the wrong kind");
    return value;
}

VcfReader* require_vcf(SEXP pointer) {
    return gsim::native::r::require<VcfReader>(pointer, "native VCF reader");
}

std::string scalar_utf8(SEXP value, const char* field) {
    if (TYPEOF(value) != STRSXP || XLENGTH(value) != 1 ||
        STRING_ELT(value, 0) == NA_STRING) {
        fail(std::string(field) + " must be one nonmissing string");
    }
    return Rf_translateCharUTF8(STRING_ELT(value, 0));
}

std::vector<std::string> strings(SEXP values, const char* field) {
    if (TYPEOF(values) != STRSXP) {
        fail(std::string(field) + " must be a character vector");
    }
    std::vector<std::string> out(static_cast<std::size_t>(XLENGTH(values)));
    for (R_xlen_t i = 0; i < XLENGTH(values); ++i) {
        if (STRING_ELT(values, i) == NA_STRING) {
            fail(std::string(field) + " must not contain missing values");
        }
        out[static_cast<std::size_t>(i)] =
            Rf_translateCharUTF8(STRING_ELT(values, i));
    }
    return out;
}

SEXP make_metadata(Metadata* value) {
    return gsim::native::r::make_owned(value);
}

SEXP make_info(const md::MetadataWriteResult& info) {
    SEXP out = PROTECT(Rf_allocVector(REALSXP, 3));
    REAL(out)[0] = static_cast<double>(info.record_count);
    REAL(out)[1] = static_cast<double>(info.bytes_written);
    REAL(out)[2] = static_cast<double>(info.maximum_record_bytes);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("record_count"));
    SET_STRING_ELT(names, 1, Rf_mkChar("bytes_written"));
    SET_STRING_ELT(names, 2, Rf_mkChar("maximum_record_bytes"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

void set_names(SEXP value, const std::vector<const char*>& labels) {
    SEXP names = PROTECT(Rf_allocVector(
        STRSXP, static_cast<R_xlen_t>(labels.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(labels.size()); ++i) {
        SET_STRING_ELT(names, i,
                       Rf_mkChar(labels[static_cast<std::size_t>(i)]));
    }
    Rf_setAttrib(value, R_NamesSymbol, names);
    UNPROTECT(1);
}

void require_r_length(std::uint64_t count, const char* field) {
    if (count > static_cast<std::uint64_t>(R_XLEN_T_MAX)) {
        fail(std::string(field) + " exceeds R vector limits");
    }
}

}  // namespace

extern "C" SEXP C_gsim_metadata_variant_create(
    SEXP chromosome, SEXP ids, SEXP cm, SEXP bp, SEXP alt, SEXP ref) {
    try {
        const auto chromosomes = strings(chromosome, "chromosome");
        const auto identifiers = strings(ids, "variant IDs");
        const auto alternates = strings(alt, "alternate alleles");
        const auto references = strings(ref, "reference alleles");
        const std::size_t count = chromosomes.size();
        if (identifiers.size() != count || alternates.size() != count ||
            references.size() != count || TYPEOF(cm) != REALSXP ||
            TYPEOF(bp) != REALSXP || XLENGTH(cm) != static_cast<R_xlen_t>(count) ||
            XLENGTH(bp) != static_cast<R_xlen_t>(count)) {
            fail("variant metadata columns must have identical lengths and numeric positions");
        }
        std::vector<md::MarkerMetadata> records;
        records.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            const double physical = REAL(bp)[static_cast<R_xlen_t>(i)];
            if (!R_FINITE(physical) || physical < 1.0 ||
                physical > 9007199254740991.0 || std::floor(physical) != physical) {
                fail("base-pair positions must be exact positive integers not exceeding 2^53-1");
            }
            records.push_back({chromosomes[i], identifiers[i],
                               REAL(cm)[static_cast<R_xlen_t>(i)],
                               static_cast<std::uint64_t>(physical),
                               alternates[i], references[i]});
        }
        return make_metadata(new Metadata(
            md::ValidatedVariantMetadata(std::move(records))));
    } catch (const std::exception& ex) {
        Rf_error("native variant metadata: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_sample_create(
    SEXP family, SEXP ids, SEXP paternal, SEXP maternal, SEXP sex) {
    try {
        const auto families = strings(family, "family IDs");
        const auto identifiers = strings(ids, "individual IDs");
        const auto fathers = strings(paternal, "paternal IDs");
        const auto mothers = strings(maternal, "maternal IDs");
        const std::size_t count = families.size();
        if (identifiers.size() != count || fathers.size() != count ||
            mothers.size() != count || TYPEOF(sex) != INTSXP ||
            XLENGTH(sex) != static_cast<R_xlen_t>(count)) {
            fail("sample metadata columns must have identical lengths and integer sex");
        }
        std::vector<md::SimulationSampleMetadata> records;
        records.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            const int value = INTEGER(sex)[static_cast<R_xlen_t>(i)];
            if (value == NA_INTEGER || value < 0 || value > 2) {
                fail("sex must be 0 unknown, 1 male, or 2 female");
            }
            records.push_back({families[i], identifiers[i], fathers[i], mothers[i],
                               static_cast<std::uint32_t>(value)});
        }
        return make_metadata(new Metadata(
            md::ValidatedSampleMetadata(std::move(records))));
    } catch (const std::exception& ex) {
        Rf_error("native sample metadata: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_write_bim(SEXP pointer, SEXP path) {
    try {
        const auto* value = require_metadata(pointer, MetadataKind::variant);
        return make_info(value->variants->write_bim(
            scalar_utf8(path, "BIM path"), false));
    } catch (const std::exception& ex) {
        Rf_error("native BIM write: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_write_fam(SEXP pointer, SEXP path) {
    try {
        const auto* value = require_metadata(pointer, MetadataKind::sample);
        return make_info(value->samples->write_fam(
            scalar_utf8(path, "FAM path"), false));
    } catch (const std::exception& ex) {
        Rf_error("native FAM write: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_read_bim(SEXP path) {
    try {
        Metadata* value = new Metadata(md::ValidatedVariantMetadata::read_bim(
            scalar_utf8(path, "BIM path")));
        SEXP pointer = PROTECT(make_metadata(value));
        const auto& records = value->variants->records();
        require_r_length(records.size(), "BIM record count");
        const R_xlen_t n = static_cast<R_xlen_t>(records.size());
        SEXP chromosome = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP ids = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP cm = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP bp = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP alt = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP ref = PROTECT(Rf_allocVector(STRSXP, n));
        for (R_xlen_t i = 0; i < n; ++i) {
            const auto& record = records[static_cast<std::size_t>(i)];
            if (record.base_pair_position > 9007199254740991ULL) {
                fail("BIM base-pair position exceeds exact R integer range");
            }
            SET_STRING_ELT(chromosome, i, Rf_mkCharCE(record.chromosome.c_str(), CE_UTF8));
            SET_STRING_ELT(ids, i, Rf_mkCharCE(record.marker_id.c_str(), CE_UTF8));
            REAL(cm)[i] = record.genetic_distance_cm;
            REAL(bp)[i] = static_cast<double>(record.base_pair_position);
            SET_STRING_ELT(alt, i, Rf_mkCharCE(record.allele1.c_str(), CE_UTF8));
            SET_STRING_ELT(ref, i, Rf_mkCharCE(record.allele2.c_str(), CE_UTF8));
        }
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 7));
        SET_VECTOR_ELT(result, 0, pointer); SET_VECTOR_ELT(result, 1, chromosome);
        SET_VECTOR_ELT(result, 2, ids); SET_VECTOR_ELT(result, 3, cm);
        SET_VECTOR_ELT(result, 4, bp); SET_VECTOR_ELT(result, 5, alt);
        SET_VECTOR_ELT(result, 6, ref);
        set_names(result, {"pointer", "chromosome", "variant_id",
                           "genetic_position_cm", "base_pair_position", "alt", "ref"});
        UNPROTECT(8);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("native BIM read: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_read_fam(SEXP path) {
    try {
        Metadata* value = new Metadata(md::ValidatedSampleMetadata::read_fam(
            scalar_utf8(path, "FAM path")));
        SEXP pointer = PROTECT(make_metadata(value));
        const auto& records = value->samples->records();
        require_r_length(records.size(), "FAM record count");
        const R_xlen_t n = static_cast<R_xlen_t>(records.size());
        SEXP family = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP ids = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP paternal = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP maternal = PROTECT(Rf_allocVector(STRSXP, n));
        SEXP sex = PROTECT(Rf_allocVector(INTSXP, n));
        for (R_xlen_t i = 0; i < n; ++i) {
            const auto& record = records[static_cast<std::size_t>(i)];
            SET_STRING_ELT(family, i, Rf_mkCharCE(record.family_id.c_str(), CE_UTF8));
            SET_STRING_ELT(ids, i, Rf_mkCharCE(record.individual_id.c_str(), CE_UTF8));
            SET_STRING_ELT(paternal, i, Rf_mkCharCE(record.paternal_id.c_str(), CE_UTF8));
            SET_STRING_ELT(maternal, i, Rf_mkCharCE(record.maternal_id.c_str(), CE_UTF8));
            INTEGER(sex)[i] = static_cast<int>(record.sex);
        }
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 6));
        SET_VECTOR_ELT(result, 0, pointer); SET_VECTOR_ELT(result, 1, family);
        SET_VECTOR_ELT(result, 2, ids); SET_VECTOR_ELT(result, 3, paternal);
        SET_VECTOR_ELT(result, 4, maternal); SET_VECTOR_ELT(result, 5, sex);
        set_names(result, {"pointer", "family_id", "individual_id",
                           "paternal_id", "maternal_id", "sex"});
        UNPROTECT(7);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("native FAM read: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_open(
    SEXP path, SEXP selected_samples, SEXP selected_chromosome, SEXP region,
    SEXP unsupported) {
    try {
        const std::string source = scalar_utf8(path, "VCF path");
        std::vector<std::string> requested_samples;
        if (selected_samples != R_NilValue) {
            requested_samples = strings(selected_samples, "selected samples");
            if (requested_samples.empty()) fail("selected samples must not be empty");
        }
        std::string chromosome;
        if (selected_chromosome != R_NilValue) {
            chromosome = scalar_utf8(selected_chromosome, "selected chromosome");
            if (chromosome.empty()) fail("selected chromosome must not be empty");
        }
        const bool has_region = region != R_NilValue;
        std::uint64_t bounds[2] = {0u, 0u};
        if (has_region) {
            if ((TYPEOF(region) != REALSXP && TYPEOF(region) != INTSXP) ||
                XLENGTH(region) != 2) {
                fail("region must be a numeric vector of length two");
            }
            for (R_xlen_t i = 0; i < 2; ++i) {
                const double value = TYPEOF(region) == REALSXP
                    ? REAL(region)[i]
                    : (INTEGER(region)[i] == NA_INTEGER
                           ? NA_REAL
                           : static_cast<double>(INTEGER(region)[i]));
                if (!std::isfinite(value) || value <= 0.0 ||
                    value != std::floor(value) || value > 9007199254740991.0) {
                    fail("region bounds must be exact positive integers");
                }
                bounds[i] = static_cast<std::uint64_t>(value);
            }
            if (chromosome.empty() || bounds[0] > bounds[1]) {
                fail("region requires a chromosome and start <= end");
            }
        }
        const std::string unsupported_value = scalar_utf8(unsupported, "unsupported");
        if (unsupported_value != "skip" && unsupported_value != "error") {
            fail("unsupported must be 'skip' or 'error'");
        }

        VcfReader* reader = new VcfReader(
            source, requested_samples, chromosome, has_region, bounds[0], bounds[1],
            unsupported_value == "skip");
        SEXP pointer = PROTECT(gsim::native::r::make_owned(reader));

        const auto& sample_values = reader->value.samples();
        const auto& variant_values = reader->value.variants();
        const auto& block_values = reader->value.chromosomes();
        const auto& import = reader->value.report();
        require_r_length(sample_values.size(), "VCF sample count");
        require_r_length(variant_values.size(), "VCF variant count");
        require_r_length(block_values.size(), "VCF chromosome count");

        SEXP samples = PROTECT(Rf_allocVector(
            STRSXP, static_cast<R_xlen_t>(sample_values.size())));
        for (R_xlen_t i = 0; i < XLENGTH(samples); ++i) {
            SET_STRING_ELT(samples, i, Rf_mkCharCE(
                sample_values[static_cast<std::size_t>(i)].c_str(), CE_UTF8));
        }

        const R_xlen_t variant_count = static_cast<R_xlen_t>(variant_values.size());
        SEXP variant_chr = PROTECT(Rf_allocVector(STRSXP, variant_count));
        SEXP variant_id = PROTECT(Rf_allocVector(STRSXP, variant_count));
        SEXP variant_bp = PROTECT(Rf_allocVector(REALSXP, variant_count));
        SEXP variant_ref = PROTECT(Rf_allocVector(STRSXP, variant_count));
        SEXP variant_alt = PROTECT(Rf_allocVector(STRSXP, variant_count));
        SEXP generated = PROTECT(Rf_allocVector(LGLSXP, variant_count));
        for (R_xlen_t i = 0; i < variant_count; ++i) {
            const auto& value = variant_values[static_cast<std::size_t>(i)];
            if (value.base_pair_position > 9007199254740991ULL) {
                fail("VCF POS exceeds exact R numeric range");
            }
            SET_STRING_ELT(variant_chr, i, Rf_mkCharCE(value.chromosome.c_str(), CE_UTF8));
            SET_STRING_ELT(variant_id, i, Rf_mkCharCE(value.variant_id.c_str(), CE_UTF8));
            REAL(variant_bp)[i] = static_cast<double>(value.base_pair_position);
            SET_STRING_ELT(variant_ref, i,
                           Rf_mkCharCE(value.reference_allele.c_str(), CE_UTF8));
            SET_STRING_ELT(variant_alt, i,
                           Rf_mkCharCE(value.alternate_allele.c_str(), CE_UTF8));
            LOGICAL(generated)[i] = value.generated_id ? TRUE : FALSE;
        }

        const R_xlen_t block_count = static_cast<R_xlen_t>(block_values.size());
        SEXP block_chr = PROTECT(Rf_allocVector(STRSXP, block_count));
        SEXP block_first = PROTECT(Rf_allocVector(REALSXP, block_count));
        SEXP block_size = PROTECT(Rf_allocVector(REALSXP, block_count));
        for (R_xlen_t i = 0; i < block_count; ++i) {
            const auto& value = block_values[static_cast<std::size_t>(i)];
            SET_STRING_ELT(block_chr, i, Rf_mkCharCE(value.chromosome.c_str(), CE_UTF8));
            REAL(block_first)[i] = static_cast<double>(value.first_variant + 1u);
            REAL(block_size)[i] = static_cast<double>(value.variant_count);
        }

        SEXP variants = PROTECT(Rf_allocVector(VECSXP, 6));
        SET_VECTOR_ELT(variants, 0, variant_chr); SET_VECTOR_ELT(variants, 1, variant_id);
        SET_VECTOR_ELT(variants, 2, variant_bp); SET_VECTOR_ELT(variants, 3, variant_ref);
        SET_VECTOR_ELT(variants, 4, variant_alt); SET_VECTOR_ELT(variants, 5, generated);
        set_names(variants, {"chromosome", "variant_id", "base_pair_position",
                             "ref", "alt", "generated_id"});

        SEXP blocks = PROTECT(Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(blocks, 0, block_chr); SET_VECTOR_ELT(blocks, 1, block_first);
        SET_VECTOR_ELT(blocks, 2, block_size);
        set_names(blocks, {"chromosome", "first_variant", "variant_count"});

        SEXP report = PROTECT(Rf_allocVector(VECSXP, 16));
        SET_VECTOR_ELT(report, 0, Rf_mkString(import.input_type.c_str()));
        const std::uint64_t counts[] = {
            import.vcf_sample_count, import.selected_sample_count,
            import.total_records_scanned, import.retained_variants,
            import.outside_selected_chromosome, import.outside_selected_region,
            import.indels, import.multiallelic_records,
            import.symbolic_or_breakend_alleles, import.other_unsupported_alleles,
            import.missing_gt, import.unphased_gt, import.non_diploid_gt,
            import.duplicate_final_ids, import.maximum_parsing_buffer_bytes};
        for (R_xlen_t i = 0; i < 15; ++i) {
            SET_VECTOR_ELT(report, i + 1,
                           Rf_ScalarReal(static_cast<double>(counts[i])));
        }
        set_names(report, {"input_type", "vcf_sample_count", "selected_sample_count",
                           "total_records_scanned", "retained_variants",
                           "outside_selected_chromosome", "outside_selected_region",
                           "indels", "multiallelic_records",
                           "symbolic_or_breakend_alleles", "other_unsupported_alleles",
                           "missing_gt", "unphased_gt", "non_diploid_gt",
                           "duplicate_final_ids", "maximum_parsing_buffer_bytes"});

        SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
        SET_VECTOR_ELT(result, 0, pointer); SET_VECTOR_ELT(result, 1, samples);
        SET_VECTOR_ELT(result, 2, variants); SET_VECTOR_ELT(result, 3, blocks);
        SET_VECTOR_ELT(result, 4, report);
        set_names(result, {"pointer", "samples", "variants", "chromosomes", "report"});
        UNPROTECT(15);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("native VCF open: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_start(SEXP pointer, SEXP chromosome_index) {
    try {
        VcfReader* reader = require_vcf(pointer);
        if (TYPEOF(chromosome_index) != INTSXP || XLENGTH(chromosome_index) != 1 ||
            INTEGER(chromosome_index)[0] == NA_INTEGER ||
            INTEGER(chromosome_index)[0] < 1) {
            fail("chromosome index must be a positive integer");
        }
        reader->value.start_chromosome(static_cast<std::uint64_t>(
            INTEGER(chromosome_index)[0] - 1));
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native VCF start: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_next(SEXP pointer) {
    try {
        VcfReader* reader = require_vcf(pointer);
        const std::uint64_t sample_count = reader->value.samples().size();
        require_r_length(sample_count, "VCF sample count");
        std::vector<std::uint8_t> h1_values;
        std::vector<std::uint8_t> h2_values;
        std::uint64_t variant = 0u;
        if (!reader->value.next(variant, h1_values, h2_values)) return R_NilValue;
        if (variant >= 9007199254740991ULL) {
            fail("VCF variant index exceeds exact R numeric range");
        }
        const R_xlen_t count = static_cast<R_xlen_t>(sample_count);
        SEXP h1 = PROTECT(Rf_allocVector(RAWSXP, count));
        SEXP h2 = PROTECT(Rf_allocVector(RAWSXP, count));
        std::copy(h1_values.begin(), h1_values.end(), RAW(h1));
        std::copy(h2_values.begin(), h2_values.end(), RAW(h2));
        SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(result, 0, Rf_ScalarReal(static_cast<double>(variant + 1u)));
        SET_VECTOR_ELT(result, 1, h1); SET_VECTOR_ELT(result, 2, h2);
        set_names(result, {"variant_index", "h1", "h2"});
        UNPROTECT(3);
        return result;
    } catch (const std::exception& ex) {
        Rf_error("native VCF next: %s", ex.what());
    }
    return R_NilValue;
}

extern "C" SEXP C_gsim_metadata_vcf_close(SEXP pointer) {
    try {
        gsim::native::r::release<VcfReader>(pointer, "native VCF reader");
        return R_NilValue;
    } catch (const std::exception& ex) {
        Rf_error("native VCF close: %s", ex.what());
    }
    return R_NilValue;
}
