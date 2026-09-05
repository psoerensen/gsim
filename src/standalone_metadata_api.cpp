#include "standalone_metadata_api.h"

#include "standalone_metadata_error.h"
#include "standalone_metadata.h"
#include "standalone_vcf.h"


#include <algorithm>
#include <exception>
#include <limits>
#include <new>
#include <string>
#include <utility>
#include <vector>

struct native_metadata_variant_metadata_handle {
  explicit native_metadata_variant_metadata_handle(
      std::vector<gsim::native::metadata::MarkerMetadata> records)
      : value(std::move(records)) {}
  explicit native_metadata_variant_metadata_handle(gsim::native::metadata::ValidatedVariantMetadata records)
      : value(std::move(records)) {}
  gsim::native::metadata::ValidatedVariantMetadata value;
};

struct native_metadata_sample_metadata_handle {
  explicit native_metadata_sample_metadata_handle(
      std::vector<gsim::native::metadata::SimulationSampleMetadata> records)
      : value(std::move(records)) {}
  explicit native_metadata_sample_metadata_handle(gsim::native::metadata::ValidatedSampleMetadata records)
      : value(std::move(records)) {}
  gsim::native::metadata::ValidatedSampleMetadata value;
};

struct native_metadata_phased_vcf_reader_handle {
  explicit native_metadata_phased_vcf_reader_handle(const std::string& source)
      : value(source) {}
  gsim::native::metadata::PhasedVcfReader value;
};

namespace {

thread_local std::string last_error;

template <typename Function>
native_metadata_status protect(Function&& function) noexcept {
  try {
    function();
    last_error.clear();
    return NATIVE_METADATA_STATUS_SUCCESS;
  } catch (const gsim::native::metadata::Error& error) {
    last_error = error.what();
    return error.code() == gsim::native::metadata::StatusCode::internal_error
               ? NATIVE_METADATA_STATUS_IO_ERROR
               : NATIVE_METADATA_STATUS_INVALID_ARGUMENT;
  } catch (const std::bad_alloc&) {
    last_error = "metadata allocation failed";
    return NATIVE_METADATA_STATUS_INTERNAL_ERROR;
  } catch (const std::exception& error) {
    last_error = error.what();
    return NATIVE_METADATA_STATUS_INTERNAL_ERROR;
  } catch (...) {
    last_error = "unknown metadata failure";
    return NATIVE_METADATA_STATUS_INTERNAL_ERROR;
  }
}

std::string required(const char* value, const char* field) {
  if (value == nullptr) {
    throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                      std::string(field) + " pointer must not be null");
  }
  return value;
}

void set_info(const gsim::native::metadata::MetadataWriteResult& value,
              native_metadata_metadata_write_info* output) {
  if (output == nullptr) {
    throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                      "metadata write-info output must not be null");
  }
  output->record_count = value.record_count;
  output->bytes_written = value.bytes_written;
  output->maximum_record_bytes = value.maximum_record_bytes;
}

}  // namespace

extern "C" {

std::uint32_t native_metadata_abi_version(void) { return 1u; }

const char* native_metadata_library_version(void) { return "gsim-native-1"; }

const char* native_metadata_last_error(void) { return last_error.c_str(); }

native_metadata_status native_metadata_variant_metadata_create(
    const char* const* chromosome, const char* const* variant_id,
    const double* genetic_position_cm, const std::uint64_t* base_pair_position,
    const char* const* alternate_allele, const char* const* reference_allele,
    std::uint64_t count, native_metadata_variant_metadata_handle** output_handle) {
  if (output_handle != nullptr) *output_handle = nullptr;
  return protect([&] {
    if (output_handle == nullptr || chromosome == nullptr ||
        variant_id == nullptr || genetic_position_cm == nullptr ||
        base_pair_position == nullptr || alternate_allele == nullptr ||
        reference_allele == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "variant metadata arrays and output must not be null");
    }
    if (count > static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_extent,
                        "variant metadata count exceeds size_t");
    }
    std::vector<gsim::native::metadata::MarkerMetadata> records;
    records.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t index = 0; index < count; ++index) {
      records.push_back({required(chromosome[index], "chromosome"),
                         required(variant_id[index], "variant ID"),
                         genetic_position_cm[index], base_pair_position[index],
                         required(alternate_allele[index], "alternate allele"),
                         required(reference_allele[index], "reference allele")});
    }
    *output_handle = new native_metadata_variant_metadata_handle(std::move(records));
  });
}

native_metadata_status native_metadata_variant_metadata_read_bim(
    const char* source_utf8, native_metadata_variant_metadata_handle** output_handle) {
  if (output_handle != nullptr) *output_handle = nullptr;
  return protect([&] {
    if (output_handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "variant metadata output must not be null");
    }
    auto value = gsim::native::metadata::ValidatedVariantMetadata::read_bim(
        required(source_utf8, "BIM source"));
    *output_handle = new native_metadata_variant_metadata_handle(std::move(value));
  });
}

native_metadata_status native_metadata_variant_metadata_count(
    const native_metadata_variant_metadata_handle* handle, std::uint64_t* output_count) {
  return protect([&] {
    if (handle == nullptr || output_count == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "variant metadata handle/count output must not be null");
    }
    *output_count = handle->value.size();
  });
}

native_metadata_status native_metadata_variant_metadata_get(
    const native_metadata_variant_metadata_handle* handle, std::uint64_t index,
    const char** chromosome, const char** variant_id,
    double* genetic_position_cm, std::uint64_t* base_pair_position,
    const char** alternate_allele, const char** reference_allele) {
  return protect([&] {
    if (handle == nullptr || chromosome == nullptr || variant_id == nullptr ||
        genetic_position_cm == nullptr || base_pair_position == nullptr ||
        alternate_allele == nullptr || reference_allele == nullptr ||
        index >= handle->value.size()) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid variant metadata record request");
    }
    const auto& value = handle->value.records()[static_cast<std::size_t>(index)];
    *chromosome = value.chromosome.c_str();
    *variant_id = value.marker_id.c_str();
    *genetic_position_cm = value.genetic_distance_cm;
    *base_pair_position = value.base_pair_position;
    *alternate_allele = value.allele1.c_str();
    *reference_allele = value.allele2.c_str();
  });
}

native_metadata_status native_metadata_variant_metadata_close(native_metadata_variant_metadata_handle* handle) {
  return protect([&] {
    if (handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "variant metadata handle must not be null");
    }
    delete handle;
  });
}

native_metadata_status native_metadata_variant_metadata_write_bim(
    const native_metadata_variant_metadata_handle* handle, const char* destination_utf8,
    std::uint32_t overwrite, native_metadata_metadata_write_info* output_info) {
  return protect([&] {
    if (handle == nullptr || overwrite > 1u) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid BIM writer handle or overwrite flag");
    }
    set_info(handle->value.write_bim(
                 required(destination_utf8, "BIM destination"), overwrite == 1u),
             output_info);
  });
}

native_metadata_status native_metadata_sample_metadata_create(
    const char* const* family_id, const char* const* individual_id,
    const char* const* paternal_id, const char* const* maternal_id,
    const std::uint32_t* sex, std::uint64_t count,
    native_metadata_sample_metadata_handle** output_handle) {
  if (output_handle != nullptr) *output_handle = nullptr;
  return protect([&] {
    if (output_handle == nullptr || family_id == nullptr ||
        individual_id == nullptr || paternal_id == nullptr ||
        maternal_id == nullptr || sex == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "sample metadata arrays and output must not be null");
    }
    if (count > static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max())) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_extent,
                        "sample metadata count exceeds size_t");
    }
    std::vector<gsim::native::metadata::SimulationSampleMetadata> records;
    records.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t index = 0; index < count; ++index) {
      records.push_back({required(family_id[index], "family ID"),
                         required(individual_id[index], "individual ID"),
                         required(paternal_id[index], "paternal ID"),
                         required(maternal_id[index], "maternal ID"), sex[index]});
    }
    *output_handle = new native_metadata_sample_metadata_handle(std::move(records));
  });
}

native_metadata_status native_metadata_sample_metadata_read_fam(
    const char* source_utf8, native_metadata_sample_metadata_handle** output_handle) {
  if (output_handle != nullptr) *output_handle = nullptr;
  return protect([&] {
    if (output_handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "sample metadata output must not be null");
    }
    auto value = gsim::native::metadata::ValidatedSampleMetadata::read_fam(
        required(source_utf8, "FAM source"));
    *output_handle = new native_metadata_sample_metadata_handle(std::move(value));
  });
}

native_metadata_status native_metadata_sample_metadata_count(
    const native_metadata_sample_metadata_handle* handle, std::uint64_t* output_count) {
  return protect([&] {
    if (handle == nullptr || output_count == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "sample metadata handle/count output must not be null");
    }
    *output_count = handle->value.size();
  });
}

native_metadata_status native_metadata_sample_metadata_get(
    const native_metadata_sample_metadata_handle* handle, std::uint64_t index,
    const char** family_id, const char** individual_id,
    const char** paternal_id, const char** maternal_id, std::uint32_t* sex) {
  return protect([&] {
    if (handle == nullptr || family_id == nullptr || individual_id == nullptr ||
        paternal_id == nullptr || maternal_id == nullptr || sex == nullptr ||
        index >= handle->value.size()) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid sample metadata record request");
    }
    const auto& value = handle->value.records()[static_cast<std::size_t>(index)];
    *family_id = value.family_id.c_str();
    *individual_id = value.individual_id.c_str();
    *paternal_id = value.paternal_id.c_str();
    *maternal_id = value.maternal_id.c_str();
    *sex = value.sex;
  });
}

native_metadata_status native_metadata_sample_metadata_close(native_metadata_sample_metadata_handle* handle) {
  return protect([&] {
    if (handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "sample metadata handle must not be null");
    }
    delete handle;
  });
}

native_metadata_status native_metadata_sample_metadata_write_fam(
    const native_metadata_sample_metadata_handle* handle, const char* destination_utf8,
    std::uint32_t overwrite, native_metadata_metadata_write_info* output_info) {
  return protect([&] {
    if (handle == nullptr || overwrite > 1u) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid FAM writer handle or overwrite flag");
    }
    set_info(handle->value.write_fam(
                 required(destination_utf8, "FAM destination"), overwrite == 1u),
             output_info);
  });
}

native_metadata_status native_metadata_phased_vcf_reader_open(
    const char* source_utf8, native_metadata_phased_vcf_reader_handle** output_handle) {
  if (output_handle != nullptr) *output_handle = nullptr;
  return protect([&] {
    if (output_handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "VCF reader output must not be null");
    }
    *output_handle = new native_metadata_phased_vcf_reader_handle(
        required(source_utf8, "VCF source"));
  });
}

native_metadata_status native_metadata_phased_vcf_reader_close(native_metadata_phased_vcf_reader_handle* handle) {
  return protect([&] {
    if (handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "VCF reader handle must not be null");
    }
    delete handle;
  });
}

native_metadata_status native_metadata_phased_vcf_reader_dimensions(
    const native_metadata_phased_vcf_reader_handle* handle, std::uint64_t* sample_count,
    std::uint64_t* variant_count, std::uint64_t* chromosome_count) {
  return protect([&] {
    if (handle == nullptr || sample_count == nullptr || variant_count == nullptr ||
        chromosome_count == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "VCF reader dimensions arguments must not be null");
    }
    *sample_count = static_cast<std::uint64_t>(handle->value.samples().size());
    *variant_count = static_cast<std::uint64_t>(handle->value.variants().size());
    *chromosome_count = static_cast<std::uint64_t>(handle->value.chromosomes().size());
  });
}

native_metadata_status native_metadata_phased_vcf_reader_sample(
    const native_metadata_phased_vcf_reader_handle* handle, std::uint64_t index,
    const char** sample_id) {
  return protect([&] {
    if (handle == nullptr || sample_id == nullptr ||
        index >= handle->value.samples().size()) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid VCF sample request");
    }
    *sample_id = handle->value.samples()[static_cast<std::size_t>(index)].c_str();
  });
}

native_metadata_status native_metadata_phased_vcf_reader_variant(
    const native_metadata_phased_vcf_reader_handle* handle, std::uint64_t index,
    const char** chromosome, const char** variant_id,
    std::uint64_t* base_pair_position, const char** reference_allele,
    const char** alternate_allele, std::uint32_t* generated_id) {
  return protect([&] {
    if (handle == nullptr || chromosome == nullptr || variant_id == nullptr ||
        base_pair_position == nullptr || reference_allele == nullptr ||
        alternate_allele == nullptr || generated_id == nullptr ||
        index >= handle->value.variants().size()) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid VCF variant request");
    }
    const auto& value = handle->value.variants()[static_cast<std::size_t>(index)];
    *chromosome = value.chromosome.c_str();
    *variant_id = value.variant_id.c_str();
    *base_pair_position = value.base_pair_position;
    *reference_allele = value.reference_allele.c_str();
    *alternate_allele = value.alternate_allele.c_str();
    *generated_id = value.generated_id ? 1u : 0u;
  });
}

native_metadata_status native_metadata_phased_vcf_reader_chromosome(
    const native_metadata_phased_vcf_reader_handle* handle, std::uint64_t index,
    const char** chromosome, std::uint64_t* first_variant,
    std::uint64_t* variant_count) {
  return protect([&] {
    if (handle == nullptr || chromosome == nullptr || first_variant == nullptr ||
        variant_count == nullptr || index >= handle->value.chromosomes().size()) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid VCF chromosome request");
    }
    const auto& value = handle->value.chromosomes()[static_cast<std::size_t>(index)];
    *chromosome = value.chromosome.c_str();
    *first_variant = value.first_variant;
    *variant_count = value.variant_count;
  });
}

native_metadata_status native_metadata_phased_vcf_reader_start_chromosome(
    native_metadata_phased_vcf_reader_handle* handle, std::uint64_t chromosome_index) {
  return protect([&] {
    if (handle == nullptr) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "VCF reader handle must not be null");
    }
    handle->value.start_chromosome(chromosome_index);
  });
}

native_metadata_status native_metadata_phased_vcf_reader_next(
    native_metadata_phased_vcf_reader_handle* handle, std::uint8_t* h1, std::uint8_t* h2,
    std::uint64_t allele_capacity, std::uint64_t* variant_index,
    std::uint32_t* has_record) {
  return protect([&] {
    if (handle == nullptr || h1 == nullptr || h2 == nullptr ||
        variant_index == nullptr || has_record == nullptr ||
        allele_capacity != handle->value.samples().size()) {
      throw gsim::native::metadata::Error(gsim::native::metadata::StatusCode::invalid_argument,
                        "invalid VCF record buffer request");
    }
    std::vector<std::uint8_t> first;
    std::vector<std::uint8_t> second;
    if (!handle->value.next(*variant_index, first, second)) {
      *has_record = 0u;
      return;
    }
    std::copy(first.begin(), first.end(), h1);
    std::copy(second.begin(), second.end(), h2);
    *has_record = 1u;
  });
}

}  // extern "C"
