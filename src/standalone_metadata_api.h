#ifndef GSIM_STANDALONE_METADATA_API_HPP
#define GSIM_STANDALONE_METADATA_API_HPP

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum native_metadata_status {
  NATIVE_METADATA_STATUS_SUCCESS = 0,
  NATIVE_METADATA_STATUS_INVALID_ARGUMENT = 1,
  NATIVE_METADATA_STATUS_IO_ERROR = 2,
  NATIVE_METADATA_STATUS_INTERNAL_ERROR = 3
} native_metadata_status;

typedef struct native_metadata_variant_metadata_handle native_metadata_variant_metadata_handle;
typedef struct native_metadata_sample_metadata_handle native_metadata_sample_metadata_handle;
typedef struct native_metadata_phased_vcf_reader_handle native_metadata_phased_vcf_reader_handle;

typedef struct native_metadata_metadata_write_info {
  uint64_t record_count;
  uint64_t bytes_written;
  uint64_t maximum_record_bytes;
} native_metadata_metadata_write_info;

uint32_t native_metadata_abi_version(void);
const char* native_metadata_library_version(void);
const char* native_metadata_last_error(void);

native_metadata_status native_metadata_variant_metadata_create(
    const char* const* chromosome,
    const char* const* variant_id,
    const double* genetic_position_cm,
    const uint64_t* base_pair_position,
    const char* const* alternate_allele,
    const char* const* reference_allele,
    uint64_t count,
    native_metadata_variant_metadata_handle** output_handle);
native_metadata_status native_metadata_variant_metadata_read_bim(
    const char* source_utf8,
    native_metadata_variant_metadata_handle** output_handle);
native_metadata_status native_metadata_variant_metadata_count(
    const native_metadata_variant_metadata_handle* handle, uint64_t* output_count);
native_metadata_status native_metadata_variant_metadata_get(
    const native_metadata_variant_metadata_handle* handle, uint64_t index,
    const char** chromosome, const char** variant_id,
    double* genetic_position_cm, uint64_t* base_pair_position,
    const char** alternate_allele, const char** reference_allele);
native_metadata_status native_metadata_variant_metadata_close(
    native_metadata_variant_metadata_handle* handle);
native_metadata_status native_metadata_variant_metadata_write_bim(
    const native_metadata_variant_metadata_handle* handle,
    const char* destination_utf8,
    uint32_t overwrite,
    native_metadata_metadata_write_info* output_info);

native_metadata_status native_metadata_sample_metadata_create(
    const char* const* family_id,
    const char* const* individual_id,
    const char* const* paternal_id,
    const char* const* maternal_id,
    const uint32_t* sex,
    uint64_t count,
    native_metadata_sample_metadata_handle** output_handle);
native_metadata_status native_metadata_sample_metadata_read_fam(
    const char* source_utf8,
    native_metadata_sample_metadata_handle** output_handle);
native_metadata_status native_metadata_sample_metadata_count(
    const native_metadata_sample_metadata_handle* handle, uint64_t* output_count);
native_metadata_status native_metadata_sample_metadata_get(
    const native_metadata_sample_metadata_handle* handle, uint64_t index,
    const char** family_id, const char** individual_id,
    const char** paternal_id, const char** maternal_id, uint32_t* sex);
native_metadata_status native_metadata_sample_metadata_close(
    native_metadata_sample_metadata_handle* handle);
native_metadata_status native_metadata_sample_metadata_write_fam(
    const native_metadata_sample_metadata_handle* handle,
    const char* destination_utf8,
    uint32_t overwrite,
    native_metadata_metadata_write_info* output_info);

native_metadata_status native_metadata_phased_vcf_reader_open(
    const char* source_utf8, native_metadata_phased_vcf_reader_handle** output_handle);
native_metadata_status native_metadata_phased_vcf_reader_close(
    native_metadata_phased_vcf_reader_handle* handle);
native_metadata_status native_metadata_phased_vcf_reader_dimensions(
    const native_metadata_phased_vcf_reader_handle* handle, uint64_t* sample_count,
    uint64_t* variant_count, uint64_t* chromosome_count);
native_metadata_status native_metadata_phased_vcf_reader_sample(
    const native_metadata_phased_vcf_reader_handle* handle, uint64_t index,
    const char** sample_id);
native_metadata_status native_metadata_phased_vcf_reader_variant(
    const native_metadata_phased_vcf_reader_handle* handle, uint64_t index,
    const char** chromosome, const char** variant_id,
    uint64_t* base_pair_position, const char** reference_allele,
    const char** alternate_allele, uint32_t* generated_id);
native_metadata_status native_metadata_phased_vcf_reader_chromosome(
    const native_metadata_phased_vcf_reader_handle* handle, uint64_t index,
    const char** chromosome, uint64_t* first_variant,
    uint64_t* variant_count);
native_metadata_status native_metadata_phased_vcf_reader_start_chromosome(
    native_metadata_phased_vcf_reader_handle* handle, uint64_t chromosome_index);
native_metadata_status native_metadata_phased_vcf_reader_next(
    native_metadata_phased_vcf_reader_handle* handle, uint8_t* h1, uint8_t* h2,
    uint64_t allele_capacity, uint64_t* variant_index, uint32_t* has_record);

#ifdef __cplusplus
}
#endif

#endif
