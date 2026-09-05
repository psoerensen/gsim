#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_gsim_hapnest_founders(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                     SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                     SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_plan(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_chromosome_keys(SEXP);
extern SEXP C_gsim_hapnest_segment_endpoint(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_copy_segment(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_hapnest_pair(SEXP, SEXP);
extern SEXP C_gsim_meiosis_materialize(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_meiosis_draw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                SEXP);
extern SEXP C_gsim_meiosis_plan(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_backend(SEXP);
extern SEXP C_gsim_gbits_pack(SEXP, SEXP);
extern SEXP C_gsim_gbits_zero(SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_unpack(SEXP);
extern SEXP C_gsim_gbits_info(SEXP);
extern SEXP C_gsim_gbits_close(SEXP);
extern SEXP C_gsim_gbits_word(SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_copy_interval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_copy_filtered(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                       SEXP, SEXP);
extern SEXP C_gsim_gbits_copy_filtered_counts(SEXP, SEXP, SEXP, SEXP, SEXP,
                                              SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_make_gamete(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                     SEXP);
extern SEXP C_gsim_gbits_decode_genotypes(SEXP, SEXP);
extern SEXP C_gsim_gbits_bed_sink_create(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_bed_sink_append(SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_bed_sink_finalize(SEXP);
extern SEXP C_gsim_gbits_bed_sink_cancel(SEXP);
extern SEXP C_gsim_gbits_bed_sink_info(SEXP);
extern SEXP C_gsim_gbits_hap_sink_create(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_hap_sink_append(SEXP, SEXP, SEXP);
extern SEXP C_gsim_gbits_hap_sink_finalize(SEXP);
extern SEXP C_gsim_gbits_hap_sink_cancel(SEXP);
extern SEXP C_gsim_gbits_hap_sink_info(SEXP);
extern SEXP C_gsim_gbits_hap_reader_open(SEXP, SEXP);
extern SEXP C_gsim_gbits_hap_reader_close(SEXP);
extern SEXP C_gsim_gbits_hap_reader_info(SEXP);
extern SEXP C_gsim_gbits_hap_reader_load(SEXP, SEXP);
extern SEXP C_gsim_gbits_bed_read_all(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_gmat_backend(SEXP);
extern SEXP C_gsim_gmat_variant_create(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                                        SEXP);
extern SEXP C_gsim_gmat_sample_create(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_gsim_gmat_write_bim(SEXP, SEXP);
extern SEXP C_gsim_gmat_write_fam(SEXP, SEXP);
extern SEXP C_gsim_gmat_read_bim(SEXP, SEXP);
extern SEXP C_gsim_gmat_read_fam(SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"C_gsim_hapnest_founders", (DL_FUNC) &C_gsim_hapnest_founders, 17},
    {"C_gsim_hapnest_plan", (DL_FUNC) &C_gsim_hapnest_plan, 13},
    {"C_gsim_hapnest_chromosome_keys", (DL_FUNC) &C_gsim_hapnest_chromosome_keys, 1},
    {"C_gsim_hapnest_segment_endpoint", (DL_FUNC) &C_gsim_hapnest_segment_endpoint, 4},
    {"C_gsim_hapnest_copy_segment", (DL_FUNC) &C_gsim_hapnest_copy_segment, 5},
    {"C_gsim_hapnest_pair", (DL_FUNC) &C_gsim_hapnest_pair, 2},
    {"C_gsim_meiosis_materialize", (DL_FUNC) &C_gsim_meiosis_materialize, 5},
    {"C_gsim_meiosis_draw", (DL_FUNC) &C_gsim_meiosis_draw, 8},
    {"C_gsim_meiosis_plan", (DL_FUNC) &C_gsim_meiosis_plan, 5},
    {"C_gsim_gbits_backend", (DL_FUNC) &C_gsim_gbits_backend, 1},
    {"C_gsim_gbits_pack", (DL_FUNC) &C_gsim_gbits_pack, 2},
    {"C_gsim_gbits_zero", (DL_FUNC) &C_gsim_gbits_zero, 3},
    {"C_gsim_gbits_unpack", (DL_FUNC) &C_gsim_gbits_unpack, 1},
    {"C_gsim_gbits_info", (DL_FUNC) &C_gsim_gbits_info, 1},
    {"C_gsim_gbits_close", (DL_FUNC) &C_gsim_gbits_close, 1},
    {"C_gsim_gbits_word", (DL_FUNC) &C_gsim_gbits_word, 3},
    {"C_gsim_gbits_copy_interval", (DL_FUNC) &C_gsim_gbits_copy_interval, 6},
    {"C_gsim_gbits_copy_filtered", (DL_FUNC) &C_gsim_gbits_copy_filtered, 8},
    {"C_gsim_gbits_copy_filtered_counts", (DL_FUNC) &C_gsim_gbits_copy_filtered_counts, 8},
    {"C_gsim_gbits_make_gamete", (DL_FUNC) &C_gsim_gbits_make_gamete, 7},
    {"C_gsim_gbits_decode_genotypes", (DL_FUNC) &C_gsim_gbits_decode_genotypes, 2},
    {"C_gsim_gbits_bed_sink_create", (DL_FUNC) &C_gsim_gbits_bed_sink_create, 5},
    {"C_gsim_gbits_bed_sink_append", (DL_FUNC) &C_gsim_gbits_bed_sink_append, 3},
    {"C_gsim_gbits_bed_sink_finalize", (DL_FUNC) &C_gsim_gbits_bed_sink_finalize, 1},
    {"C_gsim_gbits_bed_sink_cancel", (DL_FUNC) &C_gsim_gbits_bed_sink_cancel, 1},
    {"C_gsim_gbits_bed_sink_info", (DL_FUNC) &C_gsim_gbits_bed_sink_info, 1},
    {"C_gsim_gbits_hap_sink_create", (DL_FUNC) &C_gsim_gbits_hap_sink_create, 4},
    {"C_gsim_gbits_hap_sink_append", (DL_FUNC) &C_gsim_gbits_hap_sink_append, 3},
    {"C_gsim_gbits_hap_sink_finalize", (DL_FUNC) &C_gsim_gbits_hap_sink_finalize, 1},
    {"C_gsim_gbits_hap_sink_cancel", (DL_FUNC) &C_gsim_gbits_hap_sink_cancel, 1},
    {"C_gsim_gbits_hap_sink_info", (DL_FUNC) &C_gsim_gbits_hap_sink_info, 1},
    {"C_gsim_gbits_hap_reader_open", (DL_FUNC) &C_gsim_gbits_hap_reader_open, 2},
    {"C_gsim_gbits_hap_reader_close", (DL_FUNC) &C_gsim_gbits_hap_reader_close, 1},
    {"C_gsim_gbits_hap_reader_info", (DL_FUNC) &C_gsim_gbits_hap_reader_info, 1},
    {"C_gsim_gbits_hap_reader_load", (DL_FUNC) &C_gsim_gbits_hap_reader_load, 2},
    {"C_gsim_gbits_bed_read_all", (DL_FUNC) &C_gsim_gbits_bed_read_all, 4},
    {"C_gsim_gmat_backend", (DL_FUNC) &C_gsim_gmat_backend, 1},
    {"C_gsim_gmat_variant_create", (DL_FUNC) &C_gsim_gmat_variant_create, 7},
    {"C_gsim_gmat_sample_create", (DL_FUNC) &C_gsim_gmat_sample_create, 6},
    {"C_gsim_gmat_write_bim", (DL_FUNC) &C_gsim_gmat_write_bim, 2},
    {"C_gsim_gmat_write_fam", (DL_FUNC) &C_gsim_gmat_write_fam, 2},
    {"C_gsim_gmat_read_bim", (DL_FUNC) &C_gsim_gmat_read_bim, 2},
    {"C_gsim_gmat_read_fam", (DL_FUNC) &C_gsim_gmat_read_fam, 2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_gsim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
