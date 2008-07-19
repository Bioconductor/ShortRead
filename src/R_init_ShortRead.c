#include "ShortRead.h"

static const R_CallMethodDef callMethods[] = {
    /* util.c */
    {".count_lines", (DL_FUNC) &count_lines, 1},
    /* io.c */
    {".read_prb_as_character", (DL_FUNC) &read_prb_as_character, 2},
    {".read_solexa_fastq", (DL_FUNC) &read_solexa_fastq, 1},
    {".read_XStringSet_columns", (DL_FUNC) &read_XStringSet_columns, 6},
    /* alphabet */
    {".alphabet_by_cycle", (DL_FUNC) &alphabet_by_cycle, 3},
    {".alphabet_pair_by_cycle", (DL_FUNC) &alphabet_pair_by_cycle, 5},
    {".alphabet_score", (DL_FUNC) &alphabet_score, 2},
    {".alphabet_as_int", (DL_FUNC) &alphabet_as_int, 2},
    {".alphabet_order", (DL_FUNC) &alphabet_order, 1},
    {".alphabet_duplicated", (DL_FUNC) &alphabet_duplicated, 1},
    {".alphabet_rank", (DL_FUNC) &alphabet_rank, 1},
    {".read_maq_map", (DL_FUNC) &read_maq_map, 2},
    /* pileup */
    {".pileup", (DL_FUNC) &pileup, 6},
    /* readBfaToc */
    {".readBfaToc", (DL_FUNC) &readBfaToc, 1},
    {NULL, NULL, 0}
};

void
R_init_ShortRead(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
