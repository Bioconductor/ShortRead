#include "ShortRead.h"

static const R_CallMethodDef callMethods[] = {
    /* util.c */
    {".count_lines", (DL_FUNC) &count_lines, 1},
    /* io.c */
    {".read_solexa_fastq", (DL_FUNC) &read_solexa_fastq, 1},
    {".read_XStringSet_columns", (DL_FUNC) &read_XStringSet_columns, 5},
    /* alphabet */
    {".alphabet_by_cycle", (DL_FUNC) &alphabet_by_cycle, 3},
    {".alphabet_score", (DL_FUNC) &alphabet_score, 2},
    {NULL, NULL, 0}
};

void
R_init_ShortRead(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
