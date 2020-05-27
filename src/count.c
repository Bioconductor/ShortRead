#include "ShortRead.h"

/*
 * Count the number of lines ('\n') in a file.
 *
 * file: an open file stream at position 0
 *
 */
static double _count_lines(gzFile file)
{
    const int LINEBUF_SIZE = 20001;
    size_t bytes_read;
    char *buf = (char *) R_alloc(LINEBUF_SIZE + 1, sizeof(char));
    double lines = 0;

    while ((bytes_read = gzread(file, buf, LINEBUF_SIZE)) > 0) {
        char *p = buf;
        while ((p = memchr(p, '\n', (buf + bytes_read) - p))) {
            ++p;
            ++lines;
        }
    }
    return lines;
}

double _count_lines_sum(SEXP files)
{
    SEXP nlines = PROTECT(count_lines(files));
    int i;
    double nrec = 0;
    for (i = 0; i < LENGTH(files); ++i)
        nrec += REAL(nlines)[i];

    UNPROTECT(1);
    return nrec;
}

SEXP count_lines(SEXP files)
{
    int i, nfile;
    const char *filepath;
    gzFile file;
    SEXP ans = R_NilValue;

    if (!IS_CHARACTER(files))
        error("'files' must be character()");

    nfile = LENGTH(files);
    PROTECT(ans = NEW_NUMERIC(nfile));
    for (i = 0; i < nfile; ++i) {
        R_CheckUserInterrupt();
        filepath = translateChar(STRING_ELT(files, i));
        file = _fopen(filepath, "rb");
        REAL(ans)[i] = _count_lines(file);
        gzclose(file);
    }

    UNPROTECT(1);
    return ans;
}

/*
 * Count the number of records, nucleotides & quality scores. Some
 * validity checking implemented in `kseq_read()`
 */

#include <stdint.h>
#include "htslib/kseq.h"
KSEQ_INIT(gzFile, gzread)

void _count_records(const char *fname, int *result)
{
    gzFile fp;
    kseq_t *seq;
    int r, n = 0, slen = 0, qlen = 0;
    fp = gzopen(fname, "r");
    seq = kseq_init(fp);
    while ((r = kseq_read(seq)) >= 0)
        ++n, slen += seq->seq.l, qlen += seq->qual.l;
    if (r != -1) Rf_error("malformed FASTQ in record '%d'", n + 1);
    kseq_destroy(seq);
    gzclose(fp);

    result[0] = n;
    result[1] = slen;
    result[2] = qlen;
}

SEXP count_records(SEXP filename)
{
    SEXP result = PROTECT(Rf_allocVector(INTSXP, 3));

    _count_records(Rf_translateChar(Rf_asChar(filename)), INTEGER(result));

    UNPROTECT(1);
    return result;
}
