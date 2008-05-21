#include <string.h>
#include "ShortRead.h"

static const int LINEBUF_SIZE = 20001;
static const int LINES_PER_FASTQ_REC = 4;
static const int LINES_PER_FASTA_REC = 2;

/*
 * Solexa 'fastq' files consist of records, each 4 lines long. Here is
 * an example:

 @HWI-EAS88_1_1_1_1001_499
 GGACTTTGTAGGATACCCTCGCTTTCCTTCTCCTGT
 +HWI-EAS88_1_1_1_1001_499
 ]]]]]]]]]]]]Y]Y]]]]]]]]]]]]VCHVMPLAS

 * inst/extdata/s_1_sequences.txt contains 256 records
 */

/*
 * Read a solexa 's_<lane>_sequence.txt' file into CharBBuf objects.
 */
static void
_read_solexa_fastq_file(const char *fname,
                        CharBBuf *seq, CharBBuf *name, CharBBuf *qualities)
{
    FILE *file;
    char linebuf[LINEBUF_SIZE];
    int lineno, reclineno, nchar_in_buf;

    if ((file = fopen(fname, "r")) == NULL)
        error("cannot open file %s", fname);

    lineno = 0;
    while (fgets(linebuf, LINEBUF_SIZE, file) != NULL) {
        if ((reclineno = lineno % LINES_PER_FASTQ_REC) == 2) {
            lineno++;
            continue;
        }

        nchar_in_buf = _rtrim(linebuf);
        if (nchar_in_buf >= LINEBUF_SIZE - 1) { // should never be >
            fclose(file);
            error("line too long %s:%d", fname, lineno);
        } else if (nchar_in_buf == 0) {
            fclose(file);
            error("unexpected empty line %s:%d", fname, lineno);
        }
        switch(reclineno) {
        case 0:
            /* add linebuf to CharBBuf; start at char +1 to skip the
             * fastq annotation. */
            append_string_to_CharBBuf(name, linebuf+1);
            break;
        case 1:
            _solexa_to_IUPAC(linebuf);
            append_string_to_CharBBuf(seq, linebuf);
            break;
        case 3:
            append_string_to_CharBBuf(qualities, linebuf);
            break;
        default:
            error("unexpected 'reclineno'; consult maintainer");
            break;
        }
        lineno++;
    }
    fclose(file);
    if ((lineno % LINES_PER_FASTQ_REC) != 0)
        error("unexpected number of lines in file '%s'", fname);
}

SEXP
read_solexa_fastq(SEXP files)
{
    CharBBuf seq, name, qualities;
    RoSeqs roSeqs;
    int i, nfiles, nrec = 0;
    const char *fname;
    SEXP nlines, ans = R_NilValue, nms=R_NilValue, elt;

    if (!IS_CHARACTER(files))
        Rf_error("'files' must be 'character'");

    /* Count total number of lines and pre-allocated space */
    nfiles = LENGTH(files);
    PROTECT(nlines = count_lines(files));
    for (i = 0; i < nfiles; ++i) 
        nrec += INTEGER(nlines)[i];
    UNPROTECT(1);
    nrec /= LINES_PER_FASTQ_REC;
    seq = new_CharBBuf(nrec, 0);
    name = new_CharBBuf(nrec, 0);
    qualities = new_CharBBuf(nrec, 0);

    for (i = 0; i < nfiles; ++i) {
        fname = translateChar(STRING_ELT(files, i));
        _read_solexa_fastq_file(fname, &seq, &name, &qualities);
    }

    PROTECT(ans = NEW_LIST(3));
    PROTECT(nms = NEW_CHARACTER(3));

    roSeqs = new_RoSeqs_from_BBuf(seq);
    PROTECT(elt = new_XStringSet_from_RoSeqs("DNAString", roSeqs));
    SET_VECTOR_ELT(ans, 0, elt);
    SET_STRING_ELT(nms, 0, mkChar("sread"));

    roSeqs = new_RoSeqs_from_BBuf(name);
    PROTECT(elt = new_XStringSet_from_RoSeqs("BString", roSeqs));
    SET_VECTOR_ELT(ans, 1, elt);
    SET_STRING_ELT(nms, 1, mkChar("id"));
    
    roSeqs = new_RoSeqs_from_BBuf(qualities);
    PROTECT(elt = new_XStringSet_from_RoSeqs("BString", roSeqs));
    SET_VECTOR_ELT(ans, 2, elt);
    SET_STRING_ELT(nms, 2, mkChar("quality"));

    setAttrib(ans, R_NamesSymbol, nms);

    UNPROTECT(5);
    return ans;
}

int
_io_XStringSet_columns(const char *fname, const int *colidx, const int ncol,
                       const char *sep, int header, CharBBuf *sets)
{
    FILE *file;
    char *linebuf;
    int lineno, nchar_in_buf;
    char *token;

    if ((file = fopen(fname, "r")) == NULL)
        error("cannot open file %s", fname);
    linebuf = S_alloc(LINEBUF_SIZE, sizeof(char)); /* auto free'd on return */

    /* header: first line ignored, errors ignored */
    if (header == TRUE)
        fgets(linebuf, LINEBUF_SIZE, file);
    lineno = 0;
    while (fgets(linebuf, LINEBUF_SIZE, file) != NULL) {
        nchar_in_buf = _rtrim(linebuf);
        if (nchar_in_buf >= LINEBUF_SIZE - 1) { // should never be >
            fclose(file);
            error("line too long %s:%d", fname, lineno);
        } else if (nchar_in_buf == 0) {
            fclose(file);
            error("unexpected empty line %s:%d", fname, lineno);
        }
        _solexa_to_IUPAC(linebuf);

        int j = 0, cidx=0;
        char *curbuf = linebuf;
        token = strsep(&curbuf, sep);
        for (j = 0; cidx < ncol && token != NULL; ++j) {
            if (j == colidx[cidx]) {
                append_string_to_CharBBuf(&sets[cidx], token);
                cidx++;
            }
            token = strsep(&curbuf, sep);
        } 
        lineno++;
    }
    fclose(file);
    return lineno;
}

SEXP
read_XStringSet_columns(SEXP files, SEXP colIndex, SEXP colClasses,
                        SEXP sep, SEXP header)
{
    if (!IS_CHARACTER(files))
        Rf_error("'files' must be 'character(1)'");
    if (!IS_INTEGER(colIndex) || LENGTH(colIndex) == 0)
        Rf_error("'colIndex' must be 'integer' with length > 0");
    if (!IS_CHARACTER(colClasses) || LENGTH(colClasses) != LENGTH(colIndex))
        Rf_error("'colClasses' must be 'character' with length(colClasses) == length(colIndex)");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1)
        Rf_error("'sep' must be character(1)");
    if (!IS_LOGICAL(header) || LENGTH(header) != 1)
        Rf_error("'header' must be logical(1)");

    /* Count lines and pre-allocate space */
    const char *csep = translateChar(STRING_ELT(sep, 0));
    const int nfiles = LENGTH(files);
    const int *nlines = INTEGER(count_lines(files));
    int nrow = 0;
    int i, j;
    for (i = 0; i < nfiles; ++i)
        nrow += nlines[i];
    nrow -= nfiles * LOGICAL(header)[0];
    int ncol = LENGTH(colIndex);
    CharBBuf *sets = (CharBBuf*) R_alloc(sizeof(CharBBuf), ncol);
    int *colidx = (int *) R_alloc(sizeof(int), ncol);
    for (j = 0; j < ncol; ++j) {
        sets[j] = new_CharBBuf(nrow, 0);
        colidx[j] = INTEGER(colIndex)[j] - 1;
    }

    /* read columns */
    int nreads = 0;
    for (i = 0; i < nfiles; ++i) {
        const char *fname = translateChar(STRING_ELT(files, i));
        nreads += _io_XStringSet_columns(fname, colidx, ncol,
                                         csep, LOGICAL(header)[0], sets);
    }
    if (nreads != nrow)
        Rf_error("found %d reads, expected %d", nlines, nrow);

    /* formulate return value */
    SEXP ans, elt;
    RoSeqs roSeqs;
    PROTECT(ans = NEW_LIST(ncol));
    for (j = 0; j < ncol; ++j) {
        roSeqs = new_RoSeqs_from_BBuf(sets[j]);
        const char *clsName = CHAR(STRING_ELT(colClasses, j));
        PROTECT(elt = new_XStringSet_from_RoSeqs(clsName, roSeqs));
        SET_VECTOR_ELT(ans, j, elt);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}
