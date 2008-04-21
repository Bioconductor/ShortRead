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
