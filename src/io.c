#include <string.h>
#include <stdlib.h>             /* atoi */
#include "ShortRead.h"

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
 * Read a solexa .*_prb.txt file into STRING_VEC
 */

SEXP
read_prb_as_character(SEXP fname, SEXP cycles)
{
    const int MIN_SCORE = -40;
    const int NUC_PER_CYCLE = 4;
    const int SOLEXA_QBASE = 64;

    if (!IS_CHARACTER(fname) || LENGTH(fname)!=1)
        error("'fname' must be 'character(1)'");
    if (!IS_INTEGER(cycles) || LENGTH(cycles) != 1)
        error("'cycles' must be 'integer(1)'");
    const int n_reads = INTEGER(count_lines(fname))[0];
    const int n_cycles = INTEGER(cycles)[0];
    SEXP ans = PROTECT(NEW_CHARACTER(n_reads));

    FILE *file = _fopen(CHAR(STRING_ELT(fname, 0)), "r");

    int read=0, cycle=0, nuc=0, curr_max=MIN_SCORE, val;
    char *score = R_alloc(sizeof(char), n_cycles + 1);
    score[n_cycles] = '\0';
    while (fscanf(file, "%d", &val)==1) {
        if (val > curr_max) curr_max = val;
        if (++nuc == NUC_PER_CYCLE) {
            score[cycle] = (char) curr_max + SOLEXA_QBASE;
            nuc = 0;
            curr_max=MIN_SCORE;
            if (++cycle == n_cycles) {
                if (read >= n_reads)
                    error("unexpected number of reads: %d instead of %d",
                          read, n_reads);
                SET_STRING_ELT(ans, read++, mkChar(score));
                cycle = 0;
            }
        }
    }
    UNPROTECT(1);
    return ans;
}

/*
 * Read a solexa 's_<lane>_sequence.txt' file into CharAEAE objects.
 */
static void
_read_solexa_fastq_file(const char *fname,
                        CharAEAE *seq, CharAEAE *name, CharAEAE *qualities)
{
    FILE *file;
    char linebuf[LINEBUF_SIZE];
    int lineno, reclineno, nchar_in_buf;

    file = _fopen(fname, "r");
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
            /* add linebuf to CharAEAE; start at char +1 to skip the
             * fastq annotation. */
            append_string_to_CharAEAE(name, linebuf+1);
            break;
        case 1:
            _solexa_to_IUPAC(linebuf);
            append_string_to_CharAEAE(seq, linebuf);
            break;
        case 3:
            append_string_to_CharAEAE(qualities, linebuf);
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
    CharAEAE seq, name, qualities;
    RoSeqs roSeqs;
    int i, nfiles, nrec = 0;
    const char *fname;
    SEXP ans = R_NilValue, nms = R_NilValue, elt;

    if (!IS_CHARACTER(files))
        Rf_error("'files' must be 'character'");

    /* pre-allocated space */
    nfiles = LENGTH(files);
    nrec = _count_lines_sum(files) / LINES_PER_FASTQ_REC;
    seq = new_CharAEAE(nrec, 0);
    name = new_CharAEAE(nrec, 0);
    qualities = new_CharAEAE(nrec, 0);

    for (i = 0; i < nfiles; ++i) {
        R_CheckUserInterrupt();
        fname = translateChar(STRING_ELT(files, i));
        _read_solexa_fastq_file(fname, &seq, &name, &qualities);
    }

    PROTECT(ans = NEW_LIST(3));
    PROTECT(nms = NEW_CHARACTER(3));

    roSeqs = new_RoSeqs_from_CharAEAE(&seq);
    PROTECT(elt = new_XStringSet_from_RoSeqs("DNAString", &roSeqs));
    SET_VECTOR_ELT(ans, 0, elt);
    SET_STRING_ELT(nms, 0, mkChar("sread"));

    roSeqs = new_RoSeqs_from_CharAEAE(&name);
    PROTECT(elt = new_XStringSet_from_RoSeqs("BString", &roSeqs));
    SET_VECTOR_ELT(ans, 1, elt);
    SET_STRING_ELT(nms, 1, mkChar("id"));
    
    roSeqs = new_RoSeqs_from_CharAEAE(&qualities);
    PROTECT(elt = new_XStringSet_from_RoSeqs("BString", &roSeqs));
    SET_VECTOR_ELT(ans, 2, elt);
    SET_STRING_ELT(nms, 2, mkChar("quality"));

    setAttrib(ans, R_NamesSymbol, nms);

    UNPROTECT(5);
    return ans;
}

int
_io_XStringSet_columns(const char *fname, const int *colidx, int ncol,
                       const char *sep, MARK_FIELD_FUNC *mark_field,
                       int header, const char *commentChar, 
                       CharAEAE *sets, const int *toIUPAC)
{
    FILE *file;
    char *linebuf;
    int lineno;

    file = _fopen(fname, "r");
    linebuf = S_alloc(LINEBUF_SIZE, sizeof(char)); /* auto free'd on return */

    /* header: first line ignored, errors ignored */
    if (header == TRUE)
        fgets(linebuf, LINEBUF_SIZE, file);
    lineno = 0;
    while (fgets(linebuf, LINEBUF_SIZE, file) != NULL) {
        if (_linebuf_skip_p(linebuf, file, fname,
                            lineno, commentChar)) {
            lineno++;
            continue;
        }

        int j = 0, cidx=0;
        char *curr = linebuf, *next;
        for (j = 0; cidx < ncol && curr != NULL; ++j) {
            next = (*mark_field)(curr, sep);
            if (j == colidx[cidx]) {
                if (toIUPAC[cidx])
                    _solexa_to_IUPAC(curr);
                append_string_to_CharAEAE(&sets[cidx], curr);
                cidx++;
            }
            curr = next;
        } 
        lineno++;
    }
    fclose(file);
    return lineno;
}

SEXP
read_XStringSet_columns(SEXP files, SEXP colIndex, SEXP colClasses,
                        SEXP sep, SEXP header, SEXP commentChar)
{
    if (!IS_CHARACTER(files))
        Rf_error("'files' must be 'character(1)'");
    if (!IS_INTEGER(colIndex) || LENGTH(colIndex) == 0)
        Rf_error("'colIndex' must be 'integer' with length > 0");
    if (!IS_CHARACTER(colClasses) || LENGTH(colClasses) != LENGTH(colIndex))
        Rf_error("'colClasses' must be 'character' with length(colClasses) == length(colIndex)");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1)
        Rf_error("'sep' must be character(1)"); 
    /* FIXME: !nzchar(sep[1]) */
    if (!IS_LOGICAL(header) || LENGTH(header) != 1)
        Rf_error("'header' must be logical(1)");
    if (!IS_CHARACTER(commentChar) || LENGTH(commentChar) != 1)
        Rf_error("'commentChar' must be character(1)");
    if (LENGTH(STRING_ELT(commentChar, 0)) != 1)
        Rf_error("'nchar(commentChar[[1]])' must be 1 but is %d",
                 LENGTH(STRING_ELT(commentChar, 0)));

    /* Count lines and pre-allocate space */
    const char *csep = translateChar(STRING_ELT(sep, 0));
    const int nfiles = LENGTH(files);
    const int *nlines = INTEGER(count_lines(files));
    MARK_FIELD_FUNC *sep_func;  /* how to parse fields; minor efficiency */
    if (csep[0] != '\0' && csep[1] == '\0')
        sep_func = _mark_field_1;
    else
        sep_func = _mark_field_n;
    int nrow = 0;
    int i, j;
    for (i = 0; i < nfiles; ++i)
        nrow += nlines[i];
    nrow -= nfiles * LOGICAL(header)[0];
    int ncol = LENGTH(colIndex);
    CharAEAE *sets = (CharAEAE*) R_alloc(sizeof(CharAEAE), ncol);
    int *colidx = (int *) R_alloc(sizeof(int), ncol);
    int *toIUPAC = (int *) R_alloc(sizeof(int), ncol);
    for (j = 0; j < ncol; ++j) {
        sets[j] = new_CharAEAE(nrow, 0);
        colidx[j] = INTEGER(colIndex)[j] - 1;
        toIUPAC[j] = !strcmp(CHAR(STRING_ELT(colClasses, j)), "DNAString");
    }

    /* read columns */
    int nreads = 0;
    for (i = 0; i < nfiles; ++i) {
        const char *fname = translateChar(STRING_ELT(files, i));
        R_CheckUserInterrupt();
        nreads += _io_XStringSet_columns(fname, colidx, ncol,
                                         csep, sep_func,
                                         LOGICAL(header)[0],
                                         CHAR(STRING_ELT(commentChar, 0)),
                                         sets, toIUPAC);
    }
    if (nreads != nrow)
        Rf_error("found %d reads, expected %d", nreads, nrow);

    /* formulate return value */
    SEXP ans, elt;
    PROTECT(ans = NEW_LIST(ncol));
    for (j = 0; j < ncol; ++j) {
        const char *clsName = CHAR(STRING_ELT(colClasses, j));
        PROTECT(elt = _CharAEAE_to_XStringSet(sets + j, clsName));
        SET_VECTOR_ELT(ans, j, elt);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}

/*
 * _export parser
 */

typedef struct {
    SEXP ref;
    int offset;
    SEXP run;
    int *lane,
        *tile,
        *x,
        *y;
    CharAEAE read, quality;
    SEXP chromosome;
    int *position,
        *strand,
        *alignQuality,
        *filtering;
} SOLEXA_EXPORT_REC;

static const char *ELT_NMS[] = {
    "run", "lane", "tile", "x", "y", "sread", "quality",
    "chromosome", "position", "strand", "alignQuality",
    "filtering"
};
static const int N_ELTS = sizeof(ELT_NMS) / sizeof(const char*);

void 
_solexa_export_rec_set_offset(SOLEXA_EXPORT_REC *rec, int offset);

SOLEXA_EXPORT_REC *
_solexa_export_rec_new(SEXP ref, int nrec)
{
    if (LENGTH(ref) != N_ELTS)
        Rf_error("_solexa_export_rec_set_offset internal error: LENGTH(ref) != N_ELTS)");

    SEXP names = PROTECT(NEW_CHARACTER(N_ELTS));
    SET_VECTOR_ELT(ref, 0, NEW_STRING(nrec)); /* run */
    SET_VECTOR_ELT(ref, 1, NEW_INTEGER(nrec)); /* lane */
    SET_VECTOR_ELT(ref, 2, NEW_INTEGER(nrec)); /* tile */
    SET_VECTOR_ELT(ref, 3, NEW_INTEGER(nrec)); /* x */
    SET_VECTOR_ELT(ref, 4, NEW_INTEGER(nrec)); /* y */
    SET_VECTOR_ELT(ref, 7, NEW_STRING(nrec));  /* chromosome */
    SET_VECTOR_ELT(ref, 8, NEW_INTEGER(nrec)); /* position */
    SET_VECTOR_ELT(ref, 9, NEW_INTEGER(nrec)); /* strand: factor */
    SET_VECTOR_ELT(ref, 10, NEW_INTEGER(nrec)); /* alignQuality */
    SET_VECTOR_ELT(ref, 11, NEW_INTEGER(nrec)); /* filtering: factor */
    for (int i = 0; i < N_ELTS; ++i)
        SET_STRING_ELT(names, i, mkChar(ELT_NMS[i]));
    SET_ATTR(ref, R_NamesSymbol, names);

    SOLEXA_EXPORT_REC *rec = 
        (SOLEXA_EXPORT_REC *) R_alloc(1, sizeof(SOLEXA_EXPORT_REC));
    rec->ref = ref;
    rec->read = new_CharAEAE(nrec, 0);
    rec->quality = new_CharAEAE(nrec, 0);

    _solexa_export_rec_set_offset(rec, 0);

    UNPROTECT(1);
    return rec;
}

/*
 * view-like; no allocation
 */
void
_solexa_export_rec_set_offset(SOLEXA_EXPORT_REC *rec, int offset)
{
    SEXP ref = rec->ref;
    rec->offset = offset;
    rec->run = VECTOR_ELT(ref, 0);
    rec->lane = INTEGER(VECTOR_ELT(ref, 1)) + offset;
    rec->tile = INTEGER(VECTOR_ELT(ref, 2)) + offset;
    rec->x = INTEGER(VECTOR_ELT(ref, 3)) + offset;
    rec->y = INTEGER(VECTOR_ELT(ref, 4)) + offset;
    rec->chromosome = VECTOR_ELT(ref, 7);
    rec->position = INTEGER(VECTOR_ELT(ref, 8)) + offset;
    rec->strand = INTEGER(VECTOR_ELT(ref, 9)) + offset;
    rec->alignQuality = INTEGER(VECTOR_ELT(ref, 10)) + offset;
    rec->filtering = INTEGER(VECTOR_ELT(ref, 11)) + offset;
}

#define NEW_CALL(S, T, NAME, ENV, N) \
    PROTECT(S = T = allocList(N)); \
    SET_TYPEOF(T, LANGSXP); \
    SETCAR(T, findFun(install(NAME), ENV)); \
    T = CDR(T)
#define CSET_CDR(T, NAME, VALUE) \
    SETCAR(T, VALUE); \
    SET_TAG(T, install(NAME)); \
    T = CDR(T)
#define CEVAL_TO(S, ENV, GETS) \
    GETS = eval(S, ENV); \
    UNPROTECT(1)

SEXP
_AlignedRead_Solexa_make(SOLEXA_EXPORT_REC *rec)
{
    const char *FILTER_LEVELS[] = { "Y", "N" };
    SEXP ref = rec->ref;
    SEXP s, t, nmspc = PROTECT(_get_namespace("ShortRead"));

    SEXP sfq;                   /* SFastqQuality(rec->quality]) */
    SEXP quality;
    PROTECT(quality = _CharAEAE_to_XStringSet(&(rec->quality), 
                                              "BString"));
    NEW_CALL(s, t, "SFastqQuality", nmspc, 2);
    CSET_CDR(t, "quality", quality);
    CEVAL_TO(s, nmspc, sfq);
    PROTECT(sfq);

    SEXP alnq;                  /* NumericQuality() */
    NEW_CALL(s, t, "NumericQuality", nmspc, 2);
    CSET_CDR(t, "quality", VECTOR_ELT(ref, 10));
    CEVAL_TO(s, nmspc, alnq);
    PROTECT(alnq);

    SEXP adf;    /* .readAligned_SolexaExport_AlignedDataFrame(...) */
    _as_factor(VECTOR_ELT(ref, 11), FILTER_LEVELS,
               sizeof(FILTER_LEVELS) / sizeof(const char *));
    NEW_CALL(s, t, ".SolexaExport_AlignedDataFrame", nmspc, 7);
    CSET_CDR(t, "run", VECTOR_ELT(ref, 0)); 
    CSET_CDR(t, "lane", VECTOR_ELT(ref, 1)); 
    CSET_CDR(t, "tile", VECTOR_ELT(ref, 2)); 
    CSET_CDR(t, "x", VECTOR_ELT(ref, 3));
    CSET_CDR(t, "y", VECTOR_ELT(ref, 4));
    CSET_CDR(t, "filtering", VECTOR_ELT(ref, 11));
    CEVAL_TO(s, nmspc, adf);
    PROTECT(adf);

    SEXP aln;
    SEXP sread;
    PROTECT(sread = _CharAEAE_to_XStringSet(&(rec->read), 
                                            "DNAString"));
    SEXP strand_lvls = PROTECT(_get_strand_levels());
    _as_factor_SEXP(VECTOR_ELT(ref, 9), strand_lvls);
    NEW_CALL(s, t, "AlignedRead", nmspc, 8);
    CSET_CDR(t, "sread", sread);
    CSET_CDR(t, "quality", sfq); 
    CSET_CDR(t, "chromosome", VECTOR_ELT(ref, 7));
    CSET_CDR(t, "position", VECTOR_ELT(ref, 8));
    CSET_CDR(t, "strand", VECTOR_ELT(ref, 9)); 
    CSET_CDR(t, "alignQuality", alnq);
    CSET_CDR(t, "alignData", adf);
    CEVAL_TO(s, nmspc, aln);

    UNPROTECT(7);
    return aln;
}

#undef NEW_CALL
#undef CSET_CDR
#undef CEVAL_TO

int
_read_solexa_export_file(const char *fname, const char *csep,
                         const char *commentChar,
                         MARK_FIELD_FUNC *mark_func,
                         SOLEXA_EXPORT_REC *rec)
{
    const int N_FIELDS = 22;
    FILE *file;
    char linebuf[LINEBUF_SIZE], *elt[N_FIELDS];
    int lineno = 0, nrec = 0, i;

    file = _fopen(fname, "r");
    while (fgets(linebuf, LINEBUF_SIZE, file) != NULL) {
        if (_linebuf_skip_p(linebuf, file,
                            fname, lineno, commentChar)) {
            lineno++;
            continue;
        }

        /* field-ify */
        elt[0] = linebuf;
        for (i = 1; i < N_FIELDS; ++i) {
            elt[i] = (*mark_func)(elt[i-1], csep);
            if (elt[i] == elt[i-1])
                error("too few fields, %s:%d", fname, lineno);
        }
            
        SET_STRING_ELT(rec->run,
                       rec->offset + nrec, mkChar(elt[1]));
        rec->lane[lineno] = atoi(elt[2]);
        rec->tile[lineno] = atoi(elt[3]);
        rec->x[lineno] = atoi(elt[4]);
        rec->y[lineno] = atoi(elt[5]);
        /* 6: indexString, 7: pairedReadNumber */
        append_string_to_CharAEAE(&(rec->read), elt[8]);
        append_string_to_CharAEAE(&(rec->quality), elt[9]);
        SET_STRING_ELT(rec->chromosome,
                       rec->offset + nrec, mkChar(elt[10]));
        /* 11: contig */
        if (*elt[12] == '\0')
            rec->position[lineno] = NA_INTEGER;
        else
            rec->position[lineno] = atoi(elt[12]);
        if (*elt[13] == '\0')
            rec->strand[lineno] = NA_INTEGER;
        else {
            switch(*elt[13]) {
            case 'R':
                rec->strand[lineno] = 1;
                break;
            case 'F':
                rec->strand[lineno] = 2;
                break;
            default:
                error("invalid 'strand' field '%s', %s:%d",
                      *elt[13], fname, lineno);
                break;
            }
        }
        /* 14: descriptor */
        rec->alignQuality[lineno] = atoi(elt[15]);
        /* 16: pairedScore, 17: partnerCzome, 18: partnerContig
           19: partnerOffset, 20: partnerStrand */
        switch (*elt[21]) {
        case 'Y':
            rec->filtering[lineno] = 1;
            break;
        case 'N':
            rec->filtering[lineno] = 2;
            break;
        default:
            error("invalid 'filtering' field '%s', %s:%d",
                  *elt[21], fname, lineno);
            break;
        }
        lineno++;
        nrec++;
    }
    
    return nrec;
}

SEXP
read_solexa_export(SEXP files, SEXP sep, SEXP commentChar)
{
    if (!IS_CHARACTER(files))
        Rf_error("'files' must be 'character()'");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1)
        Rf_error("'sep' must be character(1)"); 
    /* FIXME: !nzchar(sep[1]) */
    if (!IS_CHARACTER(commentChar) || LENGTH(commentChar) != 1)
        Rf_error("'commentChar' must be character(1)");
    if (LENGTH(STRING_ELT(commentChar, 0)) != 1)
        Rf_error("'nchar(commentChar[[1]])' must be 1 but is %d",
                 LENGTH(STRING_ELT(commentChar, 0)));

    int nrec = _count_lines_sum(files);
    SEXP ref = PROTECT(NEW_LIST(N_ELTS));;
    SOLEXA_EXPORT_REC *rec = _solexa_export_rec_new(ref, nrec);

    const char *csep = translateChar(STRING_ELT(sep, 0));
    MARK_FIELD_FUNC *sep_func;/* how to parse fields; minor efficiency */
    if (csep[0] != '\0' && csep[1] == '\0')
        sep_func = _mark_field_1;
    else
        sep_func = _mark_field_n;

    nrec = 0;
    for (int i = 0; i < LENGTH(files); ++i) {
        R_CheckUserInterrupt();
        _solexa_export_rec_set_offset(rec, nrec);
        nrec += _read_solexa_export_file(
            CHAR(STRING_ELT(files, i)), csep,
            CHAR(STRING_ELT(commentChar, 0)),
            sep_func, rec);
    }
    SEXP aln = _AlignedRead_Solexa_make(rec);
    UNPROTECT(1);
    return aln;
}
