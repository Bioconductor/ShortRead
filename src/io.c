#include <string.h>
#include <stdlib.h>             /* atoi */
#include <zlib.h>
#include "ShortRead.h"
#include "call.h"

static const int SOLEXA_QBASE = 64;
static const int PHRED_QBASE = 33;

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

void _write_err(int i)
{
    Rf_error("failed to write record %d", i + 1);
}

char *_holder_to_char(XStringSet_holder * holder, const int i,
                     char *buf, const int width, DECODE_FUNC decode)
{
    Chars_holder chars_holder = get_elt_from_XStringSet_holder(holder, i);
    if (chars_holder.length > width)
        return NULL;
    if (decode != NULL) {
        int j;
        for (j = 0; j < chars_holder.length; ++j)
            buf[j] = decode(chars_holder.seq[j]);
    } else
        strncpy(buf, chars_holder.seq, chars_holder.length);
    buf[chars_holder.length] = '\0';
    return buf;
}

SEXP write_fastq(SEXP id, SEXP sread, SEXP quality,
                 SEXP fname, SEXP fmode, SEXP full, SEXP compress,
                 SEXP max_width)
{
    if (!(IS_S4_OBJECT(id) && strcmp(get_classname(id), "BStringSet") == 0))
        Rf_error("'%s' must be '%s'", "id", "BStringSet");
    if (!(IS_S4_OBJECT(sread) &&
          strcmp(get_classname(sread), "DNAStringSet") == 0))
        Rf_error("'%s' must be '%s'", "sread", "DNAStringSet");
    /* check in R -- C-level R_check_super... is not adequate */
    /* if (R_check_class_etc(quality, qualityClasses) < 0) */
    /*     Rf_error("'is(<%s>, \"%s\")' failed", "quality", qualityClasses[0]); */
    const int len = get_XStringSet_length(id);
    if ((len != get_XStringSet_length(sread)) ||
        (len != get_XStringSet_length(quality)))
        Rf_error("length() of %s must all be equal",
                 "'id', 'sread', 'quality'");
    if (!(IS_CHARACTER(fname) && LENGTH(fname) == 1))	/* FIXME: nzchar */
        Rf_error("'%s' must be '%s'", "file", "character(1)");
    if (!(IS_CHARACTER(fmode) && LENGTH(fmode) == 1))	/* FIXME nchar()<3 */
        Rf_error("'%s' must be '%s'", "mode", "character(1)");
    if (!(IS_LOGICAL(full) && LENGTH(full) == 1))
        Rf_error("'%s' must be '%s'", "full", "logical(1)");
    if (!(IS_LOGICAL(compress) && LENGTH(compress) == 1 &&
          LOGICAL(compress)[0] != NA_LOGICAL))
        Rf_error("'%s' must be '%s'", "compress", "logical(1) (TRUE or FALSE)");
    const int compress1 = LOGICAL(compress)[0];
    if (!(IS_INTEGER(max_width) && LENGTH(max_width) == 1 &&
          INTEGER(max_width)[0] >= 0))
        Rf_error("'%s' must be %s", "max_width", "'integer(1)', >=0");
    const int width = INTEGER(max_width)[0];

    DECODE_FUNC dnaDecoder = decoder(get_XStringSet_xsbaseclassname(sread));
    XStringSet_holder xid = hold_XStringSet(id),
        xsread = hold_XStringSet(sread), xquality = hold_XStringSet(quality);

    char *idbuf0 = (char *) R_alloc(sizeof(char), width + 1), *idbuf1,
        *readbuf = (char *) R_alloc(sizeof(char), width + 1),
        *qualbuf = (char *) R_alloc(sizeof(char), width + 1),
        *gzbuf = NULL;
    int i, gzbuf_n;
    idbuf1 = TRUE == LOGICAL(full)[0] ? idbuf0 : "";

    FILE *fout = NULL;
    gzFile gzout = NULL;

    if (compress1 == FALSE)
        fout = fopen(CHAR(STRING_ELT(fname, 0)), CHAR(STRING_ELT(fmode, 0)));
    else {
        gzout = gzopen(CHAR(STRING_ELT(fname, 0)), CHAR(STRING_ELT(fmode, 0)));
        gzbuf_n = 4 * width + 8; /* liberal */
        gzbuf = (char *) R_alloc(sizeof(char), gzbuf_n);
    }
    if ((gzout == NULL) && (fout == NULL))
        Rf_error("failed to open file '%s'", CHAR(STRING_ELT(fname, 0)));
    const char *fmt = "@%s\n%s\n+%s\n%s\n";
    int err = 0;
    for (i = 0; i < len; ++i) {
        idbuf0 = _holder_to_char(&xid, i, idbuf0, width, NULL);
        if (idbuf0 == NULL) { err = 1; break; }
        readbuf = _holder_to_char(&xsread, i, readbuf, width, dnaDecoder);
        if (readbuf == NULL) { err = 1; break; }
        qualbuf = _holder_to_char(&xquality, i, qualbuf, width, NULL);
        if (qualbuf == NULL) { err = 1; break; }
        if (compress1) {
            int n_out =
                snprintf(gzbuf, gzbuf_n, fmt, idbuf0, readbuf, idbuf1, qualbuf);
            if (n_out > gzbuf_n) {
                /* happens rarely, e.g., identifiers longer than sequence */
                gzbuf_n = n_out + 1;
                gzbuf = (char *) R_alloc(sizeof(char), gzbuf_n);
                snprintf(gzbuf, gzbuf_n, fmt, idbuf0, readbuf, idbuf1, qualbuf);
            }
            if (gzputs(gzout, gzbuf) == -1) { err = 1; break; }
        } else
            if (fprintf(fout, fmt, idbuf0, readbuf, idbuf1, qualbuf) < 0) {
                err = 1; break;
            };
    }
    if (compress1)
        gzclose(gzout);
    else
        fclose(fout);
    if (err != 0)
        _write_err(i);

    return R_NilValue;
}

/*
 * solexa/IPAR .*_int.txt.p.gz file
 */

void _count_ipar_int_recs(gzFile * file, int *n_recs, int *n_cycles)
{
    const char CYCLE_END = '#';
    const int LINEBUF_SIZE = 200001;
    size_t bytes_read = 0;
    char *buf = Calloc(LINEBUF_SIZE + 1, char);
    *n_recs = *n_cycles = 0;
    char *p = 0;
    /* records and cycles */
    while (*n_cycles == 0 && (bytes_read = gzread(file, buf, LINEBUF_SIZE)) > 0) {
        p = buf;
        while ((p = memchr(p, '\n', (buf + bytes_read) - p))) {
            ++p;
            if (*p == CYCLE_END) {
                ++p;
                *n_cycles += 1;
                break;
            } else
                *n_recs += 1;
        }
    }
    /* just cycles */
    while ((p = memchr(p, CYCLE_END, (buf + bytes_read) - p))) {
        ++p;
        *n_cycles += 1;
    }
    while ((bytes_read = gzread(file, buf, LINEBUF_SIZE)) > 0) {
        p = buf;
        while ((p = memchr(p, CYCLE_END, (buf + bytes_read) - p))) {
            ++p;
            *n_cycles += 1;
        }
    }
    Free(buf);
}

SEXP count_ipar_int_recs(SEXP fnames)
{
    int i, nfile;
    const char *filepath;
    gzFile *file;
    SEXP ans = R_NilValue, nms = R_NilValue;

    if (!IS_CHARACTER(fnames))
        error("'fnames' must be character()");

    nfile = LENGTH(fnames);
    PROTECT(ans = NEW_LIST(2));
    SET_VECTOR_ELT(ans, 0, NEW_INTEGER(nfile));
    SET_VECTOR_ELT(ans, 1, NEW_INTEGER(nfile));
    PROTECT(nms = NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("reads"));
    SET_STRING_ELT(nms, 1, mkChar("cycles"));
    setAttrib(ans, R_NamesSymbol, nms);
    for (i = 0; i < nfile; ++i) {
        R_CheckUserInterrupt();
        filepath = translateChar(STRING_ELT(fnames, i));
        file = _fopen(filepath, "rb");
        _count_ipar_int_recs(file,
                             INTEGER(VECTOR_ELT(ans, 0)) + i,
                             INTEGER(VECTOR_ELT(ans, 1)) + i);
        gzclose(file);
    }

    UNPROTECT(2);
    return ans;
}

/*
 * Read a solexa .*_prb.txt file into STRING_VEC
 */

SEXP read_prb_as_character(SEXP fname, SEXP asSolexa)
{
    const int NUC_PER_CYCLE = 4;

    if (!IS_CHARACTER(fname) || LENGTH(fname) != 1)
        error("'fname' must be 'character(1)'");
    if (!IS_LOGICAL(asSolexa) || LENGTH(asSolexa) != 1)
        error("'asSolexa' must be 'logical(1)'");
    const int n_reads = INTEGER(count_lines(fname))[0];
    const int qbase = LOGICAL(asSolexa)[0] ? SOLEXA_QBASE : PHRED_QBASE;
    SEXP ans = PROTECT(NEW_CHARACTER(n_reads));

    gzFile *file = _fopen(translateChar(STRING_ELT(fname, 0)), "rb");
    char buf[LINEBUF_SIZE + 1];
    int read = 0;
    if (gzgets(file, buf, LINEBUF_SIZE) == Z_NULL) {
        gzclose(file);
        error("could not read file '%f'", translateChar(STRING_ELT(fname, 0)));
    }
    int n_cycles = 0;
    char *quad = strtok(buf, "\t");
    while (quad != NULL) {
        n_cycles++;
        quad = strtok(NULL, "\t");
    }
    gzrewind(file);

    char *score = R_alloc(sizeof(char), n_cycles + 1);
    score[n_cycles] = '\0';

    while (gzgets(file, buf, LINEBUF_SIZE) != Z_NULL) {
        if (read >= n_reads) {
            gzclose(file);
            error("too many reads, %d expected", n_reads);
        }
        quad = strtok(buf, "\t");
        int cycle = 0;
        while (quad != NULL && cycle < n_cycles) {
            int v[4];
            int bases = sscanf(quad, " %d %d %d %d",
                               &v[0], &v[1], &v[2], &v[3]);
            if (bases != NUC_PER_CYCLE) {
                gzclose(file);
                error("%d bases observed, %d expected", bases, NUC_PER_CYCLE);
            }
            v[0] = v[0] > v[1] ? v[0] : v[1];
            v[2] = v[2] > v[3] ? v[2] : v[3];
            score[cycle++] = qbase + ((char) v[0] > v[2] ? v[0] : v[2]);
            quad = strtok(NULL, "\t");
        }
        if (cycle != n_cycles) {
            gzclose(file);
            error("%d cycles observed, %d expected", cycle, n_cycles);
        }
        SET_STRING_ELT(ans, read++, mkChar(score));
    }
    UNPROTECT(1);
    gzclose(file);
    return ans;
}

/*
 * Read a solexa 's_<lane>_sequence.txt' file into CharAEAE objects.
 */
static void _read_solexa_fastq_file(const char *fname, SEXP ans)
{
    gzFile *file;
    char linebuf[LINEBUF_SIZE];
    int lineno, reclineno, nchar_in_buf;
    _XSnap seq = VECTOR_ELT(ans, 0),
        id = VECTOR_ELT(ans, 1), qualities = VECTOR_ELT(ans, 2);

    file = _fopen(fname, "rb");
    lineno = 0;
    while (gzgets(file, linebuf, LINEBUF_SIZE) != NULL) {
        if ((reclineno = lineno % LINES_PER_FASTQ_REC) == 2) {
            lineno++;
            continue;
        }

        nchar_in_buf = _rtrim(linebuf);
        if (nchar_in_buf >= LINEBUF_SIZE - 1) {	// should never be
            gzclose(file);
            error("line too long %s:%d", fname, lineno);
        } else if ((0 == reclineno) && (0 == nchar_in_buf)) {
            gzclose(file);
            error("unexpected empty line %s:%d", fname, lineno);
        }
        switch (reclineno) {
        case 0:
            /* add linebuf to CharAEAE; start at char +1 to skip the
             * fastq annotation. */
            if (id != R_NilValue)
                _APPEND_XSNAP(id, linebuf + 1);
            break;
        case 1:
            _solexa_to_IUPAC(linebuf);
            _APPEND_XSNAP(seq, linebuf);
            break;
        case 3:
            _APPEND_XSNAP(qualities, linebuf);
            break;
        default:
            error("unexpected 'reclineno'; consult maintainer");
            break;
        }
        lineno++;
    }
    gzclose(file);
    if ((lineno % LINES_PER_FASTQ_REC) != 0)
        error("unexpected number of lines in file '%s'", fname);
}

SEXP read_solexa_fastq(SEXP files, SEXP withId)
{
    int i, nfiles, nrec = 0;
    const char *fname;
    SEXP ans = R_NilValue, nms = R_NilValue;

    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character");
    if (!IS_LOGICAL(withId) || LENGTH(withId) != 1)
        Rf_error("'%s' must be '%s'", "withId", "logical(1)");

    nfiles = LENGTH(files);
    nrec = _count_lines_sum(files) / LINES_PER_FASTQ_REC;
    PROTECT(ans = NEW_LIST(3));
    SET_VECTOR_ELT(ans, 0, _NEW_XSNAP(nrec, "DNAString"));	/* sread */
    if (LOGICAL(withId)[0] == TRUE)	/* id */
        SET_VECTOR_ELT(ans, 1, _NEW_XSNAP(nrec, "BString"));
    else
        SET_VECTOR_ELT(ans, 1, R_NilValue);
    SET_VECTOR_ELT(ans, 2, _NEW_XSNAP(nrec, "BString"));	/* quality */

    PROTECT(nms = NEW_CHARACTER(3));
    SET_STRING_ELT(nms, 0, mkChar("sread"));
    SET_STRING_ELT(nms, 1, mkChar("id"));
    SET_STRING_ELT(nms, 2, mkChar("quality"));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(1);

    for (i = 0; i < nfiles; ++i) {
        R_CheckUserInterrupt();
        fname = translateChar(STRING_ELT(files, i));
        _read_solexa_fastq_file(fname, ans);
    }
    _XSNAP_ELT(ans, 0);
    if (VECTOR_ELT(ans, 1) != R_NilValue)
        _XSNAP_ELT(ans, 1);
    _XSNAP_ELT(ans, 2);

    UNPROTECT(1);
    return ans;
}

int _io_XStringSet_columns(const char *fname, int header,
                           const char *sep, MARK_FIELD_FUNC * mark_field,
                           const int *colidx, int ncol,
                           int nrow, int skip, const char *commentChar,
                           SEXP sets, const int *toIUPAC)
{
    gzFile *file;
    char *linebuf;
    int lineno = 0, recno = 0;

    file = _fopen(fname, "rb");
    linebuf = S_alloc(LINEBUF_SIZE, sizeof(char));	/* auto free'd */

    while (skip-- > 0)
        gzgets(file, linebuf, LINEBUF_SIZE);
    if (header == TRUE)
        gzgets(file, linebuf, LINEBUF_SIZE);

    while (recno < nrow && gzgets(file, linebuf, LINEBUF_SIZE) != NULL) {
        if (_linebuf_skip_p(linebuf, file, fname, lineno, commentChar)) {
            lineno++;
            continue;
        }

        int j = 0, cidx = 0;
        char *curr = linebuf, *next;
        for (j = 0; cidx < ncol && curr != NULL; ++j) {
            next = (*mark_field) (curr, sep);
            if (j == colidx[cidx]) {
                if (toIUPAC[cidx])
                    _solexa_to_IUPAC(curr);
                _APPEND_XSNAP(VECTOR_ELT(sets, cidx), curr);
                cidx++;
            }
            curr = next;
        }
        lineno++;
        recno++;
    }
    gzclose(file);
    return recno;
}

SEXP read_XStringSet_columns(SEXP files, SEXP header, SEXP sep,
                             SEXP colIndex, SEXP colClasses,
                             SEXP nrows, SEXP skip, SEXP commentChar)
{
    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character(1)");
    if (!IS_LOGICAL(header) || LENGTH(header) != 1)
        Rf_error("'%s' must be '%s'", "header", "logical(1)");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1)
        Rf_error("'%s' must be '%s'", "sep", "character(1)");
    /* FIXME: !nzchar(sep[1]) */
    if (!IS_INTEGER(colIndex) || LENGTH(colIndex) == 0)
        Rf_error("'colIndex' must be 'integer' with length > 0");
    if (!IS_CHARACTER(colClasses) || LENGTH(colClasses) != LENGTH(colIndex))
        Rf_error("'%s' must be '%s', length(colClasses) == length(colIndex)",
                 "colClasses", "character()");
    if (!IS_INTEGER(nrows) || LENGTH(nrows) != 1)
        Rf_error("'%s' must be '%s'", "nrows", "integer(1)");
    if (!IS_INTEGER(skip) || LENGTH(skip) != 1)
        Rf_error("'%s' must be '%s'", "skiip", "integer(1)");
    if (!IS_CHARACTER(commentChar) || LENGTH(commentChar) != 1)
        Rf_error("'%s' must be '%s'", "commentChar", "character(1)");
    if (LENGTH(STRING_ELT(commentChar, 0)) != 1)
        Rf_error("'nchar(commentChar[[1]])' must be 1 but is %d",
                 LENGTH(STRING_ELT(commentChar, 0)));

    int i, j;
    /* Count lines and pre-allocate space */
    const char *csep = translateChar(STRING_ELT(sep, 0));
    const int nfiles = LENGTH(files);
    MARK_FIELD_FUNC *sep_func;  /* how to parse fields; minor efficiency */
    if (csep[0] != '\0' && csep[1] == '\0')
        sep_func = _mark_field_1;
    else
        sep_func = _mark_field_n;

    int nrow = INTEGER(nrows)[0];
    if (nrow < 0) {
        nrow = _count_lines_sum(files);
        nrow -= nfiles * (LOGICAL(header)[0] + INTEGER(skip)[0]);
    }

    int ncol = LENGTH(colIndex);
    SEXP ans = PROTECT(NEW_LIST(ncol));
    int *colidx = (int *) R_alloc(sizeof(int), ncol);
    int *toIUPAC = (int *) R_alloc(sizeof(int), ncol);
    for (j = 0; j < ncol; ++j) {
        const char *baseclass = CHAR(STRING_ELT(colClasses, j));
        SET_VECTOR_ELT(ans, j, _NEW_XSNAP(nrow, baseclass));
        colidx[j] = INTEGER(colIndex)[j] - 1;
        toIUPAC[j] = !strcmp(baseclass, "DNAString");
    }

    /* read columns */
    int nreads = 0;
    for (i = 0; i < nfiles; ++i) {
        R_CheckUserInterrupt();
        if (nreads >= nrow)
            break;
        const char *fname = translateChar(STRING_ELT(files, i));
        nreads +=
            _io_XStringSet_columns(fname, LOGICAL(header)[0], csep, sep_func,
                                   colidx, ncol, nrow - nreads,
                                   INTEGER(skip)[0],
                                   CHAR(STRING_ELT(commentChar, 0)), ans,
                                   toIUPAC);
    }

    /* formulate return value */
    for (j = 0; j < ncol; ++j)
        _XSNAP_ELT(ans, j);
    UNPROTECT(1);
    return ans;
}

/*
 * _export parser
 */

enum {
    /* fields from the _export spec */
    SLX_MACHINE = 0, SLX_RUN, SLX_LANE, SLX_TILE, SLX_X, SLX_Y,
    SLX_MULTIPLEX, SLX_PAIRID, SLX_SREAD, SLX_QUAL, SLX_CHR,
    SLX_CONTIG, SLX_POS, SLX_STRAND, SLX_ALIGNQUAL, SLX_FILT,
    /* ID, when calculated */
    SLX_ID,
    /* length of ENUM */
    SLX_ELEMENT_END
};

SEXP _AlignedRead_Solexa_make(SEXP fields)
{
    const char *FILTER_LEVELS[] = { "Y", "N" };
    SEXP s, t, nmspc = PROTECT(_get_namespace("ShortRead"));
    const Rboolean
        withMultiplexIndex = R_NilValue != VECTOR_ELT(fields, SLX_MULTIPLEX),
        withPairedReadNumber = R_NilValue != VECTOR_ELT(fields, SLX_PAIRID),
        withIds = R_NilValue != VECTOR_ELT(fields, SLX_MACHINE);

    SEXP sfq;                   /* SFastqQuality */
    NEW_CALL(s, t, "SFastqQuality", nmspc, 2);
    CSET_CDR(t, "quality", VECTOR_ELT(fields, SLX_QUAL));
    CEVAL_TO(s, nmspc, sfq);
    PROTECT(sfq);

    SEXP alnq;                  /* NumericQuality() */
    NEW_CALL(s, t, "NumericQuality", nmspc, 2);
    CSET_CDR(t, "quality", VECTOR_ELT(fields, SLX_ALIGNQUAL));
    CEVAL_TO(s, nmspc, alnq);
    PROTECT(alnq);

    /* .SolexaExport_AlignedDataFrame(...) */
    _as_factor(VECTOR_ELT(fields, SLX_FILT), FILTER_LEVELS,
               sizeof(FILTER_LEVELS) / sizeof(const char *));
    SEXP run;
    NEW_CALL(s, t, "factor", nmspc, 2);
    CSET_CDR(t, "x", VECTOR_ELT(fields, SLX_RUN));
    CEVAL_TO(s, nmspc, run);
    PROTECT(run);

    SEXP dataframe;
    NEW_CALL(s, t, "data.frame", nmspc,
             8 + withMultiplexIndex + withPairedReadNumber);
    CSET_CDR(t, "run", run);
    CSET_CDR(t, "lane", VECTOR_ELT(fields, SLX_LANE));
    CSET_CDR(t, "tile", VECTOR_ELT(fields, SLX_TILE));
    CSET_CDR(t, "x", VECTOR_ELT(fields, SLX_X));
    CSET_CDR(t, "y", VECTOR_ELT(fields, SLX_Y));
    CSET_CDR(t, "filtering", VECTOR_ELT(fields, SLX_FILT));
    CSET_CDR(t, "contig", VECTOR_ELT(fields, SLX_CONTIG));
    if (withMultiplexIndex) {
        CSET_CDR(t, "multiplexIndex", VECTOR_ELT(fields, SLX_MULTIPLEX));
    }
    if (withPairedReadNumber) {
        CSET_CDR(t, "pairedReadNumber", VECTOR_ELT(fields, SLX_PAIRID));
    }
    CEVAL_TO(s, nmspc, dataframe);
    PROTECT(dataframe);

    SEXP adf;
    NEW_CALL(s, t, ".SolexaExport_AlignedDataFrame", nmspc, 2);
    CSET_CDR(t, "data", dataframe);
    CEVAL_TO(s, nmspc, adf);
    PROTECT(adf);

    SEXP aln;
    NEW_CALL(s, t, "AlignedRead", nmspc, 8 + withIds);
    CSET_CDR(t, "sread", VECTOR_ELT(fields, SLX_SREAD));
    CSET_CDR(t, "quality", sfq);
    if (withIds) {
        CSET_CDR(t, "id", VECTOR_ELT(fields, SLX_ID));
    }
    CSET_CDR(t, "chromosome", VECTOR_ELT(fields, SLX_CHR));
    CSET_CDR(t, "position", VECTOR_ELT(fields, SLX_POS));
    CSET_CDR(t, "strand", VECTOR_ELT(fields, SLX_STRAND));
    CSET_CDR(t, "alignQuality", alnq);
    CSET_CDR(t, "alignData", adf);
    CEVAL_TO(s, nmspc, aln);

    UNPROTECT(6);
    return aln;
}

int _read_solexa_export_file(const char *fname, const char *commentChar,
                             int offset, SEXP result)
{
    const int N_FIELDS = 22;
    Rboolean
        withMultiplexIndex = R_NilValue != VECTOR_ELT(result,
                                                      SLX_MULTIPLEX),
        withPairedReadNumber = R_NilValue != VECTOR_ELT(result,
                                                        SLX_PAIRID),
        withId = R_NilValue != VECTOR_ELT(result, SLX_MACHINE);

    gzFile *file;
    char linebuf[LINEBUF_SIZE],
	**elt = (char **) R_alloc(N_FIELDS, sizeof(char*));
    int lineno = 0, irec = offset;

    SEXP machine = NULL, run = VECTOR_ELT(result, SLX_RUN);
    int *lane = INTEGER(VECTOR_ELT(result, SLX_LANE)),
        *tile = INTEGER(VECTOR_ELT(result, SLX_TILE)),
        *x = INTEGER(VECTOR_ELT(result, SLX_X)),
        *y = INTEGER(VECTOR_ELT(result, SLX_Y));
    _XSnap sread = VECTOR_ELT(result, SLX_SREAD),
        quality = VECTOR_ELT(result, SLX_QUAL);
    SEXP chromosome = VECTOR_ELT(result, SLX_CHR);
    int *position = INTEGER(VECTOR_ELT(result, SLX_POS)),
        *strand = INTEGER(VECTOR_ELT(result, SLX_STRAND)),
        *alignQuality = INTEGER(VECTOR_ELT(result, SLX_ALIGNQUAL)),
        *filtering = INTEGER(VECTOR_ELT(result, SLX_FILT));
    SEXP contig = VECTOR_ELT(result, SLX_CONTIG), multiplexIndex = NULL;
    int *pairedReadNumber = NULL;
    if (withMultiplexIndex)
        multiplexIndex = VECTOR_ELT(result, SLX_MULTIPLEX);
    if (withPairedReadNumber)
        pairedReadNumber = INTEGER(VECTOR_ELT(result, SLX_PAIRID));
    if (withId)
        machine = VECTOR_ELT(result, SLX_MACHINE);

    file = _fopen(fname, "rb");
    while (gzgets(file, linebuf, LINEBUF_SIZE) != NULL) {
        if (*linebuf == *commentChar) {
            lineno++;
            continue;
        }

        /* field-ify */
        int n_fields = _mark_field_0(linebuf, elt, N_FIELDS);
        if (n_fields != N_FIELDS) {
            gzclose(file);
            error("incorrect number of fields (%d) %s:%d",
                  n_fields, fname, lineno);
        }

        if (withId)
            SET_STRING_ELT(machine, irec, mkChar(elt[0]));
        SET_STRING_ELT(run, irec, mkChar(elt[1]));
        lane[irec] = atoi(elt[2]);
        tile[irec] = atoi(elt[3]);
        x[irec] = atoi(elt[4]);
        y[irec] = atoi(elt[5]);
        if (withMultiplexIndex) {
            SEXP idxString = *elt[6] == '\0' ? R_BlankString : mkChar(elt[6]);
            SET_STRING_ELT(multiplexIndex, irec, idxString);
        }
        if (withPairedReadNumber)
            pairedReadNumber[irec] = atoi(elt[7]);
        _APPEND_XSNAP(sread, elt[8]);
        _APPEND_XSNAP(quality, elt[9]);
        SET_STRING_ELT(chromosome, irec, mkChar(elt[10]));
        SET_STRING_ELT(contig, irec, mkChar(elt[11]));
        if (*elt[12] == '\0')
            position[irec] = NA_INTEGER;
        else
            position[irec] = atoi(elt[12]);
        if (*elt[13] == '\0')
            strand[irec] = NA_INTEGER;
        else {
            switch (*elt[13]) {
            case 'F':
                strand[irec] = 1;
                break;
            case 'R':
                strand[irec] = 2;
                break;
            default:
                gzclose(file);
                error("invalid 'strand' field '%s', %s:%d",
                      *elt[13], fname, lineno);
                break;
            }
        }
        /* 14: descriptor */
        alignQuality[irec] = atoi(elt[15]);
        /* 16: pairedScore, 17: partnerCzome, 18: partnerContig
           19: partnerOffset, 20: partnerStrand */
        switch (*elt[21]) {
        case 'Y':
            filtering[irec] = 1;
            break;
        case 'N':
            filtering[irec] = 2;
            break;
        default:
            gzclose(file);
            error("invalid 'filtering' field '%s', %s:%d",
                  *elt[21], fname, lineno);
            break;
        }
        lineno++;
        irec++;
    }

    return irec - offset;
}

int _solexa_export_make_id(SEXP result)
{
    const int
    *lane = INTEGER(VECTOR_ELT(result, SLX_LANE)),
        *tile = INTEGER(VECTOR_ELT(result, SLX_TILE)),
        *x = INTEGER(VECTOR_ELT(result, SLX_X)),
        *y = INTEGER(VECTOR_ELT(result, SLX_Y)), *pairedReadNumber = NULL;
    const SEXP
        * run = STRING_PTR(VECTOR_ELT(result, SLX_RUN)),
        *multiplexIndex = NULL,
        *machine = STRING_PTR(VECTOR_ELT(result, SLX_MACHINE));
    const Rboolean
        withMultiplexIndex = R_NilValue != VECTOR_ELT(result, SLX_MULTIPLEX),
        withPairedReadNumber = R_NilValue != VECTOR_ELT(result, SLX_PAIRID);
    if (withMultiplexIndex)
        multiplexIndex = STRING_PTR(VECTOR_ELT(result, SLX_MULTIPLEX));
    if (withPairedReadNumber)
        pairedReadNumber = INTEGER(VECTOR_ELT(result, SLX_PAIRID));

    const int nrec = LENGTH(VECTOR_ELT(result, SLX_RUN));
    char buf[LINEBUF_SIZE];
    SET_VECTOR_ELT(result, SLX_ID, _NEW_XSNAP(nrec, "BString"));

    _XSnap id = VECTOR_ELT(result, SLX_ID);
    /* FIXME: machine */
    int n = 0;
    for (int i = 0; i < nrec; ++i) {
        n = snprintf(buf, LINEBUF_SIZE,
                     "%s_%s:%d:%d:%d:%d", CHAR(machine[i]),
                     CHAR(run[i]), lane[i], tile[i], x[i], y[i]);
        if (withMultiplexIndex)
            n += snprintf(buf + n, LINEBUF_SIZE - n, "#%s",
                          CHAR(multiplexIndex[i]));
        if (withPairedReadNumber)
            n += snprintf(buf + n, LINEBUF_SIZE - n, "/%d",
                          pairedReadNumber[i]);
        if (n > LINEBUF_SIZE)
            return -1;
        _APPEND_XSNAP(id, buf);
    }
    _XSNAP_ELT(result, SLX_ID);
    return 1;
}

SEXP read_solexa_export(SEXP files, SEXP sep, SEXP commentChar, SEXP withFlags)
{
    const int N_ELTS = SLX_ELEMENT_END;

    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character()");
    if (!IS_CHARACTER(sep) || LENGTH(sep) != 1 ||
        *(CHAR(STRING_ELT(sep, 0))) != '\t')
        Rf_error("'%s' must be '%s'", "sep", "\t");
    /* FIXME: !nzchar(sep[1]) */
    if (!IS_CHARACTER(commentChar) || LENGTH(commentChar) != 1)
        Rf_error("'%s' must be '%s'", "commentChar", "character(1)");
    if (LENGTH(STRING_ELT(commentChar, 0)) != 1)
        Rf_error("'nchar(commentChar[[1]])' must be 1 but is %d",
                 LENGTH(STRING_ELT(commentChar, 0)));
    if (!IS_LOGICAL(withFlags) || LENGTH(withFlags) != 3)
        Rf_error("'%s' must be '%s'", "withFlags", "logical(3)");
    Rboolean
        withId = LOGICAL(withFlags)[0],
        withMultiplexIndex = LOGICAL(withFlags)[1],
        withPairedReadNumber = LOGICAL(withFlags)[2];

    int nrec = _count_lines_sum(files);

    SEXP result = PROTECT(NEW_LIST(N_ELTS));;
    if (withId)
        SET_VECTOR_ELT(result, SLX_MACHINE, NEW_STRING(nrec));
    SET_VECTOR_ELT(result, SLX_RUN, NEW_STRING(nrec));
    SET_VECTOR_ELT(result, SLX_LANE, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_TILE, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_X, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_Y, NEW_INTEGER(nrec));
    if (withMultiplexIndex)
        SET_VECTOR_ELT(result, SLX_MULTIPLEX, NEW_STRING(nrec));
    if (withPairedReadNumber)
        SET_VECTOR_ELT(result, SLX_PAIRID, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_SREAD, _NEW_XSNAP(nrec, "DNAString"));
    SET_VECTOR_ELT(result, SLX_QUAL, _NEW_XSNAP(nrec, "BString"));
    SET_VECTOR_ELT(result, SLX_CHR, NEW_STRING(nrec));
    SET_VECTOR_ELT(result, SLX_POS, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_STRAND, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_ALIGNQUAL, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_FILT, NEW_INTEGER(nrec));
    SET_VECTOR_ELT(result, SLX_CONTIG, NEW_STRING(nrec));

    nrec = 0;
    for (int i = 0; i < LENGTH(files); ++i) {
        R_CheckUserInterrupt();
        nrec += _read_solexa_export_file(CHAR(STRING_ELT(files, i)),
                                         CHAR(STRING_ELT(commentChar, 0)),
                                         nrec, result);
    }

    _XSNAP_ELT(result, SLX_SREAD);
    _XSNAP_ELT(result, SLX_QUAL);

    SEXP strand_lvls = PROTECT(_get_strand_levels());
    _as_factor_SEXP(VECTOR_ELT(result, SLX_STRAND), strand_lvls);

    if (withId) {
        int ok = _solexa_export_make_id(result);
        if (ok <= 0) {
            UNPROTECT(2);
            Rf_error("internal error: could not make id");
        }
    }

    SEXP aln = _AlignedRead_Solexa_make(result);
    UNPROTECT(2);
    return aln;
}
