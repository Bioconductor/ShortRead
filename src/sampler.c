#include "ShortRead.h"

SEXP _as_xstringset(SEXP width, SEXP tag, const char *classname,
                    const char *element_type)
{
    int nrec = LENGTH(width);
    SEXP start = PROTECT(NEW_INTEGER(nrec));
    int i, *s = INTEGER(start);
    const int *w = INTEGER(width);
    s[0] = 1;
    for (i = 1; i < nrec; ++i)
        s[i] = s[i - 1] + w[i - 1];
    SEXP rng = PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
    SEXP xstringset =
        PROTECT(new_XRawList_from_tag(classname, element_type, tag, rng));
    UNPROTECT(3);
    return xstringset;
}

Rbyte *_fastq_record_end(Rbyte * buf, const Rbyte * bufend)
{
    int nchr = 0;
    if (*buf++ != '@')
        Rf_error("record does not start with '@'");
    while (buf != bufend && *buf++ != '\n') ;	/* id 1 */
    while (buf != bufend && *buf != '+')	/* read */
        if (*buf++ != '\n')
            ++nchr;
    while (buf != bufend && *buf++ != '\n') ;	/* id 2 */
    while (buf != bufend && nchr)	/* qual */
        if (*buf++ != '\n')
            --nchr;
    if (0 != nchr)
        buf = NULL;
    if (buf && buf != bufend && *buf++ != '\n')
        Rf_error("internal: buf != <newline>");
    return buf;
}

SEXP sampler_rec_counter(SEXP buffer)
{
    /* full and partial records */
    if (0 == LENGTH(buffer))
        return ScalarInteger(0);

    Rbyte *buf = RAW(buffer);
    const Rbyte *bufend = buf + LENGTH(buffer);
    int n = 0;

    while (buf && buf < bufend) {
        while (buf != bufend && *buf == '\n')
            ++buf;
        if (buf != bufend) {
            buf = _fastq_record_end(buf, bufend);
            ++n;
        }
    }
    return ScalarInteger(n);
}

SEXP sampler_rec_parser(SEXP buffer, SEXP rec_n)
{
    if (0 == INTEGER(rec_n)[0])
        return NEW_LIST(0);

    const int n = INTEGER(rec_n)[0];
    Rbyte *buf = RAW(buffer);
    const Rbyte *bufend = buf + LENGTH(buffer), *prev;
    SEXP lst, elt;
    int s_idx;

    lst = PROTECT(NEW_LIST(n));

    prev = buf = RAW(buffer);
    s_idx = 0;
    while (buf && buf < bufend) {
        while (buf < bufend && *buf == '\n')
            ++buf;
        prev = buf;
        if (NULL == (buf = _fastq_record_end(buf, bufend)))
            break;
        elt = NEW_RAW(buf - prev);
        SET_VECTOR_ELT(lst, s_idx++, elt);
        memcpy((Rbyte *) RAW(elt), prev, (buf - prev) * sizeof(Rbyte));
    }
    if (n > s_idx) {
        elt = NEW_RAW(bufend - prev);
        SET_VECTOR_ELT(lst, s_idx++, elt);
        memcpy((Rbyte *) RAW(elt), prev, (bufend - prev) * sizeof(Rbyte));
    }
    if (n != s_idx)
        Rf_error("internal: sample index (%d) != 'rec_n' (%d)",
                 s_idx, n);
    UNPROTECT(1);
    return lst;
}

SEXP sampler_as_fastq(SEXP records)
{
    SEXP ans = PROTECT(NEW_LIST(3));
    SEXP lengths = PROTECT(NEW_LIST(2));
    int nrec = LENGTH(records), totlen = 0, i;

    SET_VECTOR_ELT(lengths, 0, NEW_INTEGER(nrec));	/* sread / quality */
    SET_VECTOR_ELT(lengths, 1, NEW_INTEGER(nrec));	/* id */

    for (i = 0; i < nrec; ++i)
        totlen += LENGTH(VECTOR_ELT(records, i));
    SET_VECTOR_ELT(ans, 0, NEW_RAW(totlen / 2));	/* sread */
    SET_VECTOR_ELT(ans, 1, NEW_RAW(totlen / 2));	/* quality */
    SET_VECTOR_ELT(ans, 2, NEW_RAW(totlen / 2));	/* id */

    Rbyte *sread_offset = RAW(VECTOR_ELT(ans, 0)),
        *qual_offset = RAW(VECTOR_ELT(ans, 1)),
        *id_offset = RAW(VECTOR_ELT(ans, 2));
    int *sread_len = INTEGER(VECTOR_ELT(lengths, 0)),
        *id_len = INTEGER(VECTOR_ELT(lengths, 1));

    for (i = 0; i < nrec; ++i) {
        SEXP record = VECTOR_ELT(records, i);
        Rbyte *buf = RAW(record);
        Rbyte *bufend = buf + LENGTH(record);
        Rbyte *start, *curr;
        int nchr;

        /* id */
        start = ++buf;          /* skip '@' */
        while (*buf != '\n')
            ++buf;
        memcpy(id_offset, start, (buf - start) * sizeof(Rbyte));
        *id_len++ = buf - start;
        id_offset += buf - start;

        /* read */
        while (*buf == '\n')
            ++buf;
        start = curr = buf;
        while (*buf != '+') {
            /* strip '\n' */
            while (*buf != '\n')
                *curr++ = *buf++;
            buf++;
        }
        memcpy(sread_offset, start, (curr - start) * sizeof(Rbyte));
        *sread_len++ = curr - start;
        sread_offset += curr - start;
        nchr = curr - start;

        /* second id tag */
        while (*buf != '\n')
            ++buf;
        /* quality */
        while (*buf == '\n')
            ++buf;              /* leading '\n' */
        start = curr = buf;
        while (buf != bufend && curr - start != nchr) {
            if (*buf != '\n')
                *curr++ = *buf++;
            else
                buf++;
        }
        memcpy(qual_offset, start, (curr - start) * sizeof(Rbyte));
        qual_offset += curr - start;
    }

    Rbyte *dna = RAW(VECTOR_ELT(ans, 0));
    while (dna < sread_offset) {
        Rbyte tmp = DNAencode(*dna);
        *dna++ = tmp;
    }

    SETLENGTH(VECTOR_ELT(ans, 0), sread_offset - RAW(VECTOR_ELT(ans, 0)));
    SETLENGTH(VECTOR_ELT(ans, 1), qual_offset - RAW(VECTOR_ELT(ans, 1)));
    SETLENGTH(VECTOR_ELT(ans, 2), id_offset - RAW(VECTOR_ELT(ans, 2)));

    SEXP xsset;
    xsset = _as_xstringset(VECTOR_ELT(lengths, 0), VECTOR_ELT(ans, 0),
                           "DNAStringSet", "DNAString");
    SET_VECTOR_ELT(ans, 0, xsset);;
    xsset = _as_xstringset(VECTOR_ELT(lengths, 0), VECTOR_ELT(ans, 1),
                           "BStringSet", "BString");
    SET_VECTOR_ELT(ans, 1, xsset);;
    xsset = _as_xstringset(VECTOR_ELT(lengths, 1), VECTOR_ELT(ans, 2),
                           "BStringSet", "BString");
    SET_VECTOR_ELT(ans, 2, xsset);

    SEXP nms = PROTECT(NEW_CHARACTER(3));
    SET_STRING_ELT(nms, 0, mkChar("sread"));
    SET_STRING_ELT(nms, 1, mkChar("quality"));
    SET_STRING_ELT(nms, 2, mkChar("id"));
    SET_NAMES(ans, nms);

    UNPROTECT(3);
    return ans;
}
