#include <R_ext/Utils.h>        /* R_CheckUserInterrupt */
#include <ctype.h>              /* isspace */
#include "ShortRead.h"

unsigned char _bDecode(char);
unsigned char _dnaDecode(char);
unsigned char _rnaDecode(char);

/*
 * Decode XStringSet wrappers
 */
unsigned char
_bDecode(char c)
{
    return (unsigned char) c;
}

unsigned char
_dnaDecode(char c)
{
    return (unsigned char) DNAdecode(c);
}

unsigned char
_rnaDecode(char c)
{
    return (unsigned char) RNAdecode(c);
}

DECODE_FUNC
decoder(const char* base)
{
    DECODE_FUNC decode;
    if (strcmp(base, "DNAString")==0) {
        decode = _dnaDecode;
    } else if (strcmp(base, "RNAString")==0) {
        decode = _rnaDecode;
    } else if (strcmp(base, "BString")==0) {
        decode = _bDecode;
    } else {
        Rf_error("unknown class '%s'", base);
    }
    return decode;
}

SEXP
_get_namespace(const char* pkg)
{
    SEXP fun = PROTECT(findFun(install("getNamespace"), R_GlobalEnv));
    SEXP nmspc = PROTECT(NEW_CHARACTER(1));
    SET_STRING_ELT(nmspc, 0, mkChar(pkg));
    nmspc = eval(lang2(fun, nmspc), R_GlobalEnv);
    UNPROTECT(2);
    return nmspc;
}

SEXP
_get_strand_levels()
{
    SEXP nmspc = PROTECT(_get_namespace("ShortRead"));
    SEXP ans = eval(findVar(install(".STRAND_LEVELS"), nmspc), nmspc);
    UNPROTECT(1);
    return ans;
}

int
_char_as_strand_int(const char c, const char *fname, const int lineno)
{
    int strand;
    if (c == '\0')
        strand = NA_INTEGER;
    else {
        switch (c) {
        case '-':
            strand = 1;
            break;
        case '+':
            strand = 2;
            break;
        default:
            error("invalid 'strand' field '%s', %s:%d",
                  c, fname, lineno);
            break;
        }
    }
    return strand;
}

/*
 * apply function 'with' to object 'from' in environment 'rho', e.g.,
 * becuase 'from' is an object and 'with' an accessor.
 */

SEXP
_get_SEXP(SEXP from, SEXP rho, const char *with)
{
    SEXP fun = PROTECT(findFun(install(with), rho));
    SEXP res = eval(lang2(fun, from), rho);
    UNPROTECT(1);
    return res;
}


/*
 * parse lines into fields.
 *
 * string is parsed until a character in delim is found, or end of
 * string reached.
 *
 * return value is pointer to the start of the next field, or NULL if
 * no more fields.
 */

char *
_mark_field_1(char *curr, const char *delim)
{
    char *c = curr;
    while (*c != '\0' && *c != *delim)
        c++;
    if (*c != '\0')             /* i.e., delim */
        *c++ = '\0';
    return c;
}

char *
_mark_field_n(char *curr, const char *delim)
{
    const char *d = '\0';
    while (*curr != '\0' && *curr != '\n') {
        d = delim;
        while (*d != '\0' && *d != *curr)
            ++d;
        if (*d != '\0')
            *curr = '\0';
        else
            ++curr;
    }
    if (*curr == '\n') {
        *curr = '\0';
        return '\0';
    }
    return *d == '\0' ? '\0' : curr + 1;
}

SEXP
_mark_field_test(SEXP filename, SEXP delimiters, SEXP dim)
{
    if (!IS_CHARACTER(filename)  || LENGTH(filename) !=1)
        error("'%s' must be '%s'", "filename", "character(1)");
    if (!IS_CHARACTER(delimiters) || LENGTH(delimiters) != 1)
        error("'%s' must be '%s'", "delimiters", "character(1)");
    if  (!IS_INTEGER(dim) || LENGTH(dim) != 2)
        error("'%s' must be '%s'", "dim", "integer(2)");

    SEXP ans = PROTECT(NEW_LIST(INTEGER(dim)[0]));
    int i, j;
    for (i = 0; i < INTEGER(dim)[0]; ++i)
        SET_VECTOR_ELT(ans, i, NEW_CHARACTER(INTEGER(dim)[1]));
    
#define LINEBUF_SIZE 1024
    FILE *file;
    char linebuf[LINEBUF_SIZE];
    if ((file = fopen(CHAR(STRING_ELT(filename, 0)), "rb")) == NULL)
        error("cannot open file '%s'", CHAR(STRING_ELT(filename, 0)));
    const char *delim = CHAR(STRING_ELT(delimiters, 0));

    for (i = 0; i < INTEGER(dim)[0]; ++i) {
        if (fgets(linebuf, LINEBUF_SIZE, file) == NULL)
            error("unexpected end-of-file");
        j = 0;
        char *curr = linebuf, *next;
        while (curr != NULL) {
            if (j >= INTEGER(dim)[1])
                error("too many fields");
            next = _mark_field_n(curr, delim);
            SET_STRING_ELT(VECTOR_ELT(ans, i), j, mkChar(curr));
            j++;
            curr = next;
        }
    }
#undef LINEBUF_SIZE

    UNPROTECT(1);
    return ans;
}

const int LINEBUF_SIZE = 20001;

/*
 * open and check file; signal error
 */
gzFile *
_fopen(const char *fname, const char *mode)
{
    gzFile *file = gzopen(fname, mode);
    if (file == NULL)
        error("cannot open file %s", fname);
    return file;
}

/*
 * trim & check linebuf, return 0 if processing should continue
 */
int
_linebuf_skip_p(char *linebuf, gzFile *file,
                const char *fname, int lineno, const char *commentChar)
{
    int nchar_in_buf;
    nchar_in_buf = _rtrim(linebuf);
    if (nchar_in_buf >= LINEBUF_SIZE - 1) { // should never be >
        gzclose(file);
        error("line too long %s:%d", fname, lineno);
    } else if (nchar_in_buf == 0) {
        gzclose(file);
        error("unexpected empty line %s:%d", fname, lineno);
    }
    return *linebuf == *commentChar;
}              

/*
 * Return the number of chars that remain in the buffer after we've removed
 * the right spaces ('\n', '\r', '\t', ' ', etc...)
 */
int
_rtrim(char *linebuf)
{
    int i;

    i = strlen(linebuf) - 1;
    while (i >= 0 && isspace(linebuf[i])) i--;
    linebuf[++i] = 0;
    return i;
}

/*
 * Solexa sometimes encodes an uncalled base as '.', but the
 * Biostrings standard is '-'. Convert a null-terminated character
 * string in-place.
 */
void
_solexa_to_IUPAC(char *p)
{
    while ((p = strchr(p, '.')) != NULL)
        *p = '-';
}

void
_reverse(char *linebuf)
{
    size_t len = strlen(linebuf);
    int i;
    char tmp;
    for (i = 0; i < floor(len / 2); ++i) {
        tmp = linebuf[len-i-1];
        linebuf[len-i-1] = linebuf[i];
        linebuf[i] = tmp;
    }
}

void
_reverseComplement(char *linebuf)
{
    static const int MAX_MAP = 256;
    static char map[256];
    static int init = 0;
    if (init == 0) {
        init = 1;
        for (int i = 0; i < MAX_MAP; ++i)
            map[i] = (char) i;
        map['A'] = 'T'; map['C'] = 'G'; map['G'] = 'C'; map['T'] = 'A';
        map['a'] = 't'; map['c'] = 'g'; map['g'] = 'c'; map['t'] = 'a';
    }
    _reverse(linebuf);
    for (int i = 0; i < strlen(linebuf); ++i)
        linebuf[i] = map[(int) linebuf[i]];
}

/*
 * Convert roSeq to named XStringSet
 */
SEXP
_CharAEAE_to_XStringSet(CharAEAE* aeae, const char *clsName)
{
    RoSeqs roSeqs = new_RoSeqs_from_CharAEAE(aeae);
    return new_XStringSet_from_RoSeqs(clsName, &roSeqs);
}

/*
 * Chenge vector class and attribute to represent factor
 */
void
_as_factor_SEXP(SEXP vec, SEXP lvls)
{
    SEXP cls = PROTECT(NEW_CHARACTER(1));
    SET_STRING_ELT(cls, 0, mkChar("factor"));
    SET_CLASS(vec, cls);
    SET_ATTR(vec, install("levels"), lvls);
    UNPROTECT(1);
}

void
_as_factor(SEXP vec, const char **levels, const int n_lvls)
{
    SEXP lvls = PROTECT(NEW_CHARACTER(n_lvls));
    int i;
    for (i = 0; i < n_lvls; ++i)
        SET_STRING_ELT(lvls, i, mkChar(levels[i]));
    _as_factor_SEXP(vec, lvls);
    UNPROTECT(1);
}

/*
 * Count the number of lines ('\n') in a file.
 *
 * file: an open file stream at position 0
 *
 */
static int
_count_lines(gzFile *file)
{
    const int LINEBUF_SIZE=20001;
    size_t bytes_read;
    char buf[LINEBUF_SIZE + 1];
    int lines = 0;

    while ((bytes_read = gzread(file, buf, LINEBUF_SIZE)) > 0) {
        char *p = buf;
        while ((p = memchr(p, '\n', (buf + bytes_read) - p))) {
            ++p;
            ++lines;
        }
    }
    return lines;
}

int
_count_lines_sum(SEXP files)
{
    SEXP nlines = count_lines(files);
    int i, nrec = 0;
    for (i = 0; i < LENGTH(files); ++i)
        nrec += INTEGER(nlines)[i];
    return nrec;
}

SEXP
count_lines(SEXP files)
{
    int i, nfile;
    const char *filepath;
    gzFile *file;
    SEXP ans = R_NilValue;
    
    if (!IS_CHARACTER(files))
        error("'files' must be character()");

    nfile = LENGTH(files);
    PROTECT(ans = NEW_INTEGER(nfile));
    for (i = 0; i < nfile; ++i) {
        R_CheckUserInterrupt();
        filepath = translateChar(STRING_ELT(files, i));
        file = _fopen(filepath, "rb");
        INTEGER(ans)[i] = _count_lines(file);
        gzclose(file);
    }

    UNPROTECT(1);
    return ans;
}
