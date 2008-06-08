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
_mark_field(char *curr, const char *delim)
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
    if ((file = fopen(CHAR(STRING_ELT(filename, 0)), "r")) == NULL)
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
            next = _mark_field(curr, delim);
            SET_STRING_ELT(VECTOR_ELT(ans, i), j, mkChar(curr));
            j++;
            curr = next;
        }
    }
#undef LINEBUF_SIZE

    UNPROTECT(1);
    return ans;
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
_solexa_to_IUPAC(char *linebuf)
{
    char *p = linebuf;
    while ((p = strchr(p, '.')) != NULL)
        *p = '-';
}

/*
 * Count the number of lines ('\n') in a file.
 *
 * file: an open file stream at position 0
 *
 */
static int
_count_lines(FILE *file)
{
    const int LINEBUF_SIZE=20001;
    size_t bytes_read;
    char buf[LINEBUF_SIZE + 1];
    int lines = 0;

    while ((bytes_read = fread(buf, sizeof(char),
                               LINEBUF_SIZE, file)) != 0) {
        char *p = buf;
        while ((p = memchr(p, '\n', (buf + bytes_read) - p))) {
            ++p;
            ++lines;
        }
        bytes_read += bytes_read;
    }
    return lines;
}

SEXP
count_lines(SEXP files)
{
    int i, nfile;
    const char *filepath;
    FILE *file;
    SEXP ans = R_NilValue;
    
    if (!IS_CHARACTER(files))
        error("'files' must be character()");

    nfile = LENGTH(files);
    PROTECT(ans = NEW_INTEGER(nfile));
    for (i = 0; i < nfile; ++i) {
        R_CheckUserInterrupt();
        filepath = translateChar(STRING_ELT(files, i));
        if ((file = fopen(filepath, "r")) == NULL)
            error("cannot open file '%s'", filepath);
        INTEGER(ans)[i] = _count_lines(file);
        fclose(file);
    }

    UNPROTECT(1);
    return ans;
}
