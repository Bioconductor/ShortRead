#include <R_ext/Utils.h>        /* R_CheckUserInterrupt */
#include <ctype.h>              /* isspace */
#include "ShortRead.h"


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
