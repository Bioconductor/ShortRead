#include "trim.h"
#include "IRanges_interface.h"
#include "Biostrings_interface.h"

SEXP
trim_tails(SEXP quality, SEXP k, SEXP a_map, SEXP successive)
{
    SEXP end;
    int map[256];

    const cachedXStringSet cache = cache_XStringSet(quality);
    const int len = get_XStringSet_length(quality);
    int i, j, *endp;

    end = PROTECT(NEW_INTEGER(len));
    endp = INTEGER(end);
    for (j = 0; j < Rf_length(a_map); ++j) {
        const char c = CHAR(STRING_ELT(GET_NAMES(a_map), j))[0];
        map[(int) c] = INTEGER(a_map)[j];
    }

    const int kmax = INTEGER(k)[0];
    if (!LOGICAL(successive)[0]) {
#ifdef SUPPORT_OPENMP
#pragma omp parallel for private(j)
#endif
        for (i = 0; i < len; ++i) {
            const cachedCharSeq seq =
                get_cachedXStringSet_elt(&cache, i);
            int n = 0;
            for (j = 0; j < seq.length; ++j) {
                n += map[(int) seq.seq[j]];
                if (kmax <= n)
                    break;
                }
            endp[i] = j;
        }
    } else {
        const int nbuf = INTEGER(k)[0];
        int *kbuf = (int *) R_alloc(sizeof(int), nbuf), ibuf;
#ifdef SUPPORT_OPENMP
#pragma omp parallel for private(j)
#endif
        for (i = 0; i < len; ++i) {
            const cachedCharSeq seq =
                get_cachedXStringSet_elt(&cache, i);
            int n = 0;
            for (ibuf = 0; ibuf < nbuf; ++ibuf)
                kbuf[ibuf] = 0;
            int m;
            for (j = 0; j < seq.length; ++j) {
                m = map[(int) seq.seq[j]];
                n += m - kbuf[j % nbuf];
                if (kmax <= n)
                    break;
                kbuf[j % nbuf] = m;
            }
            endp[i] = j == seq.length ? j : j - nbuf + 1L;
        }
    }

    UNPROTECT(1);
    return end;
}

