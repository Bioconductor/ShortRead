#ifndef TRIM_H
#define TRIM_H

#include <Rdefines.h>

SEXP trim_tails(SEXP quality, SEXP k, SEXP a_map, SEXP successive);
SEXP trim_tailw(SEXP quality, SEXP k, SEXP a_map, SEXP winsize);

#endif
