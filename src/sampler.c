#include "ShortRead.h"

SEXP
_as_xstringset(SEXP width, SEXP tag, 
	       const char *classname, const char *element_type)
{
  int nrec = LENGTH(width);
  SEXP start = PROTECT(NEW_INTEGER(nrec));
  int i, *s = INTEGER(start);
  const int *w = INTEGER(width);
  s[0] = 1;
  for (i = 1; i < nrec; ++i) 
    s[i] = s[i-1] + w[i-1];
  SEXP rng =
    PROTECT(new_IRanges("IRanges", start, width, R_NilValue));
  SEXP xstringset =
    PROTECT(new_XRawList_from_tag(classname, element_type, tag, rng));
  UNPROTECT(3);
  return xstringset;
}

SEXP
sampler_rec_counter(SEXP buffer, SEXP recsep)
{
  if (0 == LENGTH(buffer))
    return ScalarInteger(0);

  Rbyte *buf = RAW(buffer);
  const Rbyte *bufend = buf + LENGTH(buffer);
  Rbyte *sep = RAW(recsep);
  const Rbyte *sepend = sep + LENGTH(recsep);
  int n = 0;

  /* count */
  /* buffer starts on a record boundary, but boundary might be
     'special', e.g., start of file; gooble any initial record
     start  */
  while (buf != bufend && sep != sepend && *buf == *sep) {
    ++buf;
    ++sep;
  }
  ++n;
  sep = RAW(recsep);

  while (buf != bufend) {
    if (*buf != *sep) ++buf;
    else {
      ++buf; ++sep;
      while (buf != bufend && sep != sepend && *buf == *sep) {
	++buf;
	++sep;
      }
      if (sep == sepend) ++n;
      sep = RAW(recsep);
    }
  }
  return ScalarInteger(n);
}

SEXP
sampler_rec_parser(SEXP buffer, SEXP recsep, SEXP sample)
{
  Rbyte *buf = RAW(buffer), *prev;
  const Rbyte *bufend = buf + LENGTH(buffer);
  Rbyte *sep = RAW(recsep);
  const int seplen = LENGTH(recsep);
  const Rbyte *sepend = sep + seplen;
  SEXP lst, elt;
  int i, s_idx;

  lst = NEW_LIST(LENGTH(sample));
  if (0 == LENGTH(lst))
    return lst;
  PROTECT(lst);

  prev = buf = RAW(buffer);
  sep = RAW(recsep);
  while (buf != bufend && sep != sepend && *buf == *sep) {
    ++buf;
    ++sep;
  }
  sep = RAW(recsep);
  i = 1;			/* 1-based R indicies */
  s_idx = 0;
  while (buf != bufend) {
    if (*buf != *sep) ++buf;
    else {
      ++buf; ++sep;
      while (buf != bufend && sep != sepend && *buf == *sep) {
	++buf; 
	++sep; 
      }
      if (sep == sepend) {
	if (i == INTEGER(sample)[s_idx]) {
	  elt = NEW_RAW(buf - seplen - prev);
	  SET_VECTOR_ELT(lst, s_idx++, elt);
	  memcpy((Rbyte *) RAW(elt), prev, 
		 (buf - prev - seplen) * sizeof(Rbyte));
	}
	++i;
	prev = buf - seplen;
      }
      sep = RAW(recsep);
    }
  }
  if (i == INTEGER(sample)[s_idx]) {
    elt = NEW_RAW(bufend - prev);
    SET_VECTOR_ELT(lst, s_idx, elt);
    memcpy((Rbyte *) RAW(elt), prev, (bufend - prev) * sizeof(Rbyte));
  }
  UNPROTECT(1);
  return lst;
}

SEXP
sampler_as_fastq(SEXP records)
{
  SEXP ans = PROTECT(NEW_LIST(3));
  SEXP lengths = PROTECT(NEW_LIST(2));
  int nrec = LENGTH(records), totlen = 0, i;

  SET_VECTOR_ELT(lengths, 0, NEW_INTEGER(nrec)); /* sread / quality */
  SET_VECTOR_ELT(lengths, 1, NEW_INTEGER(nrec)); /* id */

  for (i = 0; i < nrec; ++i) 
    totlen += LENGTH(VECTOR_ELT(records, i));
  SET_VECTOR_ELT(ans, 0, NEW_RAW(totlen / 2)); /* sread */
  SET_VECTOR_ELT(ans, 1, NEW_RAW(totlen / 2)); /* quality */
  SET_VECTOR_ELT(ans, 2, NEW_RAW(totlen / 4)); /* id */

  Rbyte 
    *sread_offset = RAW(VECTOR_ELT(ans, 0)),
    *qual_offset = RAW(VECTOR_ELT(ans, 1)),
    *id_offset = RAW(VECTOR_ELT(ans, 2));
  int 
    *sread_len = INTEGER(VECTOR_ELT(lengths, 0)),
    *id_len = INTEGER(VECTOR_ELT(lengths, 1));

  for (i = 0; i < nrec; ++i) {
    SEXP record = VECTOR_ELT(records, i);
    Rbyte *buf = RAW(record);
    Rbyte *bufend = buf + LENGTH(record);
    Rbyte *start, *curr;

    /* id */
    while (*buf == '\n') ++buf;
    start = ++buf;		/* skip '@' */
    while (*buf != '\n') ++buf;
    memcpy(id_offset, start, (buf - start) * sizeof(Rbyte));
    *id_len++ = buf - start;
    id_offset += buf - start;

    /* read */
    while (*buf == '\n') ++buf;
    start = curr = buf;
    while (*buf != '+') {
      /* strip '\n' */
      while (*buf != '\n') *curr++ = *buf++;
      buf++;
    }
    memcpy(sread_offset, start, (curr - start) * sizeof(Rbyte));
    *sread_len++ = curr - start;
    sread_offset += curr - start;

    /* quality */
    while (*buf != '\n') ++buf;	/* second id tag */
    while (*buf == '\n') ++buf;	/* leading '\n' */
    start = curr = buf;
    while (buf != bufend) {
      if (*buf != '\n') *curr++ = *buf++;
      else buf++;
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
