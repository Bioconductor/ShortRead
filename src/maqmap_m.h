/* Note: This file is based on the file maqmap.h of the source code 
of Maq, version 0.7.2, which is copyright (c) Li Hang, who has 
released Maq under GPL 2.
The changes to the original allow to switch the value of MAX_READLEN,
which is a preprocessor macro in heng's code, at run-time, because
Maq used 64 until 0.6.x, then (unless the macro MAQ_LONGREADS is not set)
the value 128. */

#ifndef MAQMAP_M_H_
#define MAQMAP_M_H_

#define MAX_READLEN_OLD 64
#define MAX_READLEN_NEW 128
#define MAX_NAMELEN 36

#define MAQMAP_FORMAT_OLD 0
#define MAQMAP_FORMAT_NEW -1

#define PAIRFLAG_FF      0x01
#define PAIRFLAG_FR      0x02
#define PAIRFLAG_RF      0x04
#define PAIRFLAG_RR      0x08
#define PAIRFLAG_PAIRED  0x10
#define PAIRFLAG_DIFFCHR 0x20
#define PAIRFLAG_NOMATCH 0x40
#define PAIRFLAG_SW      0x80

#include <string.h>
#include <zlib.h>
#include "const.h"

/*
  name: read name
  size: the length of the read
  seq: read sequence (see also below)
  seq[MAX_READLEN-1]: single end mapping quality (equals to map_qual if not paired)
  map_qual: the final mapping quality
  alt_qual: the lower quality of the two ends (equals to map_qual if not paired)
  flag: status of the pair
  dist: offset of the mate (zero if not paired)
  info1: mismatches in the 24bp (higher 4 bits) and mismatches (lower 4 bits)
  info2: sum of errors of the best hit
  c[2]: count of all 0- and 1-mismatch hits on the reference
 */
template< int max_readlen > struct maqmap1_T
{
	bit8_t seq[max_readlen]; /* the last base is the single-end mapping quality. */
	bit8_t size, map_qual, info1, info2, c[2], flag, alt_qual;
	bit32_t seqid, pos;
	int dist;
	char name[MAX_NAMELEN];
};

template< int max_readlen > struct maqmap_T
{
	int format, n_ref;
	char **ref_name;
	bit64_t n_mapped_reads;
	maqmap1_T< max_readlen > *mapped_reads;
};

template< int max_readlen > 
inline int maqmap_read1( gzFile fp, maqmap1_T<max_readlen> * m1 ) {
   return gzread( fp, m1, sizeof( maqmap1_T<max_readlen> ) );
}

template< int max_readlen > maqmap_T<max_readlen> *maq_new_maqmap()
{
	maqmap_T<max_readlen> *mm = 
	   (maqmap_T<max_readlen>*)calloc(1, sizeof(maqmap_T<max_readlen>));
	mm->format = MAQMAP_FORMAT_NEW;
	return mm;
}


template< int max_readlen > void maq_delete_maqmap(maqmap_T<max_readlen> *mm)
{
	int i;
	if (mm == 0) return;
	for (i = 0; i < mm->n_ref; ++i)
		free(mm->ref_name[i]);
	free(mm->ref_name);
	free(mm->mapped_reads);
	free(mm);
}

template< int max_readlen > maqmap_T<max_readlen> *maqmap_read_header(gzFile fp)
{
	maqmap_T<max_readlen> *mm;
	int k, len;
	mm = maq_new_maqmap<max_readlen>();
	gzread(fp, &mm->format, sizeof(int));
	if (mm->format != MAQMAP_FORMAT_NEW) {
		if (mm->format > 0) {
			maq_delete_maqmap(mm);
			error("obsolete map format; use MAQ 'mapass2maq' command to convert");
		}
		assert(mm->format == MAQMAP_FORMAT_NEW);
	}
	gzread(fp, &mm->n_ref, sizeof(int));
	mm->ref_name = (char**)calloc(mm->n_ref, sizeof(char*));
	for (k = 0; k != mm->n_ref; ++k) {
		gzread(fp, &len, sizeof(int));
		mm->ref_name[k] = (char*)malloc(len * sizeof(char));
		gzread(fp, mm->ref_name[k], len);
	}
	/* read number of mapped reads */
	gzread(fp, &mm->n_mapped_reads, sizeof(bit64_t));
	return mm;
}


#endif
