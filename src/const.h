/* 

   Note: This file has been copied from the source code of Maq,
   version 0.6.5, and is hence copyright (c) Li Hang, who has released
   Maq under GPL 2.

   The test for SIZEOF_UNSIGNED_LONG == 8 has been added, and depends
   on ../configure.ac

*/

#ifndef NST_CONST_H
#define NST_CONST_H

#define MAX_ULL 0xffffffffffffffffull

#if SIZEOF_UNSIGNED_LONG == 8
typedef unsigned long bit64_t;
#else
typedef unsigned long long bit64_t;
#endif
typedef unsigned bit32_t;
typedef unsigned short bit16_t;
typedef unsigned char bit8_t;

extern bit8_t nst_nt4_table[];
extern bit8_t nst_nt16_table[];
extern char *nst_nt4_rev_table;
extern char *nst_nt16_rev_table;
extern bit8_t nst_nt16_nt4_table[];
extern int nst_nt16_count_table[];

#endif
