/* Note: This file has been copied at verbatim (except for the
present note) from the source code of Maq, version 0.6.5, and is
hence copyright (c) Li Hang, who has released Maq under GPL 2 */

#ifndef NST_CONST_H
#define NST_CONST_H

#define MAX_ULL 0xffffffffffffffffull

typedef unsigned long long bit64_t;
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
