/* Note: This file has been copied at verbatim (except for the
present note) from the source code of Maq, version 0.6.5, and is
hence copyright (c) Li Hang, who has released Maq under GPL 2 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "const.h"
#include "maqmap.h"

maqmap_t *maq_new_maqmap()
{
	maqmap_t *mm = (maqmap_t*)calloc(1, sizeof(maqmap_t));
	mm->format = MAQMAP_FORMAT_NEW;
	return mm;
}
void maq_delete_maqmap(maqmap_t *mm)
{
	int i;
	if (mm == 0) return;
	for (i = 0; i < mm->n_ref; ++i)
		free(mm->ref_name[i]);
	free(mm->ref_name);
	free(mm->mapped_reads);
	free(mm);
}
void maqmap_write_header(gzFile fp, const maqmap_t *mm)
{
	int i, len;
	gzwrite(fp, &mm->format, sizeof(int));
	gzwrite(fp, &mm->n_ref, sizeof(int));
	for (i = 0; i != mm->n_ref; ++i) {
		len = strlen(mm->ref_name[i]) + 1;
		gzwrite(fp, &len, sizeof(int));
		gzwrite(fp, mm->ref_name[i], len);
	}
	gzwrite(fp, &mm->n_mapped_reads, sizeof(bit64_t));
}
maqmap_t *maqmap_read_header(gzFile fp)
{
	maqmap_t *mm;
	int k, len;
	mm = maq_new_maqmap();
	gzread(fp, &mm->format, sizeof(int));
	if (mm->format != MAQMAP_FORMAT_NEW) {
		if (mm->format > 0) {
			fprintf(stderr, "** Obsolete map format is detected. Please use 'mapass2maq' command to convert the format.\n");
			exit(3);
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

/* mapview */

static void mapview_core(FILE *fpout, gzFile fpin, int is_verbose, int is_mm)
{
	bit32_t j;
	maqmap_t *m = maqmap_read_header(fpin);
	maqmap1_t *m1, mm1;
	m1 = &mm1;
	while (maqmap_read1(fpin, m1)) {
		fprintf(fpout, "%s\t%s\t%d\t%c\t%d\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
				m1->name, m->ref_name[m1->seqid], (m1->pos>>1) + 1,
				(m1->pos&1)? '-' : '+', m1->dist, m1->flag, m1->map_qual, (signed char)m1->seq[MAX_READLEN-1],
				m1->alt_qual, m1->info1&0xf, m1->info2, m1->c[0], m1->c[1], m1->size);
		if (is_verbose) {
			fputc('\t', fpout);
			for (j = 0; j != m1->size; ++j) {
				if (m1->seq[j] == 0) fputc('n', fpout);
				else if ((m1->seq[j]&0x3f) < 27) fputc("acgt"[m1->seq[j]>>6&3], fpout);
				else fputc("ACGT"[m1->seq[j]>>6&3], fpout);
			}
			fputc('\t', fpout);
			for (j = 0; j != m1->size; ++j)
				fputc((m1->seq[j]&0x3f) + 33, fpout);
		}
		if (is_mm) {
			bit64_t *p = (bit64_t*)(m1->seq + 55);
			fprintf(fpout, "\t%llx", *p);
		}
		fputc('\n', fpout);
	}
	maq_delete_maqmap(m);
}

int ma_mapview(int argc, char *argv[])
{
	int c, is_verbose = 1, is_mm = 0;
	while ((c = getopt(argc, argv, "bN")) >= 0) {
		switch (c) {
		case 'b': is_verbose = 0; break;
		case 'N': is_mm = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: maq mapview [-bN] <in.map>\n");
		return 1;
	}
	gzFile fp = (strcmp(argv[optind], "-") == 0)? gzdopen(STDIN_FILENO, "r") : gzopen(argv[optind], "r");
	mapview_core(stdout, fp, is_verbose, is_mm);
	gzclose(fp);
	return 0;
}
