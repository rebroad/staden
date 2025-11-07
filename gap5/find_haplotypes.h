#ifndef _FIND_HAPLOTYPES_H_
#define _FIND_HAPLOTYPES_H_

#include <tg_gio.h>

Array find_haplotypes(GapIO *io, contig_list_t *contigs, int ncontigs,
		      int pairs, float het_score, float discrep_score,
		      int min_count);

#endif /* _FIND_HAPLOTYPES_H_ */
