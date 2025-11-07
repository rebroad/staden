/*
 * ---------------------------------------------------------------------------
 * Identifies SNP sites and forms haplotypes from these by linking SNPs
 * together based on spanning readings. 
 *
 * This isn't as advanced as Gap4's solution, having no knowledge of
 * read-pairs and no reverse correlation (only 2 haplotypes, so lack of
 * presence on one implies the other), but Gap4's version is very slow.
 * Here we take a quick and dirty approach.
 *
 * - Find primary SNPs and add to haplotype_pos array
 *   Primary SNPs are the two main alleles only; b1 and b2 bases.
 *
 * - Find reads covering SNPs and (if matching b1/b2) and increment
 *   l_b1, l_b2 counters.
 *
 * - Extract primary links from haplotype_pos[] array, nulling the
 *   weak ones. (Forming haplotype-contigs.)
 *
 * FAIL: if we have mix of real SNPs and some error-rich regions, we need
 * to tease out real SNPs from errors, so daisy chaining left to right is
 * not sufficient.   Need NxN linkage matrix instead?
 *
 * Try 2:
 *
 * Foreach SNP of A/B alleles, assign them arbitrarily.
 *
 * Foreach read
 *     Foreach pair of SNP base covered, increment counters of A/B
 *     (4 reads may give 0/4 or 4/0 if all consistent, or something like
 *      1/3 2/2 or 3/1 if inconsistent.)
 *
 * Foreach SNP, compute ratio of SUM(MIN(A,B)) / SUM(MAX(A,B)).  Any SNP which
 *    has, say, > 0.5 ratio can be removed as conflicting. (Needs recomputation
 *    of neighbouring ones then? That makes it O(N^3) instead of O(N^2)?)
 *
 * Cruder try 3:
 *
 * Foreach neighbouring position (not N vs N, just neighbours), count
 * purity of matching reads to prev/next SNP.  Eg bin A/A A/B B/A B/B counts.
 * For correct cases we expect A/A+B/B or A/B+B/A to be high but not the
 * other combinations.  If left is bad and right not or vice versa then
 * maybe it is OK, but if both are bad then we probably have a false SNP.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "misc.h"
#include "find_haplotypes.h"
#include "consensus.h"
#include "align.h"

/*
 * A simple linked list of haplotypic sites, forming a doubly linked list so we can easily remove from it.
 */
typedef struct haplotype_pos {
    char b1, b2;         // bases 1 and 2
    int pos;             // pos
    int score[2];        // score based on links to prev/next.
    int same, opp, mis;  // counts of match/mismatch with next haplotype
    struct haplotype_pos *prev, *next;
} haplotype_pos;

int add_haplotype_pos(haplotype_pos **phead, haplotype_pos **ptail,
		      int pos, int b1, int b2) {
    haplotype_pos *p = calloc(1, sizeof(*p));
    if (!p)
	return -1;

    p->pos = pos;
    p->b1  = b1;
    p->b2  = b2;

    if (*ptail) {
	(*ptail)->next = p;
	p->prev = *ptail;
	*ptail = p;
    } else {
	*phead = *ptail = p;
    }

    return 0;
}

void del_haplotype_pos(haplotype_pos **phead, haplotype_pos **ptail,
		       haplotype_pos *p) {
    if (p == *phead)
	*phead = p->next;
    else
	p->prev->next = p->next;

    if (p == *ptail)
	*ptail = p->prev;
    else
	p->next->prev = p->prev;

    free(p);
}


static int find_haplotypes_single(GapIO *io, tg_rec crec, int start, int end) {
    consensus_t *cons = NULL;
    int ret = -1, i;
    haplotype_pos *phead = NULL, *ptail = NULL;
    int pass;
    rangec_t *rng = NULL;
    int nr;
    contig_t *c;

    // Accumulate a list of haplotypes
    if (!(cons = calloc(end-start+1, sizeof(*cons))))
	goto err;

    if (-1 == calculate_consensus(io, crec, start, end, cons))
	goto err;

    for (i = start; i <= end; i++) {
	if (cons[i-start].scores[6]>0) {
	    printf("Pos %5d: het %c/%c  score %d\n",
		   i,
		   "ACGT*"[cons[i-start].het_call / 5],
		   "ACGT*"[cons[i-start].het_call % 5],
		   (int)cons[i-start].scores[6]);

	    add_haplotype_pos(&phead, &ptail, i,
			      "ACGT*"[cons[i-start].het_call / 5],
			      "ACGT*"[cons[i-start].het_call % 5]);
	}
    }

    c = cache_search(io, GT_Contig, crec);
    if (!c)
	goto err;

    rng = contig_seqs_in_range(io, &c, start, end,
			       CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED, &nr);
    if (!rng)
	goto err;

    pass = 1;
    while (1) {
	int changed = -1;
	printf("\n=== Pass %d ===\n", pass);

	{
	    haplotype_pos *p;

	    for (p = phead; p; p = p->next) {
		p->same = p->opp = p->mis = p->score[0] = p->score[1] = 0;
	    }
	}

	// Accumulate haplotypes
	{
	    rangec_t *r;
	    haplotype_pos *p1, *p2;
	    int i;

	    p1 = phead;
	    for (i = 0; i < nr; i++) {
		rangec_t *r = &rng[i];
		int left, right;
		seq_t *s;
		char b;

		if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		    continue;

		s = cache_search(io, GT_Seq, r->rec);
		if (s->right < s->left)
		    continue; // ERROR: no unclipped bases.

		if ((s->len < 0) ^ r->comp) {
		    left = r->start + ABS(s->len) - (s->right-1) - 1;
		    right = r->start + ABS(s->len) - (s->left-1) - 1;
		} else {
		    left = r->start + s->left - 1;
		    right = r->start + s->right - 1;
		}

		while (p1 && p1->pos < left)
		    p1 = p1->next;
		if (!p1)
		    break;

		if (right < p1->pos)
		    continue;

		if ((s->len < 0) ^ r->comp) {
		    b = complement_base(s->seq[ABS(s->len)-1 - (p1->pos - r->start)]);
		} else {
		    b = s->seq[p1->pos - r->start];
		}
		for (p2 = p1; p2->next && p2->next->pos <= right; p2 = p2->next) {
		    char bn;
		    if ((s->len < 0) ^ r->comp) {
			bn = complement_base(s->seq[ABS(s->len)-1 - (p2->next->pos - r->start)]);
		    } else {
			bn = s->seq[p2->next->pos - r->start];
		    }

		    //printf("Pos %d, rec %c%"PRIrec", base %c(%c)\n",
		    //	   p2->pos, "+-"[(s->len < 0) ^ r->comp],
		    //	   r->rec, b, bn);
		
		    if ((p2->b1 == b && p2->next->b1 == bn) ||
			(p2->b2 == b && p2->next->b2 == bn))
			p2->same++;
		    else if ((p2->b1 == b && p2->next->b2 == bn) ||
			     (p2->b2 == b && p2->next->b1 == bn))
			p2->opp++;
		    else
			p2->mis++;

		    b = bn;
		}
	    }
	}

	// Score haplotypes
	{
	    haplotype_pos *p;
	    int p_score = 0;

	    for (p = phead; p; p = p->next) {
		int count = p->same + p->opp; // + p->mis?
		int score = ABS(p->same - p->opp) * ((2.0 * MAX(p->same, p->opp)) / count -1) -  count/2;

		score = score>0 ? sqrt(100*score) : -sqrt(100*-score);

		p->score[0] = p_score;
		p->score[1] = score;

		printf("Hap %5d %c/%c   score %5d/%5d  count %3d     %3d %3d %3d %c\n",
		       p->pos, p->b1, p->b2,
		       p_score, score, count, p->same, p->opp, p->mis,
		       " .*"[(p->score[0] < 0) + (p->score[1] < 0)]);

		p_score = score;
	    }
	}

	// Cull bad haplotypes
	{
	    haplotype_pos *p, *pn;

	    for (p = phead; p; p = pn) {
		pn = p->next;

		if ((pass == 1 &&  p->score[0] < 0 && p->score[1] < 0) ||
		    (pass == 2 && (p->score[0] < 0 || p->score[1] < 0))) {
		    del_haplotype_pos(&phead, &ptail, p);
		    changed = 1;
		}
	    }
	}

	if (changed == -1) {
	    if (pass == 1)
		pass = 2;
	    else
		break;
	}
    }

    // Dump haplotype
    {
	int allele;
	int al = 0;
	haplotype_pos *p;

	// Allele 2
	for (p = phead; p; p = p->next) {
	    putchar(al ? p->b1 : p->b2);

	    if (al) {
		char tmp = p->b1;
		p->b1 = p->b2;
		p->b2 = tmp;
	    }

	    if (p->score[1] > 0) {
		al = p->same >= p->opp ? al : 1-al;
	    } else {
		putchar(' ');
		al = 0;
	    }
	}
	putchar('\n');

	// Allele 1
	for (p = phead; p; p = p->next) {
	    putchar(p->b1);

	    if (p->score[1] <= 0)
		putchar(' ');
	}
	putchar('\n');
    }


    // Assign reads to haplotypes & build consensus
    if (0) {
	rangec_t *r;
	haplotype_pos *p1, *p2;
	int i;
	rangec_t *r1, *r2;
	int nr1 = 0, nr2 = 0;
	consensus_t *cons_allele;
	char *cons1, *cons2;

	r1 = malloc(nr * sizeof(*r1));
	r2 = malloc(nr * sizeof(*r1));

	cons_allele = malloc((end-start+1) * sizeof(*cons_allele));

	cons1 = malloc(end-start+2);
	cons2 = malloc(end-start+2);

	p1 = phead;
	for (i = 0; i < nr; i++) {
	    rangec_t *r = &rng[i];
	    int left, right;
	    seq_t *s;
	    int a1, a2, an;

	    if ((r->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		continue;

	    s = cache_search(io, GT_Seq, r->rec);
	    if (s->right < s->left)
		continue; // ERROR: no unclipped bases.

	    if ((s->len < 0) ^ r->comp) {
		left = r->start + ABS(s->len) - (s->right-1) - 1;
		right = r->start + ABS(s->len) - (s->left-1) - 1;
	    } else {
		left = r->start + s->left - 1;
		right = r->start + s->right - 1;
	    }

	    while (p1 && p1->pos < left)
		p1 = p1->next;
	    if (!p1)
		break;

	    if (right < p1->pos)
		continue;

	    a1 = a2 = an = 0;
	    for (p2 = p1; p2 && p2->pos <= right; p2 = p2->next) {
		char bn;
		if ((s->len < 0) ^ r->comp) {
		    bn = complement_base(s->seq[ABS(s->len)-1 - (p2->pos - r->start)]);
		} else {
		    bn = s->seq[p2->pos - r->start];
		}
		
		//printf("Pos %d, rec %c%"PRIrec", base %c (%c %c)\n",
		//       p2->pos, "+-"[(s->len < 0) ^ r->comp],
		//       r->rec, bn, p2->b1, p2->b2);
		
		if (p2->b1 == bn)
		    a1++;
		else if (p2->b2 == bn)
		    a2++;
		else
		    an++;
	    }

	    if (a1 / .9 > a1+a2+an) {
		r1[nr1++] = *r;
		//printf("1: %"PRIrec"\n", r->rec);
	    } else if (a2 / .9 > a1+a2+an) {
		r2[nr2++] = *r;
		//printf("2: %"PRIrec"\n", r->rec);
	    } else {
		//printf("?: %"PRIrec"\n", r->rec);
	    }
	}

	calculate_consensus_bit_het(io, crec, start, end, 0, r1, nr1, cons_allele);
	for (i = 0; i <= end-start; i++) {
	    cons1[i] = "ACGT*N"[cons_allele[i].call];
//	    if (cons1[i] == 'N')
//		cons1[i] = "ACGT*N"[cons[i].call];
	}
	cons1[i] = 0;

	calculate_consensus_bit_het(io, crec, start, end, 0, r2, nr2, cons_allele);
	for (i = 0; i <= end-start; i++) {
	    cons2[i] = "ACGT*N"[cons_allele[i].call];
//	    if (cons2[i] == 'N')
//		cons2[i] = "ACGT*N"[cons[i].call];
	}
	cons2[i] = 0;

	printf(">allele1\n%s\n>allele2\n%s\n", cons1, cons2);


	// Align cons1 to cons2 to see if any shift improves things.
	// If so it demonstrates the best way to edit our contigs.
	align_int *S = (align_int *)xcalloc((end-start+1)*2+1, sizeof(*S));
	calign(cons1, cons2, end-start+1, end-start+1,
	       NULL, NULL, NULL, NULL,
	       0, 0, 4, 1,
	       3, // job. Also see  ALIGN_GAP_[ES][12]
	       //3 | ALIGN_GAP_E1 | ALIGN_GAP_E2 | ALIGN_GAP_S1 | ALIGN_GAP_S2,
	       0, S);

	cdisplay(cons1, cons2, end-start+1, end-start+1, 3, S, 1, 1);

	if (*S > 0) {
	    printf("Shift left by %d\n", *S);
	} else if (*S < 0) {
	    printf("Shift right by %d\n", -*S);
	}

	free(S);
	free(cons_allele);
	free(cons1);
	free(cons2);

	free(r1);
	free(r2);
    }

    ret = 0;
 err:
    if (cons)
	free(cons);
    if (rng)
	free(rng);

    return ret;
}

/*
 * Splits readings into haplotypic groups and also returns haplotype consensus?
 * Works via lists? Files?
 *
 * Returns 0 for success
 *        -1 for failure
 */
int find_haplotypes(GapIO *io, contig_list_t *contigs, int ncontigs) {
    int i, err = 0;

    for (i = 0; i < ncontigs; i++) {
	err |= find_haplotypes_single(io, contigs[i].contig,
				      contigs[i].start, contigs[i].end);
    }

    return err;
}
