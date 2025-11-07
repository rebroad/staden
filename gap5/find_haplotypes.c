/*
 * ---------------------------------------------------------------------------
 * Identifies SNP sites and forms haplotypes from these by linking SNPs
 * together based on spanning readings. 
 *
 * Trim anything that maybe dodgy (quality differs to background?).
 *
 * Iteratively build up haplotypes. Haplotype is [ACGT*-]+.
 * Ie a sequence or unknown (-).
 * A sequence matches a known haplotype if it disagrees with - only.
 * Eg ---ACGTA--- and ----CGTAC-- align and will be merged.
 *
 * For efficiency leading and trailing "-" are run-length encoded.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "misc.h"
#include "find_haplotypes.h"
#include "consensus.h"
#include "align.h"
#include "array.h"
#include "interval_tree.h"
#include "dna_utils.h"

/*
 * This is an attempt to allow <short-match> to be entirely put into <long-match>
 * when it is contained within it during haplotype_str_add().  It reduces the
 * number of things to combine later on, and so speeds up the algorithm.
 *
 * Unfortunately it also gives poorer matches as we do not know what is best
 * to merge until later on.
 */
// #define ALLOW_CONTAINMENTS

/*
 * A string of haplotypic bases. This excludes non-SNP sites
 */
typedef struct haplotype_str {
    struct haplotype_str *next;
    char *snps;          // string of [ACGT*-]
    int *count;          // depth of snps[]
    int nseq;
    // start by ignoring this and use noddy implementation
    int start;           // number of leading "-"
    int end;             // (end-start) = len
    Array recs;
} haplotype_str;

/*
 * Compares 'snps' against the known haplotypes hs[0..n_hs-1].
 *
 * If it matches it incorporates it, otherwise it adds a new
 * haplotype.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int haplotype_str_add(interval_tree *it, char *snps, int start, int end,
		      tg_rec rec1, tg_rec rec2) {
    haplotype_str *tmp;
    interval_iter *iter = interval_range_iter(it, start, end);
    interval *hs, *best_hs = NULL;
    int overlap = 0;
#ifdef ALLOW_CONTAINMENTS
    int best_overlap = 0;
#endif
    int i;

    for (hs = interval_iter_next(iter); hs; hs = interval_iter_next(iter)) {
	haplotype_str *tmp = (haplotype_str *)hs->data.p;
	int i, i_end; // idx to snps
	int j, j_end; // idx to tmp->snps

#ifndef ALLOW_CONTAINMENTS
	if (start != tmp->start || end != tmp->end)
	    continue;
#endif
	// absolute positions
	i     = MAX(tmp->start, start);
	i_end = MIN(tmp->end, end);

	// relative positions to start
	j = MAX(0, i - tmp->start);
	i = MAX(0, i - start);

	j_end = MIN(tmp->end, i_end) - tmp->start;
	i_end = MIN(end,      i_end) - start;

	assert(i_end - i == j_end - j);

	overlap = 0;
	for (; i <= i_end; i++, j++) {
	    assert(snps[i] >= ' ' && snps[i] <= '~');
	    if (tmp->snps[j] != '-' && snps[i] != '-') {
		if (tmp->snps[j] == snps[i]) {
		    overlap++;
		} else {
		    break;
		}
	    }
	}

	if (i != i_end+1)
	    continue; // didn't overlap
#ifdef ALLOW_CONTAINMENTS
	if (best_overlap < overlap) {
	    best_overlap = overlap;
	    best_hs = hs;
	}
#else
	if (tmp->start == start && tmp->end == end) {
	    best_hs = hs;
	    break;
	}
#endif
    }
    interval_iter_destroy(iter);

#ifdef ALLOW_CONTAINMENTS
//    if (best_hs && (best_overlap >= end-start+1 ||
//		    best_overlap == best_hs->end - best_hs->start+1))
//		    // two way containment
    if (best_hs && best_hs->start <= start && best_hs->end >= end)  // larger
#else
    if (best_hs && best_hs->start == start && best_hs->end == end)  // exact match
#endif
    {
	// Overlaps an existing haplotype, so append.
	// NB: no attempt to do joinining made here,
	// but we're processing in left to right order
	// so it is unlikely to be needed.
	hs = best_hs;
	tmp = (haplotype_str *)best_hs->data.p;

	assert(tmp->start <= start);

//        printf("%*s%.*s Update %*s%.*s -> ", start, "", end-start+1, snps,
//	       tmp->start, "", tmp->end - tmp->start + 1, tmp->snps);

#ifdef ALLOW_CONTAINMENTS
	if (tmp->end < end) {
	    interval_tree_del(it, hs);
	    tmp->snps  = realloc(tmp->snps,  end - tmp->start+1);
	    tmp->count = realloc(tmp->count, (end - tmp->start+1)*sizeof(int));
	    memset(&tmp->count[tmp->end-tmp->start+1], 0, (end-tmp->end)*sizeof(int));
	    tmp->end = end;
	    interval_tree_add(it, tmp->start, tmp->end, (interval_data)(void *)tmp);
	}
#endif

	// Consider \0 to be the undef snp and then we can just do
	// tmp->snps[i-tmp->start] |= snps[i-start];
	for (i = start; i <= end; i++) {
	    if (snps[i-start] != '-') {
		tmp->snps[i-tmp->start] = snps[i-start];
		tmp->count[i-tmp->start]++;
	    }
	}

//	printf("%*s%.*s\n", tmp->start, "", tmp->end - tmp->start + 1, tmp->snps);
	tmp->nseq++;

//	// Maintain sorted list by nseq
//	if (last && tmp->nseq > last->nseq) {
//	    if (last_last) {
//		last_last->next = tmp;
//		last->next = tmp->next;
//		tmp->next = last;
//	    } else {
//		*hs = tmp;
//		last->next = tmp->next;
//		tmp->next = last;
//	    }
//	}

	if (rec1)
	    ArrayPush(tmp->recs, tg_rec, rec1);
	if (rec2)
	    ArrayPush(tmp->recs, tg_rec, rec2);

	return 0;
    }


    // Hasn't been merged, so start a new haplotype string
    tmp = calloc(1, sizeof(*tmp));
    tmp->snps = (char *)malloc(end-start+1);
    tmp->count = (int *)calloc(end-start+1, sizeof(int));
    tmp->start = start;
    tmp->end = end;
    tmp->nseq = 1;
    for (i = start; i <= end; i++) {
	if ((tmp->snps[i-start] = snps[i-start]) != '-')
	    tmp->count[i-start] = 1;
    }
//    printf("ADD %*s%.*s\n", tmp->start, "", tmp->end - tmp->start + 1, tmp->snps);
    interval_tree_add(it, start, end, (interval_data)(void *)tmp);

    tmp->recs = ArrayCreate(sizeof(tg_rec), 1);
    if (rec1)
	ArrayPush(tmp->recs, tg_rec, rec1);
    if (rec2)
	ArrayPush(tmp->recs, tg_rec, rec2);

    return 0;
}

void haplotype_str_free(void *vp) {
    haplotype_str *hs = (haplotype_str *)vp;
    if (hs->recs)
	ArrayDestroy(hs->recs);
    if (hs->snps)
	free(hs->snps);
    free(hs);
}

void haplotype_str_filter(interval_tree *it, int min_count) {
    interval_iter *iter = interval_range_iter(it, INT_MIN, INT_MAX);
    interval *iv;
    interval *iv_list = NULL;

    for (iv = interval_iter_next(iter); iv; iv = interval_iter_next(iter)) {
	haplotype_str *hs = (haplotype_str *)iv->data.p;
	if (hs->nseq < min_count) {
	    // Delay deletion as it breaks iterators
	    iv->u_next = iv_list;
	    iv_list = iv;
	}
    }

    while (iv_list) {
	haplotype_str *hs = (haplotype_str *)iv_list->data.p;
	iv = iv_list->u_next;

	interval_tree_del(it, iv_list);
	haplotype_str_free(hs);

	iv_list = iv;
    }

    interval_iter_destroy(iter);
}


// Sort by number of sequences (largest first)
// then start and end coordinates (smallest first).
int ivp_sort(const void *vp1, const void *vp2) {
    haplotype_str *hs1 = (haplotype_str *)(*((interval **)vp1))->data.p;
    haplotype_str *hs2 = (haplotype_str *)(*((interval **)vp2))->data.p;
    //haplotype_str *hs1 = (haplotype_str *)(*ivp1)->data.p;
    //haplotype_str *hs2 = (haplotype_str *)(*ivp2)->data.p;

    int nl1 = sqrt(hs1->end - hs1->start + 1) * hs1->nseq;
    int nl2 = sqrt(hs2->end - hs2->start + 1) * hs2->nseq;

    return (nl2 - nl1)
	 ? (nl2 - nl1)
	: ((hs1->start - hs2->start)
	   ?hs1->start - hs2->start
	   :hs1->end   - hs2->end); 
}

// Clusters a block of haplotypes in intervals head to tail, returning a
// new head/tail.
int haplotype_str_cluster_subregion(interval **head_p, interval **tail_p, int count) {
    interval **ivp, *iv;
    interval *iv_head, *iv_tail;
    interval *iv_prev, *iv_next;
    int i;

    if (count < 1)
	return 0;

    if (!head_p || !*head_p || !tail_p || !*tail_p)
	return -1;

    iv_head = *head_p;
    iv_tail = *tail_p;
    iv_prev = iv_head->u_prev; iv_head->u_prev = NULL;
    iv_next = iv_tail->u_next; iv_tail->u_next = NULL;

//    puts("::sub_start::");
//    for (iv = iv_head; iv; iv = iv->u_next) {
//	haplotype_str *hs = (haplotype_str *)iv->data.p;
//	printf("%4d %*s%.*s\n", hs->nseq, hs->start, "", hs->end - hs->start+1, hs->snps);
//    }

    // Sort it.
    ivp = malloc(count * sizeof(*ivp));
    for (count = 0, iv = iv_head; iv; iv = iv->u_next, count++)
	ivp[count] = iv;
    qsort(ivp, count, sizeof(*ivp), ivp_sort);

    iv_head = ivp[0]; iv_tail = ivp[count-1];
    for (i = 0; i < count; i++) {
	ivp[i]->u_prev = i         ? ivp[i-1] : NULL;
	ivp[i]->u_next = i+1<count ? ivp[i+1] : NULL;
    }

    /*
     * TODO:
     *
     * If we have strings like:
     *
     * 1) AGCTGACAAATGC
     * 2) AGCAAGCAAATGC
     * 3)         AATGCGGTA
     *
     * 3 matches both 1 and 2. In this case, assuming we start with 1,
     * we cannot "recruit" 3 as it could be recruited elsewhere and
     * the elsewhere read is not compatible.  Essentially this has to
     * remain as 3 contigs.
     *
     * Contrast this to:
     *
     * 1)   CTGACAAATGC
     * 2) AGCTGACAAAT
     * 3)         AATGCGGTA
     *
     * Formally:
     * For string S_a, Y_a is the set of {b,...} where S_a and S_b match
     *                 N_a is the set of {b,...} where S_a and S_b mismatch
     *
     * S_a and S_b can be merged if
     *     Y_a intersection N_b is empty &&
     *     Y_b intersection N_a is empty
     *
     * During merge of S_b into S_a:
     *     Y_a becomes Y_a union Y_b
     *     N_a becomes N_a union N_b
     */

    // Recruit overlapping nodes.
    // O(N^2) complexity, hence keep block size small
    for (iv = iv_head; iv; iv = iv->u_next) {
	haplotype_str *hs = (haplotype_str *)iv->data.p;
	interval *iv2, *next;
	int recruited;
	int iv_start = iv->start, iv_end = iv->end;

	//printf("> %4d %*s%.*s\n", hs->nseq, hs->start, "", hs->end - hs->start+1, hs->snps);

    again:
	recruited = 0;
	for (iv2 = iv->u_next; iv2; iv2 = next) {
	    haplotype_str *hs2 = (haplotype_str *)iv2->data.p;
	    int mismatch = 0, nsnp;
	    next = iv2->u_next;

	    if (iv2->start > iv_end || iv2->end < iv_start)
		continue;

	    for (i = MAX(hs->start, hs2->start); i <= MIN(hs->end, hs2->end); i++) {
		if (hs->snps[i-hs->start] != hs2->snps[i-hs2->start] &&
		    hs->snps[i-hs->start] != '-' &&
		    hs2->snps[i-hs2->start] != '-') {
		    mismatch=1;
		    break;
		}
	    }

	    //if (mismatch) printf("! %4d %*s%.*s\n", hs2->nseq, hs2->start, "", hs2->end - hs2->start+1, hs2->snps);

	    if (mismatch)
		continue;

	    //printf("+ %4d %*s%.*s\n", hs2->nseq, hs2->start, "", hs2->end - hs2->start+1, hs2->snps);
	    recruited = 1;

	    // Overlap. Merge hs2 into hs.
	    nsnp = MAX(hs->end, hs2->end) - MIN(hs->start, hs2->start)+1;
	    if (hs2->start <= hs->start) {
		hs2->snps = realloc(hs2->snps, nsnp+1);
		for (i = MAX(hs->start, hs2->start); i <= hs->end; i++) {
		    if (hs->snps[i-hs->start] != '-' || i > hs2->end)
			hs2->snps[i-hs2->start] = hs->snps[i-hs->start];
		}
		hs2->snps[nsnp] = 0;
		free(hs->snps);
		hs->snps = hs2->snps;
	    } else {
		hs->snps = realloc(hs->snps, nsnp+1);
		for (i = MAX(hs->start, hs2->start); i <= hs2->end; i++) {
		    if (hs2->snps[i-hs2->start] != '-' || i > hs->end)
			hs->snps[i-hs->start] = hs2->snps[i-hs2->start];
		}
		hs->snps[nsnp] = 0;
		free(hs2->snps);
	    }

	    hs->nseq += hs2->nseq;
	    hs->start = MIN(hs->start, hs2->start);
	    hs->end   = MAX(hs->end,   hs2->end);

	    hs2->nseq = 0;
	    hs2->snps = NULL;
	    hs2->end = hs2->start-1;

	    // merge arrays
	    ArrayConcat(hs->recs, hs2->recs);
	    ArrayDestroy(hs2->recs);
	    hs2->recs = NULL;

	    // unlink iv2
	    if (iv2->u_prev)
		iv2->u_prev->u_next = iv2->u_next;
	    else
		iv_head = iv2->u_next;

	    if (iv2->u_next)
		iv2->u_next->u_prev = iv2->u_prev;
	    else
		iv_tail = iv2->u_prev;
	}

	// Warning: don't change iv->start and iv->end directly as
	// that is modifying our interval tree structure, rather than
	// the data held within it (hs).
	//
	// Doing this breaks the tree consistency, causing _del to fail
	// in the haplotype_str_filter() function.
	iv_start = hs->start;
	iv_end   = hs->end;

	// If we're recruited more reads into our haplotype_str then try again,
	// as we may now have overlaps that we originally dismissed.  Maybe
	// this needs a limited number of passes.
	if (recruited)
	    goto again;
    }

    // This is a sub-list, so link back the lists either side of
    // the head..tail.
    if (iv_prev) {
	iv_prev->u_next = iv_head;
	iv_head->u_prev = iv_prev;
    }
    if (iv_next) {
	iv_next->u_prev = iv_tail;
	iv_tail->u_next = iv_next;
    }

    *head_p = iv_head;
    *tail_p = iv_tail;

//    puts("::sub_end::");
//    for (iv = iv_head; iv; iv = iv->u_next) {
//	haplotype_str *hs = (haplotype_str *)iv->data.p;
//	printf("%4d %*s%.*s\n", hs->nseq, hs->start, "", hs->end - hs->start+1, hs->snps);
//    }

    free(ivp);

    return 0;
}

// Merge haplotypes with best overlapping cluster
void haplotype_str_cluster(interval_tree *it) {
    interval *iv_head = NULL, *iv_tail = NULL, *iv = NULL;
    interval_iter *iter = interval_range_iter(it, INT_MIN, INT_MAX);
    int count = 0;
    interval *iv_sub_head;
    int longest_haplo = INT_MIN;

    // FIXME: Chunk into blocks of overlapping haplotypes, or maximum counts.
    // Eg 3 blocks:
    // 1111111  2222   333333
    // -----
    //  ------
    //      --
    //          ----
    //                 ----
    //                   ----


    // Produce a linear linked list.
    for (iv = interval_iter_next(iter); iv; iv = interval_iter_next(iter)) {
	if (longest_haplo == INT_MIN) {
	    longest_haplo = iv->end;
	    iv_sub_head = iv;
	} else {
	    if (iv->start > longest_haplo) {
		if (iv_head == iv_sub_head) {
		    // First sub-list
		    haplotype_str_cluster_subregion(&iv_head, &iv_tail, count);
		} else {
		    haplotype_str_cluster_subregion(&iv_sub_head, &iv_tail, count);
		}
		iv_sub_head = iv;
		longest_haplo = iv->end;
		count = 0;
	    } else {
		longest_haplo = MAX(longest_haplo, iv->end);
	    }
	}

	if ((iv->u_prev = iv_tail))
	    iv_tail->u_next = iv;
	else
	    iv_head = iv;
	iv->u_next = NULL;
	iv_tail = iv;
	
	count++;
    }

    interval_iter_destroy(iter);

    if (count == 0)
	return;

    if (iv_head == iv_sub_head)
	haplotype_str_cluster_subregion(&iv_head, &iv_tail, count);
    else
	haplotype_str_cluster_subregion(&iv_sub_head, &iv_tail, count);

    // Dump the merged data
//    puts("--");
//    for (iv = iv_head; iv; iv = iv->u_next) {
//	haplotype_str *hs = (haplotype_str *)iv->data.p;
//	printf("%4d %*s%.*s\n", hs->nseq, hs->start, "", hs->end - hs->start+1, hs->snps);
//    }
}

void haplotype_str_dump(interval_tree *it) {
    interval_iter *iter = interval_range_iter(it, INT_MIN, INT_MAX);
    interval *iv;

    for (iv = interval_iter_next(iter); iv; iv = interval_iter_next(iter)) {
	haplotype_str *hs = (haplotype_str *)iv->data.p;
	// int i;

	if (!hs->nseq)
	    continue;

	printf("%5d %*s%.*s\n",
	       hs->nseq, 
	       hs->start, "",
	       hs->end - hs->start+1, hs->snps);
//	printf("%5d ", hs->nseq);
//	for (i = 0; i < hs->nsnps; i++)
//	    putchar('!'+MIN(90,hs->count[i]));
//	putchar('\n');
    }
    puts("");

    interval_iter_destroy(iter);
}


void haplotype_str_reclist(interval_tree *it, Array rec_list) {
    interval_iter *iter = interval_range_iter(it, INT_MIN, INT_MAX);
    interval *iv;

    for (iv = interval_iter_next(iter); iv; iv = interval_iter_next(iter)) {
	haplotype_str *hs = (haplotype_str *)iv->data.p;
	if (!hs->nseq)
	    continue;

	ArrayPush(rec_list, Array, hs->recs);
	hs->recs = NULL; // avoid it being freed later
    }

    interval_iter_destroy(iter);
}

/*
 * A simple linked list of haplotypic sites, forming a doubly linked list so we can easily remove from it.
 */
typedef struct haplotype_pos {
    int pos;             // pos
    int score;           // FIXME: define this.
    struct haplotype_pos *prev, *next;
} haplotype_pos;

int add_haplotype_pos(haplotype_pos **phead, haplotype_pos **ptail, int pos) {
    haplotype_pos *p = calloc(1, sizeof(*p));
    if (!p)
	return -1;

    p->pos = pos;

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


static int find_haplotypes_single(GapIO *io, tg_rec crec, int start, int end,
				  int min_count, int pairs,
				  float het_score, float discrep_score,
				  Array rec_list) {
    consensus_t *cons = NULL;
    int ret = -1, i;
    haplotype_pos *phead = NULL, *ptail = NULL;
    rangec_t *rng = NULL;
    int nr;
    contig_t *c;
    int nsnps = 0;
    char *hstr = NULL;
    interval_tree *it = NULL;

    // Accumulate a list of haplotypes
    if (!(cons = calloc(end-start+1, sizeof(*cons))))
	goto err;

    if (-1 == calculate_consensus(io, crec, start, end, cons))
	goto err;

    if (!(it = interval_tree_create()))
	goto err;

    for (i = start; i <= end; i++) {
	if (cons[i-start].scores[6]>=het_score ||
	    cons[i-start].discrep>=discrep_score) {
	    printf("Pos %5d: het %c/%c  score %d %f\n",
		   i,
		   "ACGT*"[cons[i-start].het_call / 5],
		   "ACGT*"[cons[i-start].het_call % 5],
		   (int)cons[i-start].scores[6],
		   cons[i-start].discrep);

	    add_haplotype_pos(&phead, &ptail, i);
	    nsnps++;
	}
    }

    hstr = malloc(nsnps+1);

    c = cache_search(io, GT_Contig, crec);
    if (!c)
	goto err;

    rng = contig_seqs_in_range(io, &c, start, end,
			       CSIR_SORT_BY_X | CSIR_SORT_BY_CLIPPED, &nr);
    if (!rng)
	goto err;


    // Pair up the read-pairs, replacing r->pair_rec with -array_index
    // if pair found.
    {
	HashTable *h = HashTableCreate(nr, HASH_DYNAMIC_SIZE);
	if (!h)
	    goto err;

	for (i = 0; i < nr; i++) {
	    HashData hd;
	    HashItem *hi;
	    rangec_t *r = &rng[i];

	    hi = HashTableSearch(h, (char *)&r->pair_rec, sizeof(r->pair_rec));
	    if (hi) {
		rng[hi->data.i].pair_rec = -i;
		HashTableDel(h, hi, 0);
	    } else {
		hd.i = i;
		HashTableAdd(h, (char *)&r->rec, sizeof(r->rec), hd, NULL);
	    }
	}

	HashTableDestroy(h, 0);
    }


    // Accumulate haplotypes
    {
	haplotype_pos *p1, *p2;
	int i;
	int snp_no = 0;

	p1 = phead;
	for (i = 0; i < nr; i++) {
	    rangec_t *r = &rng[i];
	    int left, right;
	    seq_t *s;
	    char b;
	    int snp_no2;

	    // FIXME: optimise, no need to reset all the while.
	    // Fill out bits we need and remember start/end.
	    //memset(hstr, '-', nsnps); 

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
	    left = MAX(left, r->start);
	    right = MIN(right, r->end);

	    while (p1 && p1->pos < left) {
		p1 = p1->next;
		snp_no++;
	    }
	    if (!p1)
		break;

	    if (right < p1->pos)
		continue;

	    snp_no2 = snp_no;
	    for (p2 = p1; p2 && p2->pos <= right; p2 = p2->next) {
		if ((s->len < 0) ^ r->comp) {
		    b = complement_base(s->seq[ABS(s->len)-1 - (p2->pos - r->start)]);
		} else {
		    b = s->seq[p2->pos - r->start];
		}

		hstr[snp_no2++-snp_no] = b;
		assert(b >= ' ' && b <= '~');
	    }
	    hstr[snp_no2-snp_no] = 0;


	    //haplotype_str_add(it, hstr, snp_no, snp_no2-1, r->rec);
	    // Added single which helps to bulk up duplicates and get the 
	    // best haplotypes represented. (Rejected for now)

	    // Now produce pairs. These are more variable in size/distance so
	    // represented with low duplicate count, so will get added near end.
	    // However they can link haplotypes together.
	    //
	    // Do we also need a NOT join set? Ideally should resolve pairs first?
	    //
	    // FIXME brute force. Add hash or similar to speed up
	    {
//		int j;
		rangec_t *rp;

//		for (j = i+1; j < nr; j++) {
//		    rp = &rng[j];
//
//		    if (rp->rec == r->pair_rec && r->rec == rp->pair_rec)
//			break;
//		}
//
//		if (j == nr) {

		if (r->pair_rec >= 0 || !pairs) {
		    // single ended only
		    haplotype_str_add(it, hstr, snp_no, snp_no2-1, r->rec, 0);
		    continue;
		}
		rp = &rng[-r->pair_rec];

		if ((rp->flags & GRANGE_FLAG_ISMASK) != GRANGE_FLAG_ISSEQ)
		    continue;

		s = cache_search(io, GT_Seq, rp->rec);
		if (s->right < s->left)
		    continue; // ERROR: no unclipped bases.

		if ((s->len < 0) ^ rp->comp) {
		    left = rp->start + ABS(s->len) - (s->right-1) - 1;
		    right = rp->start + ABS(s->len) - (s->left-1) - 1;
		} else {
		    left = rp->start + s->left - 1;
		    right = rp->start + s->right - 1;
		}

		// FIXME: if pair overlaps, may want to pick best quality base
		// p2 is next SNP loc beyond the first end 'r'.
		// Fill out p2 to 'rp' first pos with '-'
		while (p2 && p2->pos < left) {
		    hstr[snp_no2++-snp_no] = '-';
		    p2 = p2->next;
		}

		while (p2 && p2->pos <= right) {
		    if ((s->len < 0) ^ rp->comp) {
			b = complement_base(s->seq[ABS(s->len)-1 - (p2->pos - rp->start)]);
		    } else {
			b = s->seq[p2->pos - rp->start];
		    }

		    hstr[snp_no2++-snp_no] = b;
		    p2 = p2->next;
		}
		hstr[snp_no2-snp_no] = 0;
		haplotype_str_add(it, hstr, snp_no, snp_no2-1, r->rec, rp->rec);
		
		rp->flags = GRANGE_FLAG_ISCONS; // prevent use later on
	    }
	}
    }

    //puts("=== After str_add");
    //haplotype_str_dump(it);

    haplotype_str_cluster(it);

    //puts("=== After cluster");
    //haplotype_str_dump(it);

    haplotype_str_filter(it, min_count);

    puts("=== After filter");
    haplotype_str_dump(it);

    haplotype_str_reclist(it, rec_list);

    ret = 0;
 err:

    // Free positions
    {
	haplotype_pos *curr, *next = NULL;
	for (curr = phead; curr; curr = next) {
	    next = curr->next;
	    free(curr);
	}
    }

    // Free tree / strings
    if (it) {
	interval_tree_destroy(it, haplotype_str_free);
    }

    if (cons)
	free(cons);
    if (rng)
	free(rng);
    if (hstr)
	free(hstr);

    return ret;
}

/*
 * Splits readings into haplotypic groups and also returns haplotype consensus?
 * Works via lists? Files?
 *
 * Returns an Array of Arrat of seq record numbers.
 * Returns NULL for failure.a
 */
Array find_haplotypes(GapIO *io, contig_list_t *contigs, int ncontigs,
		      int pairs, float het_score, float discrep_score,
		      int min_count) {
    int i;
    Array rec_list = ArrayCreate(sizeof(Array), 0);

    for (i = 0; i < ncontigs; i++) {
	printf("find_haplotypes =%"PRIrec"\t%d..%d\n",
	       contigs[i].contig, contigs[i].start, contigs[i].end);
	if (-1 == find_haplotypes_single(io, contigs[i].contig,
					 contigs[i].start, contigs[i].end,
					 min_count, pairs,
					 het_score, discrep_score,
					 rec_list)) {
	    // FIXME: free
	    return NULL;
	}
    }

    return rec_list;
}
