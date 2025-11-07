/*
 * An interval tree.  It's a simple alternative to an R-Tree.
 *
 * These can be implemented as a normal one dimensional tree, provided the tree
 * supports having data held in internal nodes as well as leaf nodes. 
 *
 * Each tree node has a single point and houses all data that overlaps that
 * point. The left/right nodes hold data that overlaps points to the left or
 * right of that root node.
 */

// RB_AUGMENT is used to turn the RB tree into an interval tree.
// See http://en.wikipedia.org/wiki/Interval_tree
// and http://www.mjmwired.net/kernel/Documentation/rbtree.txt
#define RB_AUGMENT(x) interval_augment((x))

#include "interval_tree.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

#include "tree.h"

// FIXME: should we just define interval to be a struct
// of a defined but fixed size, with a minimum set of values?
//
// Eg *next, nbytes, start, end.  Then we can add whatever other
// values we want into the struct, without having to have two
// memory blocks and a pointer from interval to the "payload"?
//
// (A poor mans template.)

// The nodes in our tree.
typedef struct interval_node {
    RB_ENTRY(interval_node) link;
    int start;
    int end;
    int last;
    interval *intervals;
} interval_node;

RB_HEAD(interval_t, interval_node);

// The primary tree object
struct interval_tree {
    struct interval_t tree;
};

static void interval_augment(interval_node *n) {
    interval_node *nc;

    //printf("interval_augment %p (%d..%d)\n", n, n->start, n->end);
    int last = n->end;
    n->last = last;

    if ((nc = RB_LEFT(n, link)) && last < nc->last)
	last = nc->last;
    if ((nc = RB_RIGHT(n, link)) && last < nc->last)
	last = nc->last;
    while (n && n->last < last) {
	//printf("Aug %p last=%d\n", n, last);
	n->last = last;
	n = RB_PARENT(n, link);
    }
}

static int interval_cmp(interval_node *n1, interval_node *n2) {
    return (n1->start == n2->start)
	? (n1->end - n2->end)
	: (n1->start - n2->start);
}

RB_PROTOTYPE(interval_t, interval_node, link, interval_cmp);
RB_GENERATE(interval_t, interval_node, link, interval_cmp);

/*
 * Creates a new interval tree.
 *
 * Returns pointer on success;
 *        NULL on failure.
 */
interval_tree *interval_tree_create(void) {
    interval_tree *it = malloc(sizeof(*it));
    if (!it)
	return NULL;

    RB_INIT(&it->tree);
    return it;
}

/*
 * Frees an interval tree. If free_func is non-NULL then it is
 * called for each interval_data struct held within the tree.
 */
void interval_tree_destroy(interval_tree *it, void (*free_func)(void *ptr)) {
    interval_node *next = NULL, *node;

    for (node = RB_MIN(interval_t, &it->tree); node; node = next) {
	next = RB_NEXT(interval_t, &it->tree, node);
	interval *i, *i_next;
	for (i = node->intervals; i; i = i_next) {
	    i_next = i->next;
	    if (free_func)
		free_func(i->data.p);
	    free(i);
	}
	RB_REMOVE(interval_t, &it->tree, node);
	free(node);
    }

    free(it);
}

/*
 * Adds an interval to an interval tree. 
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int interval_tree_add(interval_tree *it, int start, int end, interval_data data) {
    interval *node = malloc(sizeof(*node));
    interval_node *i_node = NULL, n;
    if (!node)
	return -1;

    node->next  = NULL;
    node->start = start;
    node->end   = end;
    node->data  = data;

    // find interval_node that immediately preceeds start.
    n.start  = start;
    n.end    = end;
    n.last   = end;
    i_node = RB_NFIND(interval_t, &it->tree, &n);
    if (i_node && i_node->start > start)
	i_node = RB_PREV(interval_t, &it_tree, i_node);
    //if (!i_node || i_node->start > end) {
    if (!i_node || i_node->start != start) {
	// No overlap, so new node
	i_node = malloc(sizeof(*i_node));
	if (!i_node) {
	    free(node);
	    return -1;
	}
	i_node->start = start;
	i_node->end   = end;
	i_node->intervals = node;
	i_node->last = end;
	RB_INSERT(interval_t, &it->tree, i_node);
    } else {
	// Overlaps, so add to intervals packed in this node.
	node->next = i_node->intervals;
	i_node->intervals = node;
	if (i_node->end < end)
	    i_node->end = end;

	while (i_node && i_node->last < end) {
	    //printf("aug %p last=%d\n", i_node, end);
	    i_node->last = end;
	    i_node = RB_PARENT(i_node, link);
	}
    }

    return 0;
}

/*
 * Updates n->end following a deletion to this node and
 * updates n->last on this node and any parent nodes that
 * may need changing.
 */
#define MAX(a,b) ((a)>(b)?(a):(b))
static void interval_recompute_last_del(interval_node *n, int end) {
    int last = INT_MIN;
    interval *i;
    interval_node *nc;

    // Removing the end item means we must recompute to see if it shrunk.
    if (n->end == end) {
	for (i = n->intervals; i; i = i->next)
	    if (last < i->end)
		last = i->end;

	n->end = last;
    } else {
	last = n->end;
    }

    if ((nc = RB_LEFT(n, link)) && last < nc->last)
	last = nc->last;
    if ((nc = RB_RIGHT(n, link)) && last < nc->last)
	last = nc->last;

    n->last = last;
    
    for (n = RB_PARENT(n, link); n; n = RB_PARENT(n, link)) {
	last = n->end;
	if ((nc = RB_LEFT(n, link)) && last < nc->last)
	    last = nc->last;
	if ((nc = RB_RIGHT(n, link)) && last < nc->last)
	    last = nc->last;

	n->last = last;
    }
}

/*
 * Removes an interval from an interval tree.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
static int counter = 0;
int interval_tree_del(interval_tree *it, interval *node) {
    interval_node *i_node, in;
    interval *n, *l = NULL;

    counter++;

    // Get the node containing this interval
    in.start = node->start;
    in.end   = node->end;
    in.last  = node->end;
    i_node = RB_NFIND(interval_t, &it->tree, &in);
    if (!i_node)
	return -1;

    // Find the interval itself
    for (n = i_node->intervals; n; l = n, n = n->next) {
	if (n == node)
	    break;
    }
    if (!n)
	return -1;

    // Remove it
    if (l) {
	l->next = n->next;
    } else {
	i_node->intervals = n->next;
    }

    if (!i_node->intervals) {
	//printf("Removing i_node %p\n", i_node);
	RB_REMOVE(interval_t, &it->tree, i_node);

	if (RB_PARENT(i_node, link))
	    interval_recompute_last_del(RB_PARENT(i_node, link), node->end);

	free(i_node);
    } else {
	interval_recompute_last_del(i_node, node->end);
    }

    free(n);

    return 0;
}


static int interval_recurse(interval_node *i_node, int start, int end,
			    int (*func)(interval *i, void *cd),
			    void *clientdata) {
    interval *i;
    interval_node *n;
    int count = 0, tmp;

    if (!i_node)
	return -1;

    if ((n = RB_LEFT(i_node, link))) {
	//printf("L: %d..%d\n", n->start, n->last);
	if (n->last >= start) {
	    tmp = interval_recurse(n, start, end, func, clientdata);
	    if (tmp < 0)
		return -1;
	    count += tmp;
	}
    }

    //printf("Node %p, %d..%d\n", i_node, i_node->start, i_node->last);
    if (end >= i_node->start && start <= i_node->end) {
	for (i = i_node->intervals; i; i = i->next) {
	    if (i->start <= end && i->end >= start) {
		//printf("%d..%d\n", i->start, i->end);
		count++;
		if (func) {
		    tmp = func(i, clientdata);
		    if (tmp < 0)
			return -1;
		    if (tmp == 0)
			return count;
		}
	    }
	}
    }

    if (i_node->start <= end && (n = RB_RIGHT(i_node, link))) {
	//printf("R: %d..%d\n", n->start, n->last);
	tmp = interval_recurse(n, start, end, func, clientdata);
	if (tmp < 0)
	    return -1;
	count += tmp;
    }

    return count;
}

/*
 * Finds a set of intervals overlapping a given range.
 * The resulting intervals are returned by executing a
 * callback function per element, if non-NULL.  The clientdata parameter is
 * passed through to the callback function as-is.
 *
 * The callback function return values are:
 *     <0 error
 *      0 do not generate more callbacks for this query
 *     >0 keep generating callbacks, if more intervals are found.
 * 
 *
 * Returns the number of items found on success.
 *         -1 on failure.
 */
int interval_range_query(interval_tree *it, int start, int end,
			 int (*func)(interval *i, void *cd),
			 void *clientdata) {
    return interval_recurse(RB_ROOT(&it->tree), start, end,
			    func, clientdata);
}

static void interval_tree_dump_(interval_node *n, int verbosity, int indent) {
    interval *i;
    int c = 0, min = INT_MAX, max = INT_MIN;
    for (i = n->intervals; i; i = i->next) {
	c++;
	if (min > i->start) min = i->start;
	if (max < i->end)   max = i->end;
    }

    printf("%*sNode %p, %d..%d, last %d, range %d..%d, count %d\n",
	   indent, "", n, n->start, n->end, n->last, min, max, c);

    assert(min == n->start);
    assert(max == n->end);
    assert(n->last >= n->end);

    if (verbosity) {
	for (i = n->intervals; i; i = i->next) {
	    printf("%*sInterval %p %d..%d\n", indent, "", i, i->start, i->end);
	}
    }

    if (RB_LEFT(n, link))
	interval_tree_dump_(RB_LEFT(n, link), verbosity, indent+2);
    if (RB_RIGHT(n, link))
	interval_tree_dump_(RB_RIGHT(n, link), verbosity, indent+2);
}

void interval_tree_dump(interval_tree *it, int verbosity) {
    interval_tree_dump_(RB_ROOT(&it->tree), verbosity, 0);
    puts("");
}


int interval_tree_check_(interval_node *n, int *end_p) {
    int err = 0;
    int start = INT_MAX, end = INT_MIN, end_l = INT_MIN, end_r = INT_MIN;
    interval *i;

    if (!n)
	return 0;

    for (i = n->intervals; i; i = i->next) {
	if (start > i->start)
	    start = i->start;
	if (end   < i->end)
	    end   = i->end;
    }

    if (start != n->start || end != n->end) {
	fprintf(stderr, "CHECK node %p: start/end mismatch\n", n);
	err |= 1;
    }

    if (RB_LEFT(n, link))
	err |= interval_tree_check_(RB_LEFT(n, link), &end_l);
    if (end < end_l)
	end = end_l;

    if (RB_RIGHT(n, link))
	err |= interval_tree_check_(RB_RIGHT(n, link), &end_r);
    if (end < end_r)
	end = end_r;

    if (end != n->last) {
	fprintf(stderr, "CHECK node %p: last mismatch\n", n);
	err |= 1;
    }

    if (end_p)
	*end_p = end;

    return err;
}

/*
 * Recursively checks an interval tree to ensure that all the internal
 * requirements hold true.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int interval_tree_check(interval_tree *it) {
    return interval_tree_check_(RB_ROOT(&it->tree), NULL);
}


struct interval_iter {
    interval_tree *tree;
    interval_node *node;
    interval      *iv;
    int start, end, done_lr;
};

interval_iter *interval_range_iter(interval_tree *it, int start, int end) {
    interval_iter *iter;

    if (!(iter = malloc(sizeof(*iter))))
	return NULL;
    
    iter->tree      = it;
    iter->node      = RB_ROOT(&it->tree);
    iter->iv        = iter->node ? iter->node->intervals : NULL;
    iter->done_lr   = 0;
    iter->start     = start;
    iter->end       = end;

    //printf("N: %p %d..%d\n", iter->node, iter->node->start, iter->node->last);

    return iter;
}

void interval_iter_destroy(interval_iter *iter) {
    if (iter)
	free(iter);
}

interval *interval_iter_next_old(interval_iter *iter) {
    interval *i;
    interval_node *n;

    if (!iter->node)
	return NULL;

    //printf("N: %p %d..%d\n", iter->node, iter->node->start, iter->node->last);

    // Current interval array
    for (i = iter->iv; i; i = i->next) {
	//printf("i: %p %d..%d\n", i, i->start, i->end);
	if (i->start <= iter->end && i->end >= iter->start) {
	    iter->iv = i->next;
	    return i;
	}
    }


    // Next node, left if able.
    if ((n = RB_LEFT(iter->node, link))) {
	if (n->last >= iter->start) {
	    //printf("L: %p %d..%d\n", n, n->start, n->last);
	    iter->node = n;
	    if (iter->end >= iter->node->start && iter->start <= iter->node->end)
		iter->iv = iter->node->intervals;
	    else
		iter->iv = NULL;
	    return interval_iter_next(iter);
	}
    }

    while (iter->node) {
	// then right, or up-right.
	//printf("iter->node = %p\n", iter->node);
	n = RB_RIGHT(iter->node, link);
	if (iter->node->start <= iter->end && n) {
	    //printf("R: %p %d..%d\n", n, n->start, n->last);
	    iter->node = n;
	    if (iter->end >= iter->node->start && iter->start <= iter->node->end)
		iter->iv = iter->node->intervals;
	    else
		iter->iv = NULL;
	    return interval_iter_next(iter);
	}

	do {
	    n = iter->node;
	    iter->node = RB_PARENT(iter->node, link);
	    //printf("U: %p\n", iter->node);
	} while (iter->node && RB_RIGHT(iter->node, link) == n);
    }

    return NULL;
}

interval *interval_iter_next(interval_iter *iter) {
    interval *i;
    interval_node *n;

    if (!iter->node)
	return NULL;

    //printf("N %p\n", iter->node);

    // Next node, left if able.
    if (!iter->done_lr && (n = RB_LEFT(iter->node, link))) {
	if (n->last >= iter->start) {
	    //printf("L: %p %d..%d\n", n, n->start, n->last);
	    iter->node = n;
	    if (iter->end >= iter->node->start && iter->start <= iter->node->end)
		iter->iv = iter->node->intervals;
	    else
		iter->iv = NULL;
	    iter->done_lr = 0;
    	    return interval_iter_next(iter);
	}
    }

    iter->done_lr = 1;

    while (iter->node) {
	// Current interval array
	for (i = iter->iv; i; i = i->next) {
	    //printf("i: %p %d..%d\n", i, i->start, i->end);
	    if (i->start <= iter->end && i->end >= iter->start) {
		iter->iv = i->next;
		return i;
	    }
	}

	// then right, or up-right.
	//printf("iter->node = %p\n", iter->node);
	n = RB_RIGHT(iter->node, link);
	if (iter->node->start <= iter->end && n) {
	    //printf("R: %p %d..%d\n", n, n->start, n->last);
	    iter->node = n;
	    if (iter->end >= iter->node->start && iter->start <= iter->node->end)
		iter->iv = iter->node->intervals;
	    else
		iter->iv = NULL;
	    iter->done_lr = 0;
	    return interval_iter_next(iter);
	}

	do {
	    n = iter->node;
	    iter->node = RB_PARENT(iter->node, link);
	    //printf("U: %p\n", iter->node);
	} while (iter->node && RB_RIGHT(iter->node, link) == n);

	if (iter->node) {
	    if (iter->end >= iter->node->start && iter->start <= iter->node->end)
		iter->iv = iter->node->intervals;
	    else
		iter->iv = NULL;
	}

	iter->done_lr = 1;
    }

    return NULL;
}

#ifdef TEST_MAIN
#define NITEMS 100000
#define RLEN   10000000
#define SLEN   10

static interval intervals[NITEMS];

int brute_force_count(int start, int end) {
    int i, c;
    for (i = c = 0; i < NITEMS; i++)
	if (intervals[i].start <= end && intervals[i].end >= start)
	    c++;

    return c;
}

int print_result(interval *i, void *prefix) {
    printf("%s%d..%d\n", (char *)prefix, i->start, i->end);
    return 1;
}

int main(int argc, char **argv) {
    interval_tree *it = interval_tree_create();
    int i, st, en, count1, count2 = 0, tc = 0;

    st = argc>1?atoi(argv[1]):500;
    en = argc>2?atoi(argv[2]):st;

    srand(0);
    for (i = 0; i < NITEMS; i++) {
	interval_data d = {0};
        int x1 = drand48()*RLEN;
        int x2 = x1 + drand48()*SLEN;
        //printf("Adding %d..%d\n", x1, x2);
	interval_tree_add(it, x1, x2, d);

	// For cross-checking
	intervals[i].start = x1;
	intervals[i].end   = x2;
	//interval_tree_dump(it, 0);
    }

    //interval_tree_dump(it, 1);

    if (0) {
	interval_range_query(it, INT_MIN, INT_MAX, print_result, "> ");
    }

    if (0) {
	interval_iter *iter = interval_range_iter(it, INT_MIN, INT_MAX);
	interval *iv;
	while (iv = interval_iter_next(iter))
	    printf("> %d..%d\n", iv->start, iv->end);
    }
    
    //interval_tree_dump(it, 1);
    //exit(0);

    puts("Checking");
    interval_tree_check(it);
    puts("Done");

    count1 = interval_range_query(it, INT_MIN, INT_MAX, NULL, NULL);

    for (i = 0; i < 10000; i++) {
	st = drand48()*(RLEN-SLEN*10);
	en = st + (drand48()*(SLEN*10));

	printf("%d: %d..%d\n", i, st, en);

	count1 = brute_force_count(st, en);

	//count1 = interval_range_query(it, st, en, NULL, NULL);
	//interval_iter *iter = interval_range_iter(it, st, en);
	//for (count1 = 0; interval_iter_next(iter); count1++);

	//count2 = interval_range_query(it, st, en, print_result, "> ");
	count2 = interval_range_query(it, st, en, NULL, NULL);
	////printf("%d..%d: %d found (brute force=%d)\n", st, en, count2, count1);
	assert(count1 == count2);

	interval_iter *iter = interval_range_iter(it, st, en);
	//interval *iv; for (count2 = 0; iv = interval_iter_next(iter); count2++) printf("iv %d..%d\n", iv->start, iv->end);
	for (count2 = 0; interval_iter_next(iter); count2++);
	//printf("%d..%d: %d found (brute force=%d)\n", st, en, count2, count1);
	assert(count1 == count2);
	interval_iter_destroy(iter);

//	iter = interval_range_iter(it, st, en);
//	interval *iv = interval_iter_next(iter);
//	interval *iv_list = NULL;
//	while (iv) {
//	    //printf("Del %d..%d\n", iv->start, iv->end);
//	    // Delay deletion as it breaks iterators.
//	    iv->u_next = iv_list;
//	    iv_list = iv;
//	    iv = interval_iter_next(iter);
//	}
//	interval_iter_destroy(iter);
//
//	while (iv_list) {
//	    iv = iv_list->u_next;
//	    if (0 != interval_tree_del(it, iv_list))
//		abort();
//	    iv_list = iv;
//	    count1--;
//	}
//
//	count2 = interval_range_query(it, INT_MIN, INT_MAX, NULL, NULL);
//	assert(count1 == count2);
//
//	interval_tree_check(it);

	tc += count1;
    }

    interval_tree_destroy(it, NULL);

    printf("Total total = %d\n", tc);

    return 0;
}
#endif
