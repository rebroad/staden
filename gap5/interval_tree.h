/*
 * A very simplistic interval tree. These are similar to the the
 * degenerate case of a 1-dimensional R-Tree.
 *
 * There is no removal system in place yet. The data is not entirely
 * sorted (although defining a sort order on a range is debatable).
 * However the purpose of how to quickly intervals which items are
 * present within a range is quickly solved.
 */

#ifndef _INTERVAL_TREE_H_
#define _INTERVAL_TREE_H_

#include <stdint.h>
#include "tree.h"

// The "payload" of an interval
typedef union {
    uint64_t i;
    void *p;
} interval_data;

// Intervals themselves
typedef struct interval {
    //RB_ENTRY(interval) link;
    struct interval *next;
    struct interval *u_next, *u_prev; // temporary user-defined lists
    int start, end;
    interval_data data;
} interval;

// Opaque data type. See interval_tree.c for implementation details.
typedef struct interval_tree interval_tree;
typedef struct interval_iter interval_iter;

/*
 * Creates a new interval tree.
 *
 * Returns pointer on success;
 *        NULL on failure.
 */
interval_tree *interval_tree_create(void);

/*
 * Frees an interval tree. If free_func is non-NULL then it is
 * called for each interval_data struct held within the tree.
 */
void interval_tree_destroy(interval_tree *it, void (*free_func)(void *ptr));

/*
 * Adds an interval to an interval tree. 
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int interval_tree_add(interval_tree *it, int start, int end, interval_data data);

/*
 * Removes an interval from an interval tree.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int interval_tree_del(interval_tree *it, interval *node);

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
			 void *clientdata);

interval_iter *interval_range_iter(interval_tree *it, int start, int end);
void interval_iter_destroy(interval_iter *iter);
interval *interval_iter_next(interval_iter *iter);

#endif /* _INTERVAL_TREE_H_ */
