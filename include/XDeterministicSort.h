#ifndef XDETERMINISTICSORT_H
#define XDETERMINISTICSORT_H

/**
 * @file XDeterministicSort.h
 * @brief A platform-independent unstable quicksort with the exact element
 * order of the Microsoft C runtime's qsort().
 *
 * The engine sorts arrays whose keys are not unique (for example the level's
 * physics object table sorted by sector).  qsort() leaves elements with equal
 * keys in an implementation-defined order, and every C library orders them
 * differently: the Microsoft runtime uses a median-of-three quicksort, glibc a
 * mergesort or introsort, bionic and musl their own smoothsort/quicksort.
 * Game logic that iterates such arrays (object creation order, physics body
 * order) then diverges between platforms even though the sort is "correct"
 * everywhere.
 *
 * The retail game runs on the Microsoft runtime, so this is the reference
 * order.  XDeterministicQSort() reproduces the Microsoft algorithm step for
 * step - median-of-three partitioning around the middle element, a cutoff of
 * eight elements below which a selection sort ("shortsort") is used, and
 * explicit recursion on the smaller partition - so that every platform yields
 * the same permutation for the same input and comparison function.  It is
 * verified against the Microsoft runtime by an exhaustive fuzz test in
 * BallanceMMO (scripts/test_det_qsort.cpp).
 *
 * The comparison function receives pointers to elements, like qsort().
 */

#include <stddef.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define XDETSORT_CUTOFF 8
#define XDETSORT_STKSIZ (8 * sizeof(void *) - 2)

typedef int (*XDeterministicSortFunc)(const void *, const void *);

static void XDeterministicSortSwap(char *a, char *b, size_t width) {
    if (a == b) return;
    while (width--) {
        char tmp = *a;
        *a++ = *b;
        *b++ = tmp;
    }
}

/* Selection sort of [lo, hi] (inclusive), moving the largest element to the
 * end of the range on every pass, exactly like the Microsoft "shortsort". */
static void XDeterministicShortSort(char *lo, char *hi, size_t width, XDeterministicSortFunc comp) {
    char *p, *max;
    while (hi > lo) {
        max = lo;
        for (p = lo + width; p <= hi; p += width) {
            if (comp(p, max) > 0) max = p;
        }
        XDeterministicSortSwap(max, hi, width);
        hi -= width;
    }
}

static void XDeterministicQSort(void *base, size_t num, size_t width, XDeterministicSortFunc comp) {
    char *lo, *hi;        /* ends of the sub-array currently sorting */
    char *mid;            /* points to the middle of the sub-array */
    char *loguy, *higuy;  /* traveling pointers for the partition step */
    size_t size;          /* size of the sub-array */
    char *lostk[XDETSORT_STKSIZ], *histk[XDETSORT_STKSIZ];
    int stkptr;           /* stack for saving sub-array to be processed */

    if (num < 2 || width == 0 || !comp) return;

    stkptr = 0;
    lo = (char *)base;
    hi = (char *)base + width * (num - 1);

recurse:
    size = (size_t)(hi - lo) / width + 1;

    if (size <= XDETSORT_CUTOFF) {
        XDeterministicShortSort(lo, hi, width, comp);
    } else {
        /* First we pick a partitioning element.  The efficiency of the
         * algorithm demands that we find one that is approximately the median
         * of the values, but also that we select one fast.  We choose the
         * median of the first, middle, and last elements, to avoid bad
         * performance in the face of already sorted data, or data that is made
         * up of multiple sorted runs appended together.  Testing shows that a
         * median-of-three algorithm provides better performance than simply
         * picking the middle element for the latter case. */
        mid = lo + (size / 2) * width;

        /* Sort the first, middle, last elements into order. */
        if (comp(lo, mid) > 0) XDeterministicSortSwap(lo, mid, width);
        if (comp(lo, hi) > 0) XDeterministicSortSwap(lo, hi, width);
        if (comp(mid, hi) > 0) XDeterministicSortSwap(mid, hi, width);

        /* We now wish to partition the array into three pieces, one consisting
         * of elements <= partition element, one of elements equal to the
         * partition element, and one of elements > than it. */
        loguy = lo;
        higuy = hi;

        for (;;) {
            if (mid > loguy) {
                do {
                    loguy += width;
                } while (loguy < mid && comp(loguy, mid) <= 0);
            }
            if (mid <= loguy) {
                do {
                    loguy += width;
                } while (loguy <= hi && comp(loguy, mid) <= 0);
            }

            do {
                higuy -= width;
            } while (higuy > mid && comp(higuy, mid) > 0);

            if (higuy < loguy) break;

            XDeterministicSortSwap(loguy, higuy, width);

            /* If the partition element was moved, follow it.  Only need to
             * check for mid == higuy, since before the swap, A[loguy] > A[mid]
             * implies loguy != mid. */
            if (mid == higuy) mid = loguy;
        }

        /* Find adjacent elements equal to the partition element.  The doubled
         * loop is to avoid calling comp(mid, mid), since some existing
         * comparison functions don't work when passed the same value for both
         * pointers. */
        higuy += width;
        if (mid < higuy) {
            do {
                higuy -= width;
            } while (higuy > mid && comp(higuy, mid) == 0);
        }
        if (mid >= higuy) {
            do {
                higuy -= width;
            } while (higuy > lo && comp(higuy, mid) == 0);
        }

        /* We've finished the partition, now we want to sort the subarrays
         * [lo, higuy] and [loguy, hi].  We do the smaller one first to
         * minimize stack usage.  We only sort arrays of length 2 or more. */
        if (higuy - lo >= hi - loguy) {
            if (lo < higuy) {
                lostk[stkptr] = lo;
                histk[stkptr] = higuy;
                ++stkptr;
            }
            if (loguy < hi) {
                lo = loguy;
                goto recurse;
            }
        } else {
            if (loguy < hi) {
                lostk[stkptr] = loguy;
                histk[stkptr] = hi;
                ++stkptr;
            }
            if (lo < higuy) {
                hi = higuy;
                goto recurse;
            }
        }
    }

    /* We have sorted the array, except for any pending sorts on the stack.
     * Check if there are any, and do them. */
    --stkptr;
    if (stkptr >= 0) {
        lo = lostk[stkptr];
        hi = histk[stkptr];
        goto recurse;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* XDETERMINISTICSORT_H */
