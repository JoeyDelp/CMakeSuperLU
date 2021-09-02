/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file ilu_spivotL.c
 * \brief Performs numerical pivoting
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory
 * June 30, 2009
 * </pre>
 */


#include <math.h>
#include <stdlib.h>
#include "superlu/slu_sdefs.h"

#ifndef SGN
#define SGN(x) ((x)>=0?1:-1)
#endif

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *   Performs the numerical pivoting on the current column of L,
 *   and the CDIV operation.
 *
 *   Pivot policy:
 *   (1) Compute thresh = u * max_(i>=j) abs(A_ij);
 *   (2) IF user specifies pivot row k and abs(A_kj) >= thresh THEN
 *	     pivot row = k;
 *	 ELSE IF abs(A_jj) >= thresh THEN
 *	     pivot row = j;
 *	 ELSE
 *	     pivot row = m;
 *
 *   Note: If you absolutely want to use a given pivot order, then set u=0.0.
 *
 *   Return value: 0	  success;
 *		   i > 0  U(i,i) is exactly zero.
 * </pre>
 */

long long
ilu_spivotL(
	const long long  jcol,     /* in */
	const double u,      /* in - diagonal pivoting threshold */
	long long	   *usepr,   /* re-use the pivot sequence given by
			      * perm_r/iperm_r */
	long long	   *perm_r,  /* may be modified */
	long long	   diagind,  /* diagonal of Pc*A*Pc' */
	long long	   *swap,    /* in/out record the row permutation */
	long long	   *iswap,   /* in/out inverse of swap, it is the same as
				perm_r after the factorization */
	long long	   *marker,  /* in */
	long long	   *pivrow,  /* in/out, as an input if *usepr!=0 */
	double	   fill_tol, /* in - fill tolerance of current column
			      * used for a singular column */
	milu_t	   milu,     /* in */
	float	   drop_sum, /* in - computed in ilu_scopy_to_ucol()
                                (MILU only) */
	GlobalLU_t *Glu,     /* modified - global LU data structures */
	SuperLUStat_t *stat  /* output */
       )
{

    long long		 n;	 /* number of columns */
    long long		 fsupc;  /* first column in the supernode */
    long long		 nsupc;  /* no of columns in the supernode */
    long long		 nsupr;  /* no of rows in the supernode */
    long long		 lptr;	 /* points to the starting subscript of the supernode */
    register long long	 pivptr;
    long long		 old_pivptr, diag, ptr0;
    register float  pivmax, rtemp;
    float	 thresh;
    float	 temp;
    float	 *lu_sup_ptr;
    float	 *lu_col_ptr;
    long long		 *lsub_ptr;
    register long long	 isub, icol, k, itemp;
    long long		 *lsub, *xlsub;
    float	 *lusup;
    long long		 *xlusup;
    flops_t	 *ops = stat->ops;
    long long		 info;

    /* Initialize pointers */
    n	       = Glu->n;
    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;
    lusup      = (float *) Glu->lusup;
    xlusup     = Glu->xlusup;
    fsupc      = (Glu->xsup)[(Glu->supno)[jcol]];
    nsupc      = jcol - fsupc;		/* excluding jcol; nsupc >= 0 */
    lptr       = xlsub[fsupc];
    nsupr      = xlsub[fsupc+1] - lptr;
    lu_sup_ptr = &lusup[xlusup[fsupc]]; /* start of the current supernode */
    lu_col_ptr = &lusup[xlusup[jcol]];	/* start of jcol in the supernode */
    lsub_ptr   = &lsub[lptr];	/* start of row indices of the supernode */

    /* Determine the largest abs numerical value for partial pivoting;
       Also search for user-specified pivot, and diagonal element. */
    pivmax = -1.0;
    pivptr = nsupc;
    diag = EMPTY;
    old_pivptr = nsupc;
    ptr0 = EMPTY;
    for (isub = nsupc; isub < nsupr; ++isub) {
        if (marker[lsub_ptr[isub]] > jcol)
            continue; /* do not overlap with a later relaxed supernode */

	switch (milu) {
	    case SMILU_1:
		rtemp = fabs(lu_col_ptr[isub] + drop_sum);
		break;
	    case SMILU_2:
	    case SMILU_3:
                /* In this case, drop_sum contains the sum of the abs. value */
		rtemp = fabs(lu_col_ptr[isub]);
		break;
	    case SILU:
	    default:
		rtemp = fabs(lu_col_ptr[isub]);
		break;
	}
	if (rtemp > pivmax) { pivmax = rtemp; pivptr = isub; }
	if (*usepr && lsub_ptr[isub] == *pivrow) old_pivptr = isub;
	if (lsub_ptr[isub] == diagind) diag = isub;
	if (ptr0 == EMPTY) ptr0 = isub;
    }

    if (milu == SMILU_2 || milu == SMILU_3) pivmax += drop_sum;

    /* Test for singularity */
    if (pivmax < 0.0) {
	fprintf(stderr, "[0]: jcol=%d, SINGULAR!!!\n", jcol);
	fflush(stderr);
	exit(1);
    }
    if ( pivmax == 0.0 ) {
	if (diag != EMPTY)
	    *pivrow = lsub_ptr[pivptr = diag];
	else if (ptr0 != EMPTY)
	    *pivrow = lsub_ptr[pivptr = ptr0];
	else {
	    /* look for the first row which does not
	       belong to any later supernodes */
	    for (icol = jcol; icol < n; icol++)
		if (marker[swap[icol]] <= jcol) break;
	    if (icol >= n) {
		fprintf(stderr, "[1]: jcol=%d, SINGULAR!!!\n", jcol);
		fflush(stderr);
		exit(1);
	    }

	    *pivrow = swap[icol];

	    /* pick up the pivot row */
	    for (isub = nsupc; isub < nsupr; ++isub)
		if ( lsub_ptr[isub] == *pivrow ) { pivptr = isub; break; }
	}
	pivmax = fill_tol;
	lu_col_ptr[pivptr] = pivmax;
	*usepr = 0;
#ifdef DEBUG
	printf("[0] ZERO PIVOT: FILL (%d, %d).\n", *pivrow, jcol);
	fflush(stdout);
#endif
	info =jcol + 1;
    } /* if (*pivrow == 0.0) */
    else {
	thresh = u * pivmax;

	/* Choose appropriate pivotal element by our policy. */
	if ( *usepr ) {
	    switch (milu) {
		case SMILU_1:
		    rtemp = fabs(lu_col_ptr[old_pivptr] + drop_sum);
		    break;
		case SMILU_2:
		case SMILU_3:
		    rtemp = fabs(lu_col_ptr[old_pivptr]) + drop_sum;
		    break;
		case SILU:
		default:
		    rtemp = fabs(lu_col_ptr[old_pivptr]);
		    break;
	    }
	    if ( rtemp != 0.0 && rtemp >= thresh ) pivptr = old_pivptr;
	    else *usepr = 0;
	}
	if ( *usepr == 0 ) {
	    /* Use diagonal pivot? */
	    if ( diag >= 0 ) { /* diagonal exists */
		switch (milu) {
		    case SMILU_1:
			rtemp = fabs(lu_col_ptr[diag] + drop_sum);
			break;
		    case SMILU_2:
		    case SMILU_3:
			rtemp = fabs(lu_col_ptr[diag]) + drop_sum;
			break;
		    case SILU:
		    default:
			rtemp = fabs(lu_col_ptr[diag]);
			break;
		}
		if ( rtemp != 0.0 && rtemp >= thresh ) pivptr = diag;
	    }
	    *pivrow = lsub_ptr[pivptr];
	}
	info = 0;

	/* Reset the diagonal */
	switch (milu) {
	    case SMILU_1:
		lu_col_ptr[pivptr] += drop_sum;
		break;
	    case SMILU_2:
	    case SMILU_3:
		lu_col_ptr[pivptr] += SGN(lu_col_ptr[pivptr]) * drop_sum;
		break;
	    case SILU:
	    default:
		break;
	}

    } /* else */

    /* Record pivot row */
    perm_r[*pivrow] = jcol;
    if (jcol < n - 1) {
	register long long t1, t2, t;
	t1 = iswap[*pivrow]; t2 = jcol;
	if (t1 != t2) {
	    t = swap[t1]; swap[t1] = swap[t2]; swap[t2] = t;
	    t1 = swap[t1]; t2 = t;
	    t = iswap[t1]; iswap[t1] = iswap[t2]; iswap[t2] = t;
	}
    } /* if (jcol < n - 1) */

    /* Interchange row subscripts */
    if ( pivptr != nsupc ) {
	itemp = lsub_ptr[pivptr];
	lsub_ptr[pivptr] = lsub_ptr[nsupc];
	lsub_ptr[nsupc] = itemp;

	/* Interchange numerical values as well, for the whole snode, such 
	 * that L is indexed the same way as A.
	 */
	for (icol = 0; icol <= nsupc; icol++) {
	    itemp = pivptr + icol * nsupr;
	    temp = lu_sup_ptr[itemp];
	    lu_sup_ptr[itemp] = lu_sup_ptr[nsupc + icol*nsupr];
	    lu_sup_ptr[nsupc + icol*nsupr] = temp;
	}
    } /* if */

    /* cdiv operation */
    ops[FACT] += nsupr - nsupc;
    temp = 1.0 / lu_col_ptr[nsupc];
    for (k = nsupc+1; k < nsupr; k++) lu_col_ptr[k] *= temp;

    return info;
}
