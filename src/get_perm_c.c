/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file get_perm_c.c
 * \brief Matrix permutation operations
 *
 * <pre>
 * -- SuperLU routine (version 3.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 1, 2008
 * </pre>
 */
#include "superlu/slu_ddefs.h"
#include "superlu/colamd.h"

extern long long  genmmd_(long long *, long long *, long long *, long long *, long long *, long long *, long long *, 
		    long long *, long long *, long long *, long long *, long long *);

void
get_colamd(
	   const long long m,  /* number of rows in matrix A. */
	   const long long n,  /* number of columns in matrix A. */
	   const long long nnz,/* number of nonzeros in matrix A. */
	   long long *colptr,  /* column pointer of size n+1 for matrix A. */
	   long long *rowind,  /* row indices of size nz for matrix A. */
	   long long *perm_c   /* out - the column permutation vector. */
	   )
{
    long long Alen, *A, i, info, *p;
    double knobs[COLAMD_KNOBS];
    long long stats[COLAMD_STATS];

    Alen = colamd_recommended(nnz, m, n);

    colamd_set_defaults(knobs);

    if (!(A = (long long *) SUPERLU_MALLOC(Alen * sizeof(long long))) )
        ABORT("Malloc fails for A[]");
    if (!(p = (long long *) SUPERLU_MALLOC((n+1) * sizeof(long long))) )
        ABORT("Malloc fails for p[]");
    for (i = 0; i <= n; ++i) p[i] = colptr[i];
    for (i = 0; i < nnz; ++i) A[i] = rowind[i];
    info = colamd(m, n, Alen, A, p, knobs, stats);
    if ( info == FALSE ) ABORT("COLAMD failed");

    for (i = 0; i < n; ++i) perm_c[p[i]] = i;

    SUPERLU_FREE(A);
    SUPERLU_FREE(p);
}
/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * Form the structure of A'*A. A is an m-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'*A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (ata_colptr, ata_rowind).
 *
 * This routine is modified from GETATA routine by Tim Davis.
 * The complexity of this algorithm is: SUM_{i=1,m} r(i)^2,
 * i.e., the sum of the square of the row counts.
 *
 * Questions
 * =========
 *     o  Do I need to withhold the *dense* rows?
 *     o  How do I know the number of nonzeros in A'*A?
 * </pre>
 */
void
getata(
       const long long m,      /* number of rows in matrix A. */
       const long long n,      /* number of columns in matrix A. */
       const long long nz,     /* number of nonzeros in matrix A */
       long long *colptr,      /* column pointer of size n+1 for matrix A. */
       long long *rowind,      /* row indices of size nz for matrix A. */
       long long *atanz,       /* out - on exit, returns the actual number of
                            nonzeros in matrix A'*A. */
       long long **ata_colptr, /* out - size n+1 */
       long long **ata_rowind  /* out - size *atanz */
       )
{
    register long long i, j, k, col, num_nz, ti, trow;
    long long *marker, *b_colptr, *b_rowind;
    long long *t_colptr, *t_rowind; /* a column oriented form of T = A' */

    if ( !(marker = (long long*) SUPERLU_MALLOC((SUPERLU_MAX(m,n)+1)*sizeof(long long))) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (long long*) SUPERLU_MALLOC((m+1) * sizeof(long long))) )
	ABORT("SUPERLU_MALLOC t_colptr[]");
    if ( !(t_rowind = (long long*) SUPERLU_MALLOC(nz * sizeof(long long))) )
	ABORT("SUPERLU_MALLOC fails for t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < m; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < m; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}

    
    /* ----------------------------------------------------------------
       compute B = T * A, where column j of B is:

       Struct (B_*j) =    UNION   ( Struct (T_*k) )
                        A_kj != 0

       do not include the diagonal entry
   
       ( Partition A as: A = (A_*1, ..., A_*n)
         Then B = T * A = (T * A_*1, ..., T * A_*n), where
         T * A_*j = (T_*1, ..., T_*m) * A_*j.  )
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    num_nz++;
		}
	    }
	}
    }
    *atanz = num_nz;
    
    /* Allocate storage for A'*A */
    if ( !(*ata_colptr = (long long*) SUPERLU_MALLOC( (n+1) * sizeof(long long)) ) )
	ABORT("SUPERLU_MALLOC fails for ata_colptr[]");
    if ( *atanz ) {
	if ( !(*ata_rowind = (long long*) SUPERLU_MALLOC( *atanz * sizeof(long long)) ) )
	    ABORT("SUPERLU_MALLOC fails for ata_rowind[]");
    }
    b_colptr = *ata_colptr; /* aliasing */
    b_rowind = *ata_rowind;
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	b_colptr[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    b_rowind[num_nz++] = trow;
		}
	    }
	}
    }
    b_colptr[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
}


/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'+A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (b_colptr, b_rowind).
 * </pre>
 */
void
at_plus_a(
	  const long long n,      /* number of columns in matrix A. */
	  const long long nz,     /* number of nonzeros in matrix A */
	  long long *colptr,      /* column pointer of size n+1 for matrix A. */
	  long long *rowind,      /* row indices of size nz for matrix A. */
	  long long *bnz,         /* out - on exit, returns the actual number of
                               nonzeros in matrix A'*A. */
	  long long **b_colptr,   /* out - size n+1 */
	  long long **b_rowind    /* out - size *bnz */
	  )
{
    register long long i, j, k, col, num_nz;
    long long *t_colptr, *t_rowind; /* a column oriented form of T = A' */
    long long *marker;

    if ( !(marker = (long long*) SUPERLU_MALLOC( n * sizeof(long long)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (long long*) SUPERLU_MALLOC( (n+1) * sizeof(long long)) ) )
	ABORT("SUPERLU_MALLOC fails for t_colptr[]");
    if ( !(t_rowind = (long long*) SUPERLU_MALLOC( nz * sizeof(long long)) ) )
	ABORT("SUPERLU_MALLOC fails t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}


    /* ----------------------------------------------------------------
       compute B = A + T, where column j of B is:

       Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

       do not include the diagonal entry
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}
    }
    *bnz = num_nz;
    
    /* Allocate storage for A+A' */
    if ( !(*b_colptr = (long long*) SUPERLU_MALLOC( (n+1) * sizeof(long long)) ) )
	ABORT("SUPERLU_MALLOC fails for b_colptr[]");
    if ( *bnz) {
      if ( !(*b_rowind = (long long*) SUPERLU_MALLOC( *bnz * sizeof(long long)) ) )
	ABORT("SUPERLU_MALLOC fails for b_rowind[]");
    }
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	(*b_colptr)[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}
    }
    (*b_colptr)[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
}

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * GET_PERM_C obtains a permutation matrix Pc, by applying the multiple
 * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A'.
 * or using approximate minimum degree column ordering by Davis et. al.
 * The LU factorization of A*Pc tends to have less fill than the LU 
 * factorization of A.
 *
 * Arguments
 * =========
 *
 * ispec   (input) long long
 *         Specifies the type of column ordering to reduce fill:
 *         = 1: minimum degree on the structure of A^T * A
 *         = 2: minimum degree on the structure of A^T + A
 *         = 3: approximate minimum degree for unsymmetric matrices
 *         If ispec == 0, the natural ordering (i.e., Pc = I) is returned.
 * 
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A 
 *         can be: Stype = NC; Dtype = _D; Mtype = GE. In the future,
 *         more general A can be handled.
 *
 * perm_c  (output) long long*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 * </pre>
 */
void
get_perm_c(long long ispec, SuperMatrix *A, long long *perm_c)
{
    NCformat *Astore = A->Store;
    long long m, n, bnz = 0, *b_colptr, i;
    long long delta, maxint, nofsub, *invp;
    long long *b_rowind, *dhead, *qsize, *llist, *marker;
    double t, SuperLU_timer_();
    
    m = A->nrow;
    n = A->ncol;

    t = SuperLU_timer_();
    switch ( ispec ) {
    case (NATURAL): /* Natural ordering */
	for (i = 0; i < n; ++i) perm_c[i] = i;
#if ( PRNTlevel>=1 )
	printf("Use natural column ordering.\n");
#endif
	return;
    case (MMD_ATA): /* Minimum degree ordering on A'*A */
	getata(m, n, Astore->nnz, Astore->colptr, Astore->rowind,
		     &bnz, &b_colptr, &b_rowind);
#if ( PRNTlevel>=1 )
	printf("Use minimum degree ordering on A'*A.\n");
#endif
	t = SuperLU_timer_() - t;
	/*printf("Form A'*A time = %8.3f\n", t);*/
	break;
    case (MMD_AT_PLUS_A): /* Minimum degree ordering on A'+A */
	if ( m != n ) ABORT("Matrix is not square");
	at_plus_a(n, Astore->nnz, Astore->colptr, Astore->rowind,
		  &bnz, &b_colptr, &b_rowind);
#if ( PRNTlevel>=1 )
	printf("Use minimum degree ordering on A'+A.\n");
#endif
	t = SuperLU_timer_() - t;
	/*printf("Form A'+A time = %8.3f\n", t);*/
	break;
    case (COLAMD): /* Approximate minimum degree column ordering. */
	get_colamd(m, n, Astore->nnz, Astore->colptr, Astore->rowind, perm_c);
#if ( PRNTlevel>=1 )
	printf(".. Use approximate minimum degree column ordering.\n");
#endif
	return; 
    default:
	ABORT("Invalid ISPEC");
    }

    if ( bnz != 0 ) {
	t = SuperLU_timer_();

	/* Initialize and allocate storage for GENMMD. */
	delta = 0; /* DELTA is a parameter to allow the choice of nodes
		      whose degree <= min-degree + DELTA. */
	maxint = 9223372036854775807; /* 2**63 - 1 */
	invp = (long long *) SUPERLU_MALLOC((n+delta)*sizeof(long long));
	if ( !invp ) ABORT("SUPERLU_MALLOC fails for invp.");
	dhead = (long long *) SUPERLU_MALLOC((n+delta)*sizeof(long long));
	if ( !dhead ) ABORT("SUPERLU_MALLOC fails for dhead.");
	qsize = (long long *) SUPERLU_MALLOC((n+delta)*sizeof(long long));
	if ( !qsize ) ABORT("SUPERLU_MALLOC fails for qsize.");
	llist = (long long *) SUPERLU_MALLOC(n*sizeof(long long));
	if ( !llist ) ABORT("SUPERLU_MALLOC fails for llist.");
	marker = (long long *) SUPERLU_MALLOC(n*sizeof(long long));
	if ( !marker ) ABORT("SUPERLU_MALLOC fails for marker.");

	/* Transform adjacency list into 1-based indexing required by GENMMD.*/
	for (i = 0; i <= n; ++i) ++b_colptr[i];
	for (i = 0; i < bnz; ++i) ++b_rowind[i];
	
	genmmd_(&n, b_colptr, b_rowind, perm_c, invp, &delta, dhead, 
		qsize, llist, marker, &maxint, &nofsub);

	/* Transform perm_c into 0-based indexing. */
	for (i = 0; i < n; ++i) --perm_c[i];

	SUPERLU_FREE(invp);
	SUPERLU_FREE(dhead);
	SUPERLU_FREE(qsize);
	SUPERLU_FREE(llist);
	SUPERLU_FREE(marker);
	SUPERLU_FREE(b_rowind);

	t = SuperLU_timer_() - t;
	/*  printf("call GENMMD time = %8.3f\n", t);*/

    } else { /* Empty adjacency structure */
	for (i = 0; i < n; ++i) perm_c[i] = i;
    }

    SUPERLU_FREE(b_colptr);
}
