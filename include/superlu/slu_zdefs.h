/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_zdefs.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * 
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
#ifndef __SUPERLU_zSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_zSP_DEFS

/*
 * File name:		zsp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

/* Define my integer type int_t */
typedef long long int_t; /* default */

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"
#include "slu_dcomplex.h"


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
zgssv(superlu_options_t *, SuperMatrix *, long long *, long long *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, long long *);
extern void
zgssvx(superlu_options_t *, SuperMatrix *, long long *, long long *, long long *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, long long, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, long long *);
    /* ILU */
extern void
zgsisv(superlu_options_t *, SuperMatrix *, long long *, long long *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, long long *);
extern void
zgsisx(superlu_options_t *, SuperMatrix *, long long *, long long *, long long *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, long long, SuperMatrix *, SuperMatrix *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, long long *);


/*! \brief Supernodal LU factor related */
extern void
zCreate_CompCol_Matrix(SuperMatrix *, long long, long long, long long, doublecomplex *,
		       long long *, long long *, Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_CompRow_Matrix(SuperMatrix *, long long, long long, long long, doublecomplex *,
		       long long *, long long *, Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
zCreate_Dense_Matrix(SuperMatrix *, long long, long long, doublecomplex *, long long,
		     Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Matrix(SuperMatrix *, long long, long long, long long, doublecomplex *, 
		         long long *, long long *, long long *, long long *, long long *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
zCopy_Dense_Matrix(long long, long long, doublecomplex *, long long, doublecomplex *, long long);

extern void    countnz (const long long, long long *, long long *, long long *, GlobalLU_t *);
extern void    ilu_countnz (const long long, long long *, long long *, GlobalLU_t *);
extern void    fixupL (const long long, const long long *, GlobalLU_t *);

extern void    zallocateA (long long, long long, doublecomplex **, long long **, long long **);
extern void    zgstrf (superlu_options_t*, SuperMatrix*,
                       long long, long long, long long*, void *, long long, long long *, long long *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, long long *);
extern long long     zsnode_dfs (const long long, const long long, const long long *, const long long *,
			     const long long *, long long *, long long *, GlobalLU_t *);
extern long long     zsnode_bmod (const long long, const long long, const long long, doublecomplex *,
                              doublecomplex *, GlobalLU_t *, SuperLUStat_t*);
extern void    zpanel_dfs (const long long, const long long, const long long, SuperMatrix *,
			   long long *, long long *, doublecomplex *, long long *, long long *, long long *,
			   long long *, long long *, long long *, long long *, GlobalLU_t *);
extern void    zpanel_bmod (const long long, const long long, const long long, const long long,
                           doublecomplex *, doublecomplex *, long long *, long long *,
			   GlobalLU_t *, SuperLUStat_t*);
extern long long     zcolumn_dfs (const long long, const long long, long long *, long long *, long long *, long long *,
			   long long *, long long *, long long *, long long *, long long *, GlobalLU_t *);
extern long long     zcolumn_bmod (const long long, const long long, doublecomplex *,
			   doublecomplex *, long long *, long long *, long long,
                           GlobalLU_t *, SuperLUStat_t*);
extern long long     zcopy_to_ucol (long long, long long, long long *, long long *, long long *,
                              doublecomplex *, GlobalLU_t *);         
extern long long     zpivotL (const long long, const double, long long *, long long *, 
                         long long *, long long *, long long *, GlobalLU_t *, SuperLUStat_t*);
extern void    zpruneL (const long long, const long long *, const long long, const long long,
			  const long long *, const long long *, long long *, GlobalLU_t *);
extern void    zreadmt (long long *, long long *, long long *, doublecomplex **, long long **, long long **);
extern void    zGenXtrue (long long, long long, doublecomplex *, long long);
extern void    zFillRHS (trans_t, long long, doublecomplex *, long long, SuperMatrix *,
			  SuperMatrix *);
extern void    zgstrs (trans_t, SuperMatrix *, SuperMatrix *, long long *, long long *,
                        SuperMatrix *, SuperLUStat_t*, long long *);
/* ILU */
extern void    zgsitrf (superlu_options_t*, SuperMatrix*, long long, long long, long long*,
		        void *, long long, long long *, long long *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, long long *);
extern long long     zldperm(long long, long long, long long, long long [], long long [], doublecomplex [],
                        long long [],	double [], double []);
extern long long     ilu_zsnode_dfs (const long long, const long long, const long long *, const long long *,
			       const long long *, long long *, GlobalLU_t *);
extern void    ilu_zpanel_dfs (const long long, const long long, const long long, SuperMatrix *,
			       long long *, long long *, doublecomplex *, double *, long long *, long long *,
			       long long *, long long *, long long *, long long *, GlobalLU_t *);
extern long long     ilu_zcolumn_dfs (const long long, const long long, long long *, long long *, long long *,
				long long *, long long *, long long *, long long *, long long *,
				GlobalLU_t *);
extern long long     ilu_zcopy_to_ucol (long long, long long, long long *, long long *, long long *,
                                  doublecomplex *, long long, milu_t, double, long long,
                                  doublecomplex *, long long *, GlobalLU_t *, double *);
extern long long     ilu_zpivotL (const long long, const double, long long *, long long *, long long, long long *,
			    long long *, long long *, long long *, double, milu_t,
                            doublecomplex, GlobalLU_t *, SuperLUStat_t*);
extern long long     ilu_zdrop_row (superlu_options_t *, long long, long long, double,
                              long long, long long *, double *, GlobalLU_t *, 
                              double *, double *, long long);


/*! \brief Driver related */

extern void    zgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, long long *);
extern void    zlaqgs (SuperMatrix *, double *, double *, double,
                        double, double, char *);
extern void    zgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, long long *);
extern double   zPivotGrowth(long long, SuperMatrix *, long long *, 
                            SuperMatrix *, SuperMatrix *);
extern void    zgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, long long *, long long *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, long long *);

extern long long     sp_ztrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, doublecomplex *, SuperLUStat_t*, long long *);
extern long long     sp_zgemv (char *, doublecomplex, SuperMatrix *, doublecomplex *,
			long long, doublecomplex, doublecomplex *, long long);

extern long long     sp_zgemm (char *, char *, long long, long long, long long, doublecomplex,
			SuperMatrix *, doublecomplex *, long long, doublecomplex, 
			doublecomplex *, long long);
extern         double dmach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern long long     zLUMemInit (fact_t, void *, long long, long long, long long, long long, long long,
                            double, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, long long **, doublecomplex **);
extern void    zSetRWork (long long, long long, doublecomplex *, doublecomplex **, doublecomplex **);
extern void    zLUWorkFree (long long *, doublecomplex *, GlobalLU_t *);
extern long long     zLUMemXpand (long long, long long, MemType, long long *, GlobalLU_t *);

extern doublecomplex  *doublecomplexMalloc(long long);
extern doublecomplex  *doublecomplexCalloc(long long);
extern double  *doubleMalloc(long long);
extern double  *doubleCalloc(long long);
extern long long     zmemory_usage(const long long, const long long, const long long, const long long);
extern long long     zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern long long     ilu_zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    zreadhb(FILE *, long long *, long long *, long long *, doublecomplex **, long long **, long long **);
extern void    zreadrb(long long *, long long *, long long *, doublecomplex **, long long **, long long **);
extern void    zreadtriple(long long *, long long *, long long *, doublecomplex **, long long **, long long **);
extern void    zreadMM(FILE *, long long *, long long *, long long *, doublecomplex **, long long **, long long **);
extern void    zCompRow_to_CompCol(long long, long long, long long, doublecomplex*, long long*, long long*,
		                   doublecomplex **, long long **, long long **);
extern void    zfill (doublecomplex *, long long, doublecomplex);
extern void    zinf_norm_error (long long, SuperMatrix *, doublecomplex *);
extern double  dqselect(long long, double *, long long);


/*! \brief Routines for debugging */
extern void    zPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    zPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    zPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    zprint_lu_col(char *, long long, long long, long long *, GlobalLU_t *);
extern long long     print_double_vec(char *, long long, double *);
extern void    zcheck_tempv(long long, doublecomplex *);

/*! \brief BLAS */

extern long long zgemm_(const char*, const char*, const long long*, const long long*, const long long*,
                  const doublecomplex*, const doublecomplex*, const long long*, const doublecomplex*,
		  const long long*, const doublecomplex*, doublecomplex*, const long long*);
extern long long ztrsv_(char*, char*, char*, long long*, doublecomplex*, long long*,
                  doublecomplex*, long long*);
extern long long ztrsm_(char*, char*, char*, char*, long long*, long long*,
                  doublecomplex*, doublecomplex*, long long*, doublecomplex*, long long*);
extern long long zgemv_(char *, long long *, long long *, doublecomplex *, doublecomplex *a, long long *,
                  doublecomplex *, long long *, doublecomplex *, doublecomplex *, long long *);

extern void zusolve(long long, long long, doublecomplex*, doublecomplex*);
extern void zlsolve(long long, long long, doublecomplex*, doublecomplex*);
extern void zmatvec(long long, long long, long long, doublecomplex*, doublecomplex*, doublecomplex*);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_zSP_DEFS */

