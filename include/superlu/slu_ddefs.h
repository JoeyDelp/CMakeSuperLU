/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_ddefs.h
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
#ifndef __SUPERLU_dSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_dSP_DEFS

/*
 * File name:		dsp_defs.h
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


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
dgssv(superlu_options_t *, SuperMatrix *, long long *, long long *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, long long *);
extern void
dgssvx(superlu_options_t *, SuperMatrix *, long long *, long long *, long long *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, long long, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, long long *);
    /* ILU */
extern void
dgsisv(superlu_options_t *, SuperMatrix *, long long *, long long *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, long long *);
extern void
dgsisx(superlu_options_t *, SuperMatrix *, long long *, long long *, long long *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, long long, SuperMatrix *, SuperMatrix *, double *, double *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, long long *);


/*! \brief Supernodal LU factor related */
extern void
dCreate_CompCol_Matrix(SuperMatrix *, long long, long long, long long, double *,
		       long long *, long long *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompRow_Matrix(SuperMatrix *, long long, long long, long long, double *,
		       long long *, long long *, Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
dCreate_Dense_Matrix(SuperMatrix *, long long, long long, double *, long long,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix(SuperMatrix *, long long, long long, long long, double *, 
		         long long *, long long *, long long *, long long *, long long *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_Dense_Matrix(long long, long long, double *, long long, double *, long long);

extern void    countnz (const long long, long long *, long long *, long long *, GlobalLU_t *);
extern void    ilu_countnz (const long long, long long *, long long *, GlobalLU_t *);
extern void    fixupL (const long long, const long long *, GlobalLU_t *);

extern void    dallocateA (long long, long long, double **, long long **, long long **);
extern void    dgstrf (superlu_options_t*, SuperMatrix*,
                       long long, long long, long long*, void *, long long, long long *, long long *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, long long *);
extern long long     dsnode_dfs (const long long, const long long, const long long *, const long long *,
			     const long long *, long long *, long long *, GlobalLU_t *);
extern long long     dsnode_bmod (const long long, const long long, const long long, double *,
                              double *, GlobalLU_t *, SuperLUStat_t*);
extern void    dpanel_dfs (const long long, const long long, const long long, SuperMatrix *,
			   long long *, long long *, double *, long long *, long long *, long long *,
			   long long *, long long *, long long *, long long *, GlobalLU_t *);
extern void    dpanel_bmod (const long long, const long long, const long long, const long long,
                           double *, double *, long long *, long long *,
			   GlobalLU_t *, SuperLUStat_t*);
extern long long     dcolumn_dfs (const long long, const long long, long long *, long long *, long long *, long long *,
			   long long *, long long *, long long *, long long *, long long *, GlobalLU_t *);
extern long long     dcolumn_bmod (const long long, const long long, double *,
			   double *, long long *, long long *, long long,
                           GlobalLU_t *, SuperLUStat_t*);
extern long long     dcopy_to_ucol (long long, long long, long long *, long long *, long long *,
                              double *, GlobalLU_t *);         
extern long long     dpivotL (const long long, const double, long long *, long long *, 
                         long long *, long long *, long long *, GlobalLU_t *, SuperLUStat_t*);
extern void    dpruneL (const long long, const long long *, const long long, const long long,
			  const long long *, const long long *, long long *, GlobalLU_t *);
extern void    dreadmt (long long *, long long *, long long *, double **, long long **, long long **);
extern void    dGenXtrue (long long, long long, double *, long long);
extern void    dFillRHS (trans_t, long long, double *, long long, SuperMatrix *,
			  SuperMatrix *);
extern void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, long long *, long long *,
                        SuperMatrix *, SuperLUStat_t*, long long *);
/* ILU */
extern void    dgsitrf (superlu_options_t*, SuperMatrix*, long long, long long, long long*,
		        void *, long long, long long *, long long *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, long long *);
extern long long     dldperm(long long, long long, long long, long long [], long long [], double [],
                        long long [],	double [], double []);
extern long long     ilu_dsnode_dfs (const long long, const long long, const long long *, const long long *,
			       const long long *, long long *, GlobalLU_t *);
extern void    ilu_dpanel_dfs (const long long, const long long, const long long, SuperMatrix *,
			       long long *, long long *, double *, double *, long long *, long long *,
			       long long *, long long *, long long *, long long *, GlobalLU_t *);
extern long long     ilu_dcolumn_dfs (const long long, const long long, long long *, long long *, long long *,
				long long *, long long *, long long *, long long *, long long *,
				GlobalLU_t *);
extern long long     ilu_dcopy_to_ucol (long long, long long, long long *, long long *, long long *,
                                  double *, long long, milu_t, double, long long,
                                  double *, long long *, GlobalLU_t *, double *);
extern long long     ilu_dpivotL (const long long, const double, long long *, long long *, long long, long long *,
			    long long *, long long *, long long *, double, milu_t,
                            double, GlobalLU_t *, SuperLUStat_t*);
extern long long     ilu_ddrop_row (superlu_options_t *, long long, long long, double,
                              long long, long long *, double *, GlobalLU_t *, 
                              double *, double *, long long);


/*! \brief Driver related */

extern void    dgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, long long *);
extern void    dlaqgs (SuperMatrix *, double *, double *, double,
                        double, double, char *);
extern void    dgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, long long *);
extern double   dPivotGrowth(long long, SuperMatrix *, long long *, 
                            SuperMatrix *, SuperMatrix *);
extern void    dgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, long long *, long long *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, long long *);

extern long long     sp_dtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, double *, SuperLUStat_t*, long long *);
extern long long     sp_dgemv (char *, double, SuperMatrix *, double *,
			long long, double, double *, long long);

extern long long     sp_dgemm (char *, char *, long long, long long, long long, double,
			SuperMatrix *, double *, long long, double, 
			double *, long long);
extern         double dmach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern long long     dLUMemInit (fact_t, void *, long long, long long, long long, long long, long long,
                            double, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, long long **, double **);
extern void    dSetRWork (long long, long long, double *, double **, double **);
extern void    dLUWorkFree (long long *, double *, GlobalLU_t *);
extern long long     dLUMemXpand (long long, long long, MemType, long long *, GlobalLU_t *);

extern double  *doubleMalloc(long long);
extern double  *doubleCalloc(long long);
extern long long     dmemory_usage(const long long, const long long, const long long, const long long);
extern long long     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern long long     ilu_dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    dreadhb(FILE *, long long *, long long *, long long *, double **, long long **, long long **);
extern void    dreadrb(long long *, long long *, long long *, double **, long long **, long long **);
extern void    dreadtriple(long long *, long long *, long long *, double **, long long **, long long **);
extern void    dreadMM(FILE *, long long *, long long *, long long *, double **, long long **, long long **);
extern void    dCompRow_to_CompCol(long long, long long, long long, double*, long long*, long long*,
		                   double **, long long **, long long **);
extern void    dfill (double *, long long, double);
extern void    dinf_norm_error (long long, SuperMatrix *, double *);
extern double  dqselect(long long, double *, long long);


/*! \brief Routines for debugging */
extern void    dPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    dPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    dPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    dprint_lu_col(char *, long long, long long, long long *, GlobalLU_t *);
extern long long     print_double_vec(char *, long long, double *);
extern void    dcheck_tempv(long long, double *);

/*! \brief BLAS */

extern long long dgemm_(const char*, const char*, const long long*, const long long*, const long long*,
                  const double*, const double*, const long long*, const double*,
		  const long long*, const double*, double*, const long long*);
extern long long dtrsv_(char*, char*, char*, long long*, double*, long long*,
                  double*, long long*);
extern long long dtrsm_(char*, char*, char*, char*, long long*, long long*,
                  double*, double*, long long*, double*, long long*);
extern long long dgemv_(char *, long long *, long long *, double *, double *a, long long *,
                  double *, long long *, double *, double *, long long *);

extern void dusolve(long long, long long, double*, double*);
extern void dlsolve(long long, long long, double*, double*);
extern void dmatvec(long long, long long, long long, double*, double*, double*);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_dSP_DEFS */

