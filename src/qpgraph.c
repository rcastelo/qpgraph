/*
  qpgraph package - this C code implements functions to learn qp-graphs from
                    data and utilities for GGM model inference and simulation
 
  Copyright (C) 2011 R. Castelo and A. Roverato
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, you can obtain one via WWW at
  http://www.gnu.org/copyleft/gpl.html, or by writing to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/



#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include "cliquer.h"

/* constants */

#define E2I(v,w) (v > w ? ((int) (((double) ( v * (v - 1))) / 2.0)) + w : ((int) (((double) (w * (w - 1))) / 2.0)) + v) 
#define UTE2I(v,w) (v > w ? ((int) (((double) ( v * (v - 1))) / 2.0)) + w + v : ((int) (((double) (w * (w - 1))) / 2.0)) + v + w) 

#define USE_COMPLETE_OBS 1
#define USE_EM           2

#define RETURN_TYPE_PVALUE 1
#define RETURN_TYPE_STATN  2
#define RETURN_TYPE_ALL    3

/* datatype definitions */

typedef struct tag_clique_t {
  union {
    set_t vts; /* vertex set - vertices stored using cliquer sets */
    int * vta; /* vertex array - vertices stored as an array of integers */
  } u;
  int n;

  struct tag_clique_t* next;
} clique_t;
                                                                                                
typedef struct {
  clique_t* first;
  clique_t* last;
  int       n;
} clique_set_t;

typedef struct {
  double* Es_com;
  double* Ess_com;
  int*    n_com;
} com_stats_t;

typedef struct {
  double* ssd;
  double* K;
  double* h;
  double* m;
  double* bar_y;
} suf_stats_t;

/* function prototypes */

static SEXP installAttrib(SEXP, SEXP, SEXP);

static SEXP installAttrib(SEXP vec, SEXP name, SEXP val)
{
    SEXP s, t;

    if(TYPEOF(vec) == CHARSXP)
	error("cannot set attribute on a CHARSXP");
    PROTECT(vec);
    PROTECT(name);
    PROTECT(val);
    for (s = ATTRIB(vec); s != R_NilValue; s = CDR(s)) {
	if (TAG(s) == name) {
	    SETCAR(s, val);
	    UNPROTECT(3);
	    return val;
	}
    }
    s = allocList(1);
    SETCAR(s, val);
    SET_TAG(s, name);
    if (ATTRIB(vec) == R_NilValue)
	SET_ATTRIB(vec, s);
    else {
	t = nthcdr(ATTRIB(vec), length(ATTRIB(vec)) - 1);
	SETCDR(t, s);
    }
    UNPROTECT(3);
    return val;
}

/* from Mutils.h */
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, int length)
{
  SEXP val = allocVector(type, length);

  SET_SLOT(obj, nm, val);
  return val;
}


extern void R_FlushConsole(void);
extern void R_CheckUserInterrupt(void);
#ifdef Win32
extern void R_ProcessEvents(void);
#endif
#ifdef HAVE_AQUA
extern void R_ProcessEvents(void);
#endif

static SEXP
qp_fast_nrr(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP qR, SEXP restrictQR,
            SEXP fixQR, SEXP nTestsR, SEXP alphaR, SEXP pairup_i_nointR,
            SEXP pairup_j_nointR, SEXP pairup_ij_intR, SEXP exactTest,
            SEXP verboseR, SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP env);

static SEXP
qp_fast_nrr_identicalQs(SEXP XR, SEXP qR, SEXP restrictQR, SEXP fixQR, SEXP nTests,
                        SEXP alpha, SEXP pairup_i_noint, SEXP pairup_j_noint,
                        SEXP pairup_ij_int, SEXP verbose, SEXP startTimeR,
                        SEXP nAdj2estimateTimeR, SEXP env);

static SEXP
qp_fast_nrr_par(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP qR,
                SEXP restrictQR, SEXP fixQR, SEXP nTestsR, SEXP alphaR,
                SEXP pairup_i_nointR, SEXP pairup_j_nointR,
                SEXP pairup_ij_intR, SEXP exactTest, SEXP verboseR,
                SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP myRankR,
                SEXP clSzeR, SEXP masterNode, SEXP env);

static SEXP
qp_fast_nrr_identicalQs_par(SEXP XR, SEXP qR, SEXP restrictQR, SEXP fixQR,
                            SEXP nTestsR, SEXP alphaR, SEXP pairup_i_nointR,
                            SEXP pairup_j_nointR, SEXP pairup_ij_intR, SEXP verboseR,
                            SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP myRankR, SEXP clSzeR,
                            SEXP masterNode, SEXP env);

static SEXP
qp_fast_edge_nrr(SEXP XR, SEXP S, SEXP n_varR, SEXP nR, SEXP iR, SEXP jR, SEXP qR,
                 SEXP restrictQR, SEXP fixQR, SEXP nTestsR, SEXP alphaR);

static SEXP
qp_fast_edge_nrr_hmgm(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP ssdR,
                      SEXP mapX2ssdR, SEXP iR, SEXP jR, SEXP qR, SEXP restrictQR,
                      SEXP fixQR, SEXP nTestsR, SEXP alphaR, SEXP exactTest);

static SEXP
qp_fast_edge_nrr_hmgm_sml(SEXP XR, SEXP cumsum_sByChrR, SEXP sR, SEXP gLevelsR,
                          SEXP XEPR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP ssdR,
                          SEXP mapX2ssdR, SEXP iR, SEXP jR, SEXP qR, SEXP restrictQR,
                          SEXP fixQR, SEXP nTestsR, SEXP alphaR, SEXP exactTest);

static double
qp_edge_nrr(double* X, double* S, int p, int n, int i, int j, int q,
            int* restrictQ, int n_rQ, int* fixQ, int n_fQ, int nTests, double alpha);

static double
qp_edge_nrr_identicalQs(double* S, int n_var, int* Qs, double* Qinvs, int N, int i,
                        int j, int q, int nTests, double alpha);

static double
qp_edge_nrr_hmgm(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y,
                 int n_Y, double* ucond_ssd, int* mapX2ucond_ssd, int i, int j,
                 int q, int* restrictQ, int n_rQ, int* fixQ, int n_fQ, int nTests,
                 double alpha, int exactTest);

static double
qp_edge_nrr_hmgm_sml(SEXP X, int* cumsum_sByChr, int s, int gLevels, double* XEP1q,
                     int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
                     double* ucond_ssd, int* mapX2ucond_ssd, int i, int j, int q,
                     int* restrictQ, int n_rQ, int* fixQ, int n_fQ, int nTests,
                     double alpha, int exactTest);

static void
sampleQs(int T, int q, int v_i, int v_j, int p, int* restrictQ, int* fixQ,
         int n_fQ, int* y);

static SEXP
qp_fast_ci_test_std(SEXP SR, SEXP n_varR, SEXP NR, SEXP iR, SEXP jR, SEXP C);
static double
qp_ci_test_std(double* S, int n_var, int N, int i, int j, int* C, int q, double*);

static SEXP
qp_fast_ci_test_opt(SEXP SR, SEXP n_varR, SEXP NR, SEXP iR, SEXP jR, SEXP C);
static double
qp_ci_test_opt(double* S, int n_var, int N, int i, int j, int* C, int q, double*,
               double*);

static SEXP
qp_fast_ci_test_hmgm(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP ssdR,
                     SEXP mapX2ssdR, SEXP iR, SEXP jR, SEXP QR, SEXP exactTest,
                     SEXP use, SEXP tol);

static double
lr_complete_obs(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y,
                int n_Y, double* ucond_ssd, int* mapX2ucond_ssd, int total_n_Y,
                int i, int j, int* Q, int q, int* n_co);

com_stats_t
new_com_stats(int n_joint_levels, int n_Y);

void
free_com_stats(com_stats_t cs);

com_stats_t
stat_com(double* X, int p, int n, int* missing_mask, int n_mis, int* Is, int n_Is,
         int *Y, int n_Y, int* n_levels, int n_joint_levels);

static double
lr_em(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
      int i, int j, int* Q, int q, double tol);

static double
qp_ci_test_hmgm(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y,
                int n_Y, double* ucond_ssd, int* mapX2ucond_ssd, int i, int j, int* Q,
                int q, int use, double tol, double* df, double* a, double* b, int* n_co);

static double
qp_ci_test_hmgm_sml(SEXP Xsml, int* cumsum_sByChr, int s, int gLevels, double* XEP1q,
                    int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
                    double* ucond_ssd, int* mapX2ucond_ssd, int i, int j, int* Q, int q,
                    int use, double tol, double* df, double* a, double* b, int* n_co);

static SEXP
qp_fast_all_ci_tests(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP QR,
                     SEXP pairup_i_nointR, SEXP pairup_j_nointR, SEXP pairup_ij_intR,
                     SEXP exactTest, SEXP use, SEXP tol, SEXP return_type, SEXP verboseR,
                     SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP env);

static SEXP
qp_fast_all_ci_tests_par(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP QR,
                         SEXP pairup_i_nointR, SEXP pairup_j_nointR,
                         SEXP pairup_ij_intR, SEXP exactTest, SEXP use, SEXP tol,
                         SEXP return_type, SEXP verboseR, SEXP startTimeR,
                         SEXP nAdj2estimateTimeR, SEXP myRankR, SEXP clSzeR,
                         SEXP masterNode, SEXP env);

boolean
cliquer_cb_add_clique_to_list(set_t clique, graph_t* g, clique_options* opts);

void
add_clique_vts(clique_set_t* cset, set_t clique);

void
add_clique_vta(clique_set_t* cset, int* clique, int n);

void
destroy_cliques_vts(clique_set_t* cset);

void
destroy_cliques_vta(clique_set_t* cset);

void
init_cliques_list(clique_set_t* cset);

static SEXP
qp_fast_cliquer_get_cliques(SEXP I, SEXP clqspervtx, SEXP verbose);

Rboolean
is_maximal_clique(int* I, int n, int* clq, int cs, set_t noclq);

static SEXP
qp_fast_update_cliques_removing(SEXP I, SEXP clqlstR, SEXP vR, SEXP wR, SEXP verbose);

static SEXP
qp_clique_number_lb(SEXP I, SEXP return_vertices, SEXP approx_iter, SEXP verbose);

static SEXP
qp_clique_number_os(SEXP I, SEXP return_vertices, SEXP verbose);

int
clique_number_os(int n, int* I, int verbose);

static SEXP
qp_fast_pac_se(SEXP Shat, SEXP I);

static SEXP
qp_fast_ipf(SEXP vv, SEXP cliq, SEXP tol, SEXP verbose);

static SEXP
qp_fast_htf(SEXP S, SEXP A, SEXP tol, SEXP verbose);

static void
fast_ipf_step(int n, double* Vf, double* Vn, int* a, int csze);

static void
cov2cor(double* R, double* V, int n);

static void
matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
 
static void
matinv(double* inv, double* M, int n, int p);

static double
symmatlogdet(double* M, int n, int* sign);

static void
matsumf(double* R, int nrow, int ncol, double* M, double* N, double factor);

static void
matscalarprod(double* R, int nrow, int ncol, double* M, double* N);

static void
mattran(double* T, double* M, int nrow, int ncol);

static void
matsubm(double* subM, double* M, int n, int* subrows,
        int nsubrows, int* subcols, int nsubcols);

static void
symmatsubm(double* subM, double* M, int n, int* subrows, int nsubrows);

static double
matmxab(double* M, int nrow, int ncol);

static void
matrepm(double* M, int n, int* subrows, int nsubrows,
        int* subcols, int nsubcols, double* N);

static void
setdiff(int n, int m, int* a, int* b);

void
i2e(int i, int* e_i, int* e_j);

int
e2i(int e_i, int e_j, int* i);

int
missing_obs(double* X, int p, int n, int* Y, int n_Y, int* idx_obs, int n_idx_obs);

int
find_missing_obs(double* X, int p, int n, int* Y, int n_Y, int* idx_obs,
                 int n_idx_obs, int* missing_mask);

void
calculate_means(double* X, int p, int n, int* Y, int n_Y, int* idx_obs,
                int n_idx_obs, int* missing_mask, int n_mis, double* meanv);

int
ssd(double* X, int p, int n, int* Y, int n_Y, int* idx_obs, int n_idx_obs,
    int corrected, int* missing_mask, double* ssd_mat);

static SEXP
qp_fast_cov_upper_triangular(SEXP XR, SEXP corrected);

static SEXP
qp_fast_rnd_graph(SEXP pR, SEXP dR, SEXP excludeR, SEXP verbose);

void
calculate_xtab(double* X, int p, int n, int* I, int n_I, int* n_levels, int* xtab);

int
ssd_A(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
      int* excobs_mask, int* missing_mask, double* ssd_A);

/* Global variables */

SEXP Matrix_DimNamesSym,
     Matrix_DimSym,
     Matrix_uploSym,
     Matrix_xSym,
     SsdMatrix_ssdSym,
     SsdMatrix_nSym,

     qpgraph_NS; /* the qpgraph namespace ('environment') */

int* global_xtab; /* for cross-classifying joint levels of discrete variables */

/* for comparing cross-classified joint levels */
int
indirect_int_cmp(const void *a, const void *b) {
  return ( global_xtab[*(int*)a] - global_xtab[*(int*)b] );
}

/* R-function register */

static R_CallMethodDef
callMethods[] = {
  {"qp_fast_nrr", (DL_FUNC) &qp_fast_nrr, 17},
  {"qp_fast_nrr_identicalQs", (DL_FUNC) &qp_fast_nrr_identicalQs, 13},
  {"qp_fast_nrr_par", (DL_FUNC) &qp_fast_nrr_par, 20},
  {"qp_fast_nrr_identicalQs_par", (DL_FUNC) &qp_fast_nrr_identicalQs_par, 16},
  {"qp_fast_edge_nrr", (DL_FUNC) &qp_fast_edge_nrr, 11},
  {"qp_fast_edge_nrr_hmgm", (DL_FUNC) &qp_fast_edge_nrr_hmgm, 14},
  {"qp_fast_edge_nrr_hmgm_sml", (DL_FUNC) &qp_fast_edge_nrr_hmgm_sml, 18},
  {"qp_fast_ci_test_std", (DL_FUNC) &qp_fast_ci_test_std, 6},
  {"qp_fast_ci_test_opt", (DL_FUNC) &qp_fast_ci_test_opt, 6},
  {"qp_fast_ci_test_hmgm", (DL_FUNC) &qp_fast_ci_test_hmgm, 12},
  {"qp_fast_all_ci_tests", (DL_FUNC) &qp_fast_all_ci_tests, 16},
  {"qp_fast_all_ci_tests_par", (DL_FUNC) &qp_fast_all_ci_tests_par, 19},
  {"qp_fast_cliquer_get_cliques", (DL_FUNC) &qp_fast_cliquer_get_cliques, 3},
  {"qp_fast_update_cliques_removing", (DL_FUNC) &qp_fast_update_cliques_removing, 5},
  {"qp_clique_number_lb", (DL_FUNC) &qp_clique_number_lb, 4},
  {"qp_clique_number_os", (DL_FUNC) &qp_clique_number_os, 3},
  {"qp_fast_pac_se", (DL_FUNC) &qp_fast_pac_se, 2},
  {"qp_fast_ipf", (DL_FUNC) &qp_fast_ipf, 4},
  {"qp_fast_htf", (DL_FUNC) &qp_fast_htf, 4},
  {"qp_fast_cov_upper_triangular", (DL_FUNC) &qp_fast_cov_upper_triangular, 2},
  {"qp_fast_rnd_graph", (DL_FUNC) &qp_fast_rnd_graph, 4},
  {NULL}
};

void
R_init_qpgraph(DllInfo* info) {

  R_registerRoutines(info,NULL,callMethods,NULL,0);

  /* from the Matrix package init.c */
  Matrix_DimNamesSym = install("Dimnames");
  Matrix_DimSym = install("Dim");
  Matrix_uploSym = install("uplo");
  Matrix_xSym = install("x");
  SsdMatrix_ssdSym = install("ssd");
  SsdMatrix_nSym = install("n");

  qpgraph_NS = R_FindNamespace(mkString("qpgraph"));
  if (qpgraph_NS == R_UnboundValue)
    error("missing 'qpgraph' namespace: should never happen");

  GetRNGstate(); /* initialize the R-builtin RNG */

}



/*
  FUNCTION: qp_fast_nrr
  PURPOSE: compute for each pair of vertices indexed by the columns of the
           matrix X the non-rejection rate. Vertex pairs may be restricted
           by using the pairup_* arguments
  RETURNS: matrix of non-rejection rate values in terms of number of non-rejected
           (accepted) tests for each pair of vertices
*/

static SEXP
qp_fast_nrr(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP qR, SEXP restrictQR,
            SEXP fixQR, SEXP nTestsR, SEXP alphaR, SEXP pairup_i_nointR,
            SEXP pairup_j_nointR, SEXP pairup_ij_intR, SEXP exactTest,
            SEXP verboseR, SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP env) {
  int     N;
  int     n_var;
  int     q;
  int     nTests;
  double  alpha;
  double* S = NULL;
  double* ssdMat = NULL;
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int*    restrictQ = NULL;
  int     n_rQ = 0;
  int*    fixQ = NULL;
  int     n_fQ = 0;
  int     isMatrix_restrictQ = FALSE;
  int     work_with_margin = FALSE;
  int     l_ini = length(pairup_i_nointR);
  int     l_jni = length(pairup_j_nointR);
  int     l_int = length(pairup_ij_intR);
  int*    I = NULL;
  int*    Y = NULL;
  int*    mapX2ssd = NULL;
  int*    pairup_i_noint = INTEGER(pairup_i_nointR);
  int*    pairup_j_noint = INTEGER(pairup_j_nointR);
  int*    pairup_ij_int = INTEGER(pairup_ij_intR);
  int*    pairup_ij_noint = NULL;
  int     i,j,k;
  int     n_adj, pct, ppct;
  SEXP    nrrR;
  double* nrr;
  int     verbose;
  double  startTime, elapsedTime;
  int     nAdjEtime;
  SEXP    pb=NULL;

  N         = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  n_var     = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  q         = INTEGER(qR)[0];
  nTests    = INTEGER(nTestsR)[0];
  alpha     = REAL(alphaR)[0];
  verbose   = INTEGER(verboseR)[0];
  startTime = REAL(startTimeR)[0];
  nAdjEtime = INTEGER(nAdj2estimateTimeR)[0];

  if (q > n_var-2)
    error("q=%d > p-2=%d",q,n_var-2);

  if (q < 0)
    error("q=%d < 0",q);

  if (q > N-3)
    error("q=%d > n-3=%d", q, N-3);

  if (nTests < 1)
    error("nTests=%d < 1", nTests);

  if (alpha < 0.0 || alpha > 1.0)
    error("significance level alpha is %.2f and it should lie in the interval [0, 1]\n", alpha);


  if (n_I == 0) {
    if (!missing_obs(REAL(XR), n_var, N, NULL, n_var, NULL, N)) {
      S = ssdMat = Calloc((n_var*(n_var+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, N, NULL, n_var, NULL, N, TRUE, NULL, S);
    } else
      work_with_margin = TRUE;
  } else {
    I = Calloc(n_I, int);
    for (i=0; i < n_I; i++)
      I[i] = INTEGER(IR)[i]-1;

    Y = Calloc(n_Y, int);
    for (i=0; i < n_Y; i++)
      Y[i] = INTEGER(YR)[i]-1;

    if (!missing_obs(REAL(XR), n_var, N, NULL, n_var, NULL, N)) {
      mapX2ssd = Calloc(n_var, int);
      for (i=0; i < n_var; i++) {
        j = 0;
        while (j < n_Y && i != Y[j])
          j++;

        mapX2ssd[i] = j;
      }

      S = ssdMat = Calloc((n_Y*(n_Y+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, N, Y, n_Y, NULL, N, FALSE, NULL, ssdMat);
    } else
      work_with_margin = TRUE;
  }

  if (restrictQR != R_NilValue) {
    if (isMatrix(restrictQR)) {
      if (n_I == 0)
        error("restrict.Q as a matrix can only be employed for restricting conditioning of discrete variables\n");
      isMatrix_restrictQ = TRUE;
      restrictQ = Calloc(n_var, int);
    } else {
      n_rQ = length(restrictQR);
      restrictQ = Calloc(n_rQ, int);
      for (i=0; i < n_rQ; i++)
        restrictQ[i] = INTEGER(restrictQR)[i] - 1;
    }
  }

  if (fixQR != R_NilValue) {
    n_fQ = length(fixQR);
    fixQ = Calloc(n_fQ, int);
    for (i=0; i < n_fQ; i++)
      fixQ[i] = INTEGER(fixQR)[i] - 1;
  }

  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  PROTECT(nrrR = allocVector(REALSXP, (n_var*(n_var+1))/2)); /* diagonal should be allocated */
  nrr = REAL(nrrR);

  for (i=0;i < (n_var*(n_var+1))/2;i++)
    nrr[i] = NA_REAL;

  n_adj = l_int * (l_jni + l_ini) + l_ini * l_jni + l_int * (l_int - 1) / 2;

  elapsedTime = 0.0;
  if (startTime > 0.0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = procTime[2] - startTime;
    startTime = procTime[2];
    UNPROTECT(2); /* call procTimeR */
  }

  ppct = -1;
  k = 0;
  if (verbose && startTime == 0) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("txtProgressBar")); t=CDR(t);
    SETCAR(t, ScalarInteger(3));
    SET_TAG(t, install("style"));
    PROTECT(pb = eval(s, R_GlobalEnv));
    UNPROTECT(1); /* t s */
  }

  /* intersection variables against ij-exclusive variables */
  for (i=0; i < l_int; i++) {
    int i2 = pairup_ij_int[i] - 1;

    for (j=0; j < l_ini + l_jni; j++) {
      int j2 = pairup_ij_noint[j] - 1;

      if (restrictQ != NULL) {
        if (isMatrix_restrictQ) {
          int l;

          n_rQ = 0;
          for (l=0; l < n_var; l++)
            if (INTEGER(restrictQR)[i2 + l*n_var] || INTEGER(restrictQR)[j2 + l*n_var])
              restrictQ[n_rQ++] = l;
        }
      }

      nrr[UTE2I(i2, j2)] = n_I == 0 ? qp_edge_nrr(REAL(XR), S, n_var, N, i2, j2, q,
                                                  restrictQ, n_rQ, fixQ, n_fQ, nTests, alpha) :
                                      qp_edge_nrr_hmgm(REAL(XR), n_var, N, I, n_I,
                                                       INTEGER(n_levelsR), Y, n_Y,
                                                       ssdMat, mapX2ssd, i2, j2, q,
                                                       restrictQ, n_rQ, fixQ, n_fQ,
                                                       nTests, alpha, INTEGER(exactTest)[0]);
      k++;
      if (startTime > 0 && k == nAdjEtime)
        break;
      pct = (int) ((k * 100) / n_adj);
      if (pct != ppct) {
        if (verbose && startTime == 0) {
          SEXP s, t;
          PROTECT(t = s = allocList(3));
          SET_TYPEOF(s, LANGSXP);
          SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
          SETCAR(t, pb);
          SET_TAG(t, install("pb")); t=CDR(t);
          SETCAR(t, ScalarReal(((double) pct) / 100.0));
          SET_TAG(t, install("value"));
          eval(s, R_GlobalEnv);
          UNPROTECT(1); /* t s */
        }
        R_CheckUserInterrupt();
#ifdef Win32
        R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
    if (startTime > 0 && k == nAdjEtime)
      break;
  }

  if (l_ini + l_jni > 0)
    Free(pairup_ij_noint);

  /* i-exclusive variables against j-exclusive variables */
  if (startTime == 0 || k < nAdjEtime) {
    for (i=0; i < l_ini; i++) {
      int i2 = pairup_i_noint[i] - 1;

      for (j=0; j < l_jni; j++) {
        int j2 = pairup_j_noint[j] - 1;

        if (restrictQ != NULL) {
          if (isMatrix_restrictQ) {
            int l;

            n_rQ = 0;
            for (l=0; l < n_var; l++)
              if (INTEGER(restrictQR)[i2 + l*n_var] || INTEGER(restrictQR)[j2 + l*n_var])
                restrictQ[n_rQ++] = l;
          }
        }

        nrr[UTE2I(i2, j2)] = n_I == 0 ? qp_edge_nrr(REAL(XR), S, n_var, N, i2, j2, q,
                                                    restrictQ, n_rQ, fixQ, n_fQ, nTests, alpha) :
                                        qp_edge_nrr_hmgm(REAL(XR), n_var, N, I, n_I,
                                                         INTEGER(n_levelsR), Y, n_Y,
                                                         ssdMat, mapX2ssd, i2, j2, q,
                                                         restrictQ, n_rQ, fixQ, n_fQ,
                                                         nTests, alpha, INTEGER(exactTest)[0]);
        k++;
        if (startTime > 0 && k == nAdjEtime)
          break;
        pct = (int) ((k * 100) / n_adj);
        if (pct != ppct) {
          if (verbose && startTime == 0) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
            SETCAR(t, pb);
            SET_TAG(t, install("pb")); t=CDR(t);
            SETCAR(t, ScalarReal(((double) pct) / 100.0));
            SET_TAG(t, install("value"));
            eval(s, R_GlobalEnv);
            UNPROTECT(1); /* t s */
          }
          R_CheckUserInterrupt();
#ifdef Win32
          R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
          R_ProcessEvents();
#endif
          ppct = pct;
        }
      }
      if (startTime > 0 && k == nAdjEtime)
        break;
    }
  }

  /* intersection variables against themselves (avoiding pairing the same) */
  if (startTime == 0 || k < nAdjEtime) {
    for (i = 0; i < l_int-1; i++) {
      int i2 = pairup_ij_int[i] - 1;

      for (j = i+1; j < l_int; j++) {
        int j2 = pairup_ij_int[j] - 1;

        if (restrictQ != NULL) {
          if (isMatrix_restrictQ) {
            int l;

            n_rQ = 0;
            for (l=0; l < n_var; l++)
              if (INTEGER(restrictQR)[i2 + l*n_var] || INTEGER(restrictQR)[j2 + l*n_var])
                restrictQ[n_rQ++] = l;
          }
        }

        nrr[UTE2I(i2, j2)] = n_I == 0 ? qp_edge_nrr(REAL(XR), S, n_var, N, i2, j2, q,
                                                    restrictQ, n_rQ, fixQ, n_fQ, nTests, alpha) :
                                        qp_edge_nrr_hmgm(REAL(XR), n_var, N, I, n_I,
                                                         INTEGER(n_levelsR), Y, n_Y,
                                                         ssdMat, mapX2ssd, i2, j2, q,
                                                         restrictQ, n_rQ, fixQ, n_fQ,
                                                         nTests, alpha, INTEGER(exactTest)[0]);
        k++;
        if (startTime > 0 && k == nAdjEtime)
          break;
        pct = (int) ((k * 100) / n_adj);
        if (pct != ppct) {
          if (verbose && startTime == 0) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
            SETCAR(t, pb);
            SET_TAG(t, install("pb")); t=CDR(t);
            SETCAR(t, ScalarReal(((double) pct) / 100.0));
            SET_TAG(t, install("value"));
            eval(s, R_GlobalEnv);
            UNPROTECT(1); /* t s */
          }
          R_CheckUserInterrupt();
#ifdef Win32
          R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
          R_ProcessEvents();
#endif
          ppct = pct;
        }
      }
      if (startTime > 0 && k == nAdjEtime)
        break;
    }
  }

  if (!work_with_margin)
    Free(S); /* = Free(ssdMat) */

  if (n_I > 0) {
    if (!work_with_margin)
      Free(mapX2ssd);
    Free(Y);
    Free(I);
  }

  if (restrictQR != R_NilValue)
    Free(restrictQ);

  if (fixQR != R_NilValue)
    Free(fixQ);

  if (verbose && startTime == 0) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("close")); t=CDR(t);
    SETCAR(t, pb);
    eval(s, R_GlobalEnv);
    UNPROTECT(2); /* t s pb */
  }

  UNPROTECT(1);   /* nrrR */

  if (startTime > 0) {
    SEXP procTimeR;
    double* procTime;
    SEXP nm;
    int* estimatedTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = elapsedTime + ((procTime[2] - startTime) / (double) k) * (double) n_adj;
    UNPROTECT(2); /* call procTimeR */

    PROTECT(nrrR = allocVector(INTSXP, 4));
    PROTECT(nm = allocVector(STRSXP, 4));
    estimatedTime = INTEGER(nrrR);
    estimatedTime[0] = (int) (elapsedTime / (24.0*3600.0));
    estimatedTime[1] = (int) ((elapsedTime - estimatedTime[0]*24.0*3600.0) / 3600.0);
    estimatedTime[2] = (int) ((elapsedTime - estimatedTime[0]*24.0*3600.0 -
                               estimatedTime[1]*3600.0) / 60.0);
    estimatedTime[3] = (int) (elapsedTime - estimatedTime[0]*24.0*3600.0 -
                              estimatedTime[1]*3600.0 - estimatedTime[2]*60.0 + 1.0);
    SET_STRING_ELT(nm, 0, mkChar("days"));
    SET_STRING_ELT(nm, 1, mkChar("hours"));
    SET_STRING_ELT(nm, 2, mkChar("minutes"));
    SET_STRING_ELT(nm, 3, mkChar("seconds"));
    setAttrib(nrrR, R_NamesSymbol, nm);

    UNPROTECT(2); /* nrrR nm */
  }

  return nrrR;
}



/*
  FUNCTION: qp_fast_nrr_identicalQs
  PURPOSE: compute for each pair of vertices indexed by the rows (columns)
           of the matrix S the non-rejection rate using a common set of Q sets
           for all vertex pairs considered. Vertex pairs may be restricted
           by using the pairup_* arguments
  RETURNS: matrix of non-rejection rate values in terms of number of non-rejected
           (accepted) tests for each pair of vertices
*/

static SEXP
qp_fast_nrr_identicalQs(SEXP XR, SEXP qR, SEXP restrictQR, SEXP fixQR, SEXP nTestsR,
                        SEXP alphaR, SEXP pairup_i_nointR, SEXP pairup_j_nointR,
                        SEXP pairup_ij_intR, SEXP verboseR, SEXP startTimeR,
                        SEXP nAdj2estimateTimeR, SEXP env) {
  int     N;
  int     n_var;
  int     q;
  int     nTests;
  double  alpha;
  double* S = NULL;
  int*    restrictQ = NULL;
  int     n_rQ = length(restrictQR);
  int*    fixQ = NULL;
  int     n_fQ = length(fixQR);
  int     l_ini = length(pairup_i_nointR);
  int     l_jni = length(pairup_j_nointR);
  int     l_int = length(pairup_ij_intR);
  int*    pairup_i_noint = INTEGER(pairup_i_nointR);
  int*    pairup_j_noint = INTEGER(pairup_j_nointR);
  int*    pairup_ij_int = INTEGER(pairup_ij_intR);
  int*    pairup_ij_noint = NULL;
  int     i,j,k;
  int     n_adj,pct,ppct;
  int*    q_by_T_samples;
  int*    Q;
  double* Qmat;
  double* Qinv;
  SEXP    nrrR;
  double* nrr;
  int     verbose;
  double  startTime, elapsedTime;
  int     nAdjEtime;
  SEXP    pb=NULL;

  N         = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  n_var     = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  q         = INTEGER(qR)[0];
  nTests    = INTEGER(nTestsR)[0];
  alpha     = REAL(alphaR)[0];
  verbose   = INTEGER(verboseR)[0];
  startTime = REAL(startTimeR)[0];
  nAdjEtime = INTEGER(nAdj2estimateTimeR)[0];

  if (q > n_var-2)
    error("q=%d > n.var-2=%d",q,n_var-2);

  if (q < 0)
    error("q=%d < 0",q);

  if (q > N-3)
    error("q=%d > n-3=%d", q, N-3);

  if (nTests < 1)
    error("nTests=%d < 1", nTests);

  if (alpha < 0.0 || alpha > 1.0)
    error("significance level alpha is %.2f and it should lie in the interval [0, 1]\n", alpha);

  if (missing_obs(REAL(XR), n_var, N, NULL, n_var, NULL, N))
    error("Missing values present in the data. The current setting identicalQs=TRUE speeds up calculations as long as data is complete. Please set identicalQs=FALSE or discard variables and/or observations with missing values from the input data.\n");

  S = Calloc((n_var*(n_var+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
  ssd(REAL(XR), n_var, N, NULL, n_var, NULL, N, TRUE, NULL, S);

  if (n_rQ > 0) {
    restrictQ = Calloc(n_rQ, int);
    for (i=0; i < n_rQ; i++)
      restrictQ[i] = INTEGER(restrictQR)[i] - 1;
  }

  if (n_fQ > 0) {
    fixQ = Calloc(n_fQ, int);
    for (i=0; i < n_fQ; i++)
      fixQ[i] = INTEGER(fixQR)[i] - 1;
  }

  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  /* sample the Q sets and pre-calculate the inverse matrices */

  q_by_T_samples = Calloc(q * nTests, int);

  if (n_rQ == 0)
    sampleQs(nTests, q, -1, -1, n_var, NULL, fixQ, n_fQ, q_by_T_samples);
  else
    sampleQs(nTests, q, -1, -1, n_rQ, restrictQ, fixQ, n_fQ, q_by_T_samples);

  Qmat = Calloc(q*q, double);
  Qinv = Calloc(q*q*nTests, double);

  /* ** DEPRECATED NON-TRIANGULAR COVARIANCE MATRIX MANIPULATION **
  for (i=0; i < nTests; i++) {
    Q = (int*) (q_by_T_samples+i*q);
    for (j=0; j < q; j++) {
      for (k=0; k < j; k++)
        Qmat[j + k*q] = Qmat[k + j*q] = S[Q[j] + Q[k] * n_var];
      Qmat[j + j*q] = S[Q[j] + Q[j] * n_var];
    }
    matinv((double*) (Qinv+i*q*q), Qmat, q);
  }
  */

  for (i=0; i < nTests; i++) {
    Q = (int*) (q_by_T_samples+i*q);
    for (j=0; j < q; j++) {
      for (k=0; k < j; k++)
        Qmat[j + k*q] = Qmat[k + j*q] = S[UTE2I(Q[j], Q[k])];
      Qmat[j + j*q] = S[UTE2I(Q[j], Q[j])];
    }
    matinv((double*) (Qinv+i*q*q), Qmat, q, 0);
  }
  Free(Qmat);
  
  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  PROTECT(nrrR = allocVector(REALSXP, (n_var*(n_var-1))/2+n_var));
  nrr = REAL(nrrR);

  for (i=0;i < (n_var*(n_var-1))/2+n_var;i++)
    nrr[i] = NA_REAL;

  n_adj = l_int * (l_jni + l_ini) + l_ini * l_jni + l_int * (l_int - 1) / 2;

  elapsedTime = 0.0;
  if (startTime > 0.0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = procTime[2] - startTime;
    startTime = procTime[2];
    UNPROTECT(2); /* call procTimeR */
  }

  ppct = -1;
  k = 0;
  if (verbose && startTime == 0) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("txtProgressBar")); t=CDR(t);
    SETCAR(t, ScalarInteger(3));
    SET_TAG(t, install("style"));
    PROTECT(pb = eval(s, R_GlobalEnv));
    UNPROTECT(1); /* t s */
  }

  /* intersection variables against ij-exclusive variables */
  for (i=0; i < l_int; i++) {
    int i2 = pairup_ij_int[i] - 1;

    for (j=0; j < l_ini + l_jni; j++) {
      int j2 = pairup_ij_noint[j] - 1;

      nrr[UTE2I(i2, j2)] = qp_edge_nrr_identicalQs(S, n_var, q_by_T_samples, Qinv,
                                                   N, i2, j2, q, nTests, alpha);
      k++;
      if (startTime > 0 && k == nAdjEtime)
        break;
      pct = (int) ((k * 100) / n_adj);
      if (pct != ppct) {
        if (verbose && startTime == 0) {
          SEXP s, t;
          PROTECT(t = s = allocList(3));
          SET_TYPEOF(s, LANGSXP);
          SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
          SETCAR(t, pb);
          SET_TAG(t, install("pb")); t=CDR(t);
          SETCAR(t, ScalarReal(((double) pct) / 100.0));
          SET_TAG(t, install("value"));
          eval(s, R_GlobalEnv);
          UNPROTECT(1); /* t s */
        }
        R_CheckUserInterrupt();
#ifdef Win32
        R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
    if (startTime > 0 && k == nAdjEtime)
      break;
  }

  if (l_ini + l_jni > 0)
    Free(pairup_ij_noint);

  /* i-exclusive variables against j-exclusive variables */
  if (startTime == 0 || k < nAdjEtime) {
    for (i=0; i < l_ini; i++) {
      int i2 = pairup_i_noint[i] - 1;

      for (j=0; j < l_jni; j++) {
        int j2 = pairup_j_noint[j] - 1;

        nrr[UTE2I(i2, j2)] = qp_edge_nrr_identicalQs(S, n_var, q_by_T_samples, Qinv,
                                                     N, i2, j2, q, nTests, alpha);
        k++;
        if (startTime > 0 && k == nAdjEtime)
          break;
        pct = (int) ((k * 100) / n_adj);
        if (pct != ppct) {
          if (verbose && startTime == 0) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
            SETCAR(t, pb);
            SET_TAG(t, install("pb")); t=CDR(t);
            SETCAR(t, ScalarReal(((double) pct) / 100.0));
            SET_TAG(t, install("value"));
            eval(s, R_GlobalEnv);
            UNPROTECT(1); /* t s */
          }
          R_CheckUserInterrupt();
#ifdef Win32
          R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
          R_ProcessEvents();
#endif
          ppct = pct;
        }
      }
      if (startTime > 0 && k == nAdjEtime)
        break;
    }
  }

  /* intersection variables against themselves (avoiding pairing the same) */
  if (startTime == 0 || k < nAdjEtime) {
    for (i = 0; i < l_int-1; i++) {
      int i2 = pairup_ij_int[i] - 1;

      for (j = i+1; j < l_int; j++) {
        int j2 = pairup_ij_int[j] - 1;

        nrr[UTE2I(i2, j2)] = qp_edge_nrr_identicalQs(S, n_var, q_by_T_samples, Qinv,
                                                     N, i2, j2, q, nTests, alpha);
        k++;
        if (startTime > 0 && k == nAdjEtime)
          break;
        pct = (int) ((k * 100) / n_adj);
        if (pct != ppct) {
          if (verbose && startTime == 0) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
            SETCAR(t, pb);
            SET_TAG(t, install("pb")); t=CDR(t);
            SETCAR(t, ScalarReal(((double) pct) / 100.0));
            SET_TAG(t, install("value"));
            eval(s, R_GlobalEnv);
            UNPROTECT(1); /* t s */
          }
          R_CheckUserInterrupt();
#ifdef Win32
          R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
          R_ProcessEvents();
#endif
          ppct = pct;
        }
      }
      if (startTime > 0 && k == nAdjEtime)
        break;
    }
  }

  Free(S);
  Free(Qinv);

  if (restrictQR != R_NilValue)
    Free(restrictQ);

  if (fixQR != R_NilValue)
    Free(fixQ);

  if (verbose && startTime == 0) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("close")); t=CDR(t);
    SETCAR(t, pb);
    eval(s, R_GlobalEnv);
    UNPROTECT(2); /* t s pb */
  }

  UNPROTECT(1);   /* nrrR */

  if (startTime > 0) {
    SEXP procTimeR;
    double* procTime;
    SEXP nm;
    int* estimatedTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = elapsedTime + ((procTime[2] - startTime) / (double) k) * (double) n_adj;
    UNPROTECT(2); /* call procTimeR */

    PROTECT(nrrR = allocVector(INTSXP, 4));
    PROTECT(nm = allocVector(STRSXP, 4));
    estimatedTime = INTEGER(nrrR);
    estimatedTime[0] = (int) (elapsedTime / (24.0*3600.0));
    estimatedTime[1] = (int) ((elapsedTime - estimatedTime[0]*24.0*3600.0) / 3600.0);
    estimatedTime[2] = (int) ((elapsedTime - estimatedTime[0]*24.0*3600.0 -
                               estimatedTime[1]*3600.0) / 60.0);
    estimatedTime[3] = (int) (elapsedTime - estimatedTime[0]*24.0*3600.0 -
                              estimatedTime[1]*3600.0 - estimatedTime[2]*60.0 + 1.0);
    SET_STRING_ELT(nm, 0, mkChar("days"));
    SET_STRING_ELT(nm, 1, mkChar("hours"));
    SET_STRING_ELT(nm, 2, mkChar("minutes"));
    SET_STRING_ELT(nm, 3, mkChar("seconds"));
    setAttrib(nrrR, R_NamesSymbol, nm);

    UNPROTECT(2); /* nrrR nm */
  }

  return nrrR;
}



/*
  FUNCTION: qp_fast_nrr_par
  PURPOSE: compute for each pair of vertices indexed by the rows (columns)
           of the matrix S the non-rejection rate. Vertex pairs may be restricted
           by using the pairup_* arguments. This function should be called only
           within a parallel environment running in a cluster where arguments
           myRankR and clSzeR tell how many nodes form the cluster (clSzeR) and
           which is the node running the function (myRankR)
  RETURNS: matrix of non-rejection rate values in terms of number of non-rejected
           (accepted) tests for each pair of vertices
*/

static SEXP
qp_fast_nrr_par(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP qR,
                SEXP restrictQR, SEXP fixQR, SEXP nTestsR, SEXP alphaR,
                SEXP pairup_i_nointR, SEXP pairup_j_nointR,
                SEXP pairup_ij_intR, SEXP exactTest, SEXP verboseR,
                SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP myRankR,
                SEXP clSzeR, SEXP masterNode, SEXP env) {
  int     N;
  int     n_var;
  int     q;
  int     nTests;
  double  alpha;
  double* S = NULL;
  double* ssdMat = NULL;
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int*    restrictQ = NULL;
  int     n_rQ = 0;
  int*    fixQ = NULL;
  int     n_fQ = 0;
  int     isMatrix_restrictQ = FALSE;
  int     l_ini = length(pairup_i_nointR);
  int     l_jni = length(pairup_j_nointR);
  int     l_int = length(pairup_ij_intR);
  int*    I = NULL;
  int*    Y = NULL;
  int*    mapX2ssd = NULL;
  int*    pairup_i_noint = INTEGER(pairup_i_nointR);
  int*    pairup_j_noint = INTEGER(pairup_j_nointR);
  int*    pairup_ij_int = INTEGER(pairup_ij_intR);
  int*    pairup_ij_noint = NULL;
  int     i,j,k,n_adj,n_adj_this_proc,pct,ppct;
  SEXP    nrrR, idxR;
  SEXP    result, result_names;
  double* nrr;
  int*    idx;
  int     verbose;
  int     myrank;
  int     clsze;
  int     firstAdj, lastAdj;
  double  startTime, elapsedTime;
  int     nAdjEtime;
  SEXP    progressReport,progressReportType,
          progressReportValue,progressReportSuccess,
          progressReportTag,progressReport_names;

  PROTECT(progressReport = allocVector(VECSXP,4));
  SET_VECTOR_ELT(progressReport,0,progressReportType = allocVector(STRSXP,1));
  SET_VECTOR_ELT(progressReport,1,progressReportValue = allocVector(INTSXP,1));
  SET_VECTOR_ELT(progressReport,2,progressReportSuccess = allocVector(LGLSXP,1));
  SET_VECTOR_ELT(progressReport,3,progressReportTag = allocVector(STRSXP,1));
  PROTECT(progressReport_names = allocVector(STRSXP,4));
  SET_STRING_ELT(progressReport_names,0,mkChar("type"));
  SET_STRING_ELT(progressReport_names,1,mkChar("value"));
  SET_STRING_ELT(progressReport_names,2,mkChar("success"));
  SET_STRING_ELT(progressReport_names,3,mkChar("tag"));
  setAttrib(progressReport,R_NamesSymbol,progressReport_names);
  SET_STRING_ELT(VECTOR_ELT(progressReport,0), 0, mkChar("VALUE"));
  INTEGER(VECTOR_ELT(progressReport,1))[0] = 0;
  LOGICAL(VECTOR_ELT(progressReport,2))[0] = TRUE;
  SET_STRING_ELT(VECTOR_ELT(progressReport,3), 0, mkChar("UPDATE"));

  N         = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  n_var     = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  q         = INTEGER(qR)[0];
  nTests    = INTEGER(nTestsR)[0];
  alpha     = REAL(alphaR)[0];
  verbose   = INTEGER(verboseR)[0];
  startTime = REAL(startTimeR)[0];
  nAdjEtime = INTEGER(nAdj2estimateTimeR)[0];
  myrank    = INTEGER(myRankR)[0];
  clsze     = INTEGER(clSzeR)[0];

  if (q > n_var-2)
    error("q=%d > p-2=%d",q,n_var-2);

  if (q < 0)
    error("q=%d < 0",q);

  if (q > N-3)
    error("q=%d > n-3=%d", q, N-3);

  if (nTests < 1)
    error("nTests=%d < 1", nTests);

  if (alpha < 0.0 || alpha > 1.0)
    error("significance level alpha is %.2f and it should lie in the interval [0, 1]\n", alpha);

  if (n_I == 0) {
    if (!missing_obs(REAL(XR), n_var, N, NULL, n_var, NULL, N)) {
      S = ssdMat = Calloc((n_var*(n_var+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, N, NULL, n_var, NULL, N, TRUE, NULL, S);
    }
  } else {
    I = Calloc(n_I, int);
    for (i=0; i < n_I; i++)
      I[i] = INTEGER(IR)[i]-1;

    Y = Calloc(n_Y, int);
    for (i=0; i < n_Y; i++)
      Y[i] = INTEGER(YR)[i]-1;

    if (!missing_obs(REAL(XR), n_var, N, NULL, n_var, NULL, N)) {
      mapX2ssd = Calloc(n_var, int);
      for (i=0; i < n_var; i++) {
        j = 0;
        while (j < n_Y && i != Y[j])
          j++;

        mapX2ssd[i] = j;
      }

      S = ssdMat = Calloc((n_Y*(n_Y+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, N, Y, n_Y, NULL, N, FALSE, NULL, ssdMat);
    }
  }

  if (restrictQR != R_NilValue) {
    if (isMatrix(restrictQR)) {
      isMatrix_restrictQ = TRUE;
      restrictQ = Calloc(n_var, int);
    } else {
      n_rQ = length(restrictQR);
      restrictQ = Calloc(n_rQ, int);
      for (i=0; i < n_rQ; i++)
        restrictQ[i] = INTEGER(restrictQR)[i] - 1;
    }
  }

  if (fixQR != R_NilValue) {
    n_fQ = length(fixQR);
    fixQ = Calloc(n_fQ, int);
    for (i=0; i < n_fQ; i++)
      fixQ[i] = INTEGER(fixQR)[i] - 1;
  }

  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  n_adj = l_int * (l_jni + l_ini) + l_ini * l_jni + l_int * (l_int - 1) / 2;

  firstAdj = (myrank-1) * (n_adj / clsze);
  lastAdj  = myrank * (n_adj / clsze);

  if (myrank == clsze)
    lastAdj += n_adj - lastAdj;

  lastAdj--;

  n_adj_this_proc = lastAdj - firstAdj + 1;

  PROTECT(result = allocVector(VECSXP,2));
  SET_VECTOR_ELT(result, 0, nrrR = allocVector(REALSXP, lastAdj-firstAdj+1));
  SET_VECTOR_ELT(result, 1, idxR = allocVector(INTSXP, lastAdj-firstAdj+1));
  PROTECT(result_names = allocVector(STRSXP, 2));
  SET_STRING_ELT(result_names, 0, mkChar("nrr"));
  SET_STRING_ELT(result_names, 1, mkChar("idx"));
  setAttrib(result, R_NamesSymbol, result_names);
  nrr = REAL(VECTOR_ELT(result, 0));
  idx = INTEGER(VECTOR_ELT(result, 1));

  elapsedTime = 0.0;
  if (startTime > 0.0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    /* initialize 'idx' so that the R code copying the result works as
     * in a normal execution */
    for (k=0; k < n_adj_this_proc; k++)
      idx[k] = firstAdj + k + 1;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = procTime[2] - startTime;
    startTime = procTime[2];
    UNPROTECT(2); /* call procTimeR */
  }

  k = firstAdj;
  ppct = -1;

  if (k < l_int * (l_ini + l_jni)) {
    int j_first = k % (l_ini + l_jni);

    /* intersection variables against ij-exclusive variables */
    for (i=((int) (k/(l_ini + l_jni))); i < l_int && k <= lastAdj; i++) {
      int i2 = pairup_ij_int[i] - 1;

      for (j=j_first; j < l_ini + l_jni && k <= lastAdj; j++) {
        int j2 = pairup_ij_noint[j] - 1;

        if (restrictQ != NULL) {
          if (isMatrix_restrictQ) {
            int m;

            n_rQ = 0;
            for (m=0; m < n_var; m++)
              if (INTEGER(restrictQR)[i2 + m*n_var] || INTEGER(restrictQR)[j2 + m*n_var])
                restrictQ[n_rQ++] = m;
          }
        }

        nrr[k-firstAdj] = n_I == 0 ? qp_edge_nrr(REAL(XR), S, n_var, N, i2, j2, q, restrictQ,
                                                 n_rQ, fixQ, n_fQ, nTests, alpha) :
                                     qp_edge_nrr_hmgm(REAL(XR), n_var, N, I, n_I,
                                                      INTEGER(n_levelsR), Y, n_Y,
                                                      ssdMat, mapX2ssd, i2, j2, q,
                                                      restrictQ, n_rQ, fixQ, n_fQ,
                                                      nTests, alpha, INTEGER(exactTest)[0]);
        idx[k-firstAdj] = UTE2I(i2, j2) + 1;
        k++;
        if (startTime > 0 && k-firstAdj == nAdjEtime)
          break;
        if (verbose && startTime == 0) {
          pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
          if (pct != ppct) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("sendData")); t=CDR(t);
            SETCAR(t, masterNode);
            SET_TAG(t, install("node")); t=CDR(t);
            INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
            SETCAR(t, progressReport);
            SET_TAG(t, install("data"));
            eval(s, env);
            UNPROTECT(1); /* t s */
          }
          ppct = pct;
        }
      }
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      j_first = 0;
    }
  }

  if (l_ini + l_jni > 0)
    Free(pairup_ij_noint);

  if (k <= lastAdj && k < l_int * (l_ini + l_jni) + l_ini * l_jni &&
      (startTime == 0 || k-firstAdj < nAdjEtime)) {
    int i_first = ((int) ((k - l_int * (l_ini + l_jni)) / l_jni));
    int j_first = (k - l_int * (l_ini + l_jni)) % l_jni;

    /* i-exclusive variables against j-exclusive variables */
    for (i=i_first; i < l_ini && k <= lastAdj; i++) {
      int i2 = pairup_i_noint[i] - 1;

      for (j=j_first; j < l_jni && k <= lastAdj; j++) {
        int j2 = pairup_j_noint[j] - 1;

        if (restrictQ != NULL) {
          if (isMatrix_restrictQ) {
            int m;

            n_rQ = 0;
            for (m=0; m < n_var; m++)
              if (INTEGER(restrictQR)[i2 + m*n_var] || INTEGER(restrictQR)[j2 + m*n_var])
                restrictQ[n_rQ++] = m;
          }
        }

        nrr[k-firstAdj] = n_I == 0 ? qp_edge_nrr(REAL(XR), S, n_var, N, i2, j2, q, restrictQ,
                                                 n_rQ, fixQ, n_fQ, nTests, alpha) :
                                     qp_edge_nrr_hmgm(REAL(XR), n_var, N, I, n_I,
                                                      INTEGER(n_levelsR), Y, n_Y,
                                                      ssdMat, mapX2ssd, i2, j2, q,
                                                      restrictQ, n_rQ, fixQ, n_fQ,
                                                      nTests, alpha, INTEGER(exactTest)[0]);
        idx[k-firstAdj] = UTE2I(i2, j2) + 1;
        k++;
        if (startTime > 0 && k-firstAdj == nAdjEtime)
          break;
        if (verbose && startTime == 0) {
          pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
          if (pct != ppct) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("sendData")); t=CDR(t);
            SETCAR(t, masterNode);
            SET_TAG(t, install("node")); t=CDR(t);
            INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
            SETCAR(t, progressReport);
            SET_TAG(t, install("data"));
            eval(s, env);
            UNPROTECT(1); /* t s */
          }
          ppct = pct;
        }
      }
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      j_first = 0;
    }
  }

  if (k <= lastAdj && (startTime == 0 || k-firstAdj < nAdjEtime)) {
    int i_first = k - l_int * (l_ini + l_jni) - l_ini * l_jni;
    int l;

    /* intersection variables against themselves (avoiding pairing the same) */
    for (l = i_first; l < (l_int * (l_int - 1)) / 2 && k <= lastAdj; l++) {
      int i,j,i2,j2;
      i2e(l, &i, &j);

      i2 = pairup_ij_int[i] - 1;
      j2 = pairup_ij_int[j] - 1;

      if (restrictQ != NULL) {
        if (isMatrix_restrictQ) {
          int m;

          n_rQ = 0;
          for (m=0; m < n_var; m++)
            if (INTEGER(restrictQR)[i2 + m*n_var] || INTEGER(restrictQR)[j2 + m*n_var])
              restrictQ[n_rQ++] = m;
        }
      }

      nrr[k-firstAdj] = n_I == 0 ?  qp_edge_nrr(REAL(XR), S, n_var, N, i2, j2, q, restrictQ,
                                                n_rQ, fixQ, n_fQ, nTests, alpha) :
                                    qp_edge_nrr_hmgm(REAL(XR), n_var, N, I, n_I,
                                                     INTEGER(n_levelsR), Y, n_Y,
                                                     ssdMat, mapX2ssd, i2, j2, q,
                                                     restrictQ, n_rQ, fixQ, n_fQ,
                                                     nTests, alpha, INTEGER(exactTest)[0]);
      idx[k-firstAdj] = UTE2I(i2, j2) + 1;
      k++;
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      if (verbose && startTime == 0) {
        pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
        if (pct != ppct) {
          SEXP s, t;
          PROTECT(t = s = allocList(3));
          SET_TYPEOF(s, LANGSXP);
          SETCAR(t, install("sendData")); t=CDR(t);
          SETCAR(t, masterNode);
          SET_TAG(t, install("node")); t=CDR(t);
          INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
          SETCAR(t, progressReport);
          SET_TAG(t, install("data"));
          eval(s, env);
          UNPROTECT(1); /* t s */
        }
        ppct = pct;
      }
    }
  }

  Free(S); /* = Free(ssdMat) */

  if (n_I > 0) {
    if (restrictQR != R_NilValue)
      Free(restrictQ);
    Free(mapX2ssd);
    Free(Y);
    Free(I);
  }

  if (restrictQR != R_NilValue)
    Free(restrictQ);

  if (fixQR != R_NilValue)
    Free(fixQ);

  if (startTime > 0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = elapsedTime + ((procTime[2] - startTime) / (double) (k-firstAdj)) * (double) n_adj_this_proc;
    UNPROTECT(2); /* call procTimeR */

    nrr[0] = elapsedTime; /* store in the first position of the nrr vector the estimated time */
  }

  UNPROTECT(4);   /* result result_names progressReport progressReport_names */

  return result;
}



/*
  FUNCTION: qp_fast_nrr_identicalQs_par
  PURPOSE: compute for each pair of vertices indexed by the rows (columns)
           of the matrix S the non-rejection rate using a common set of Q sets
           for all vertex pairs considered. Vertex pairs may be restricted
           by using the pairup_* arguments
  RETURNS: matrix of non-rejection rate values in terms of number of non-rejected
           (accepted) tests for each pair of vertices
*/

static SEXP
qp_fast_nrr_identicalQs_par(SEXP XR, SEXP qR, SEXP restrictQR, SEXP fixQR,
                            SEXP nTestsR, SEXP alphaR, SEXP pairup_i_nointR,
                            SEXP pairup_j_nointR, SEXP pairup_ij_intR,
                            SEXP verboseR, SEXP startTimeR, SEXP nAdj2estimateTimeR,
                            SEXP myRankR, SEXP clSzeR, SEXP masterNode, SEXP env) {
  int     N;
  int     n_var;
  int     q;
  int     nTests;
  double  alpha;
  double* S;
  int*    restrictQ = NULL;
  int     n_rQ = length(restrictQR);
  int*    fixQ = NULL;
  int     n_fQ = length(fixQR);
  int     l_ini = length(pairup_i_nointR);
  int     l_jni = length(pairup_j_nointR);
  int     l_int = length(pairup_ij_intR);
  int*    pairup_i_noint = INTEGER(pairup_i_nointR);
  int*    pairup_j_noint = INTEGER(pairup_j_nointR);
  int*    pairup_ij_int = INTEGER(pairup_ij_intR);
  int*    pairup_ij_noint = NULL;
  int     i,j,k,n_adj,n_adj_this_proc,pct,ppct;
  int*    q_by_T_samples;
  int*    Q;
  double* Qmat;
  double* Qinv;
  SEXP    nrrR, idxR;
  SEXP    result, result_names;
  double* nrr;
  int*    idx;
  int     verbose;
  int     myrank;
  int     clsze;
  int     firstAdj, lastAdj;
  double  startTime, elapsedTime;
  int     nAdjEtime;
  SEXP    progressReport,progressReportType,
          progressReportValue,progressReportSuccess,
          progressReportTag,progressReport_names;

  PROTECT(progressReport = allocVector(VECSXP,4));
  SET_VECTOR_ELT(progressReport,0,progressReportType = allocVector(STRSXP,1));
  SET_VECTOR_ELT(progressReport,1,progressReportValue = allocVector(INTSXP,1));
  SET_VECTOR_ELT(progressReport,2,progressReportSuccess = allocVector(LGLSXP,1));
  SET_VECTOR_ELT(progressReport,3,progressReportTag = allocVector(STRSXP,1));
  PROTECT(progressReport_names = allocVector(STRSXP,4));
  SET_STRING_ELT(progressReport_names,0,mkChar("type"));
  SET_STRING_ELT(progressReport_names,1,mkChar("value"));
  SET_STRING_ELT(progressReport_names,2,mkChar("success"));
  SET_STRING_ELT(progressReport_names,3,mkChar("tag"));
  setAttrib(progressReport,R_NamesSymbol,progressReport_names);
  SET_STRING_ELT(VECTOR_ELT(progressReport,0), 0, mkChar("VALUE"));
  INTEGER(VECTOR_ELT(progressReport,1))[0] = 0;
  LOGICAL(VECTOR_ELT(progressReport,2))[0] = TRUE;
  SET_STRING_ELT(VECTOR_ELT(progressReport,3), 0, mkChar("UPDATE"));

  N         = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  n_var     = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  q         = INTEGER(qR)[0];
  nTests    = INTEGER(nTestsR)[0];
  alpha     = REAL(alphaR)[0];
  verbose   = INTEGER(verboseR)[0];
  startTime = REAL(startTimeR)[0];
  nAdjEtime = INTEGER(nAdj2estimateTimeR)[0];
  myrank    = INTEGER(myRankR)[0];
  clsze     = INTEGER(clSzeR)[0];

  if (q > n_var-2)
    error("q=%d > n.var-2=%d",q,n_var-2);

  if (q < 0)
    error("q=%d < 0",q);

  if (q > N-3)
    error("q=%d > N-3=%d", q, N-3);

  if (missing_obs(REAL(XR), n_var, N, NULL, n_var, NULL, N))
    error("Missing values present in the data. The current setting identicalQs=TRUE speeds up calculations as long as data is complete. Please set identicalQs=FALSE or discard variables and/or observations with missing values from the input data.\n");

  S = Calloc((n_var*(n_var+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
  ssd(REAL(XR), n_var, N, NULL, n_var, NULL, N, TRUE, NULL, S);

  if (n_rQ > 0) {
    restrictQ = Calloc(n_rQ, int);
    for (i=0; i < n_rQ; i++)
      restrictQ[i] = INTEGER(restrictQR)[i] - 1;
  }

  if (n_fQ > 0) {
    fixQ = Calloc(n_fQ, int);
    for (i=0; i < n_fQ; i++)
      fixQ[i] = INTEGER(fixQR)[i] - 1;
  }

  /* sample the Q sets and pre-calculate the inverse matrices */

  q_by_T_samples = Calloc(q * nTests, int);

  if (n_rQ == 0)
    sampleQs(nTests, q, -1, -1, n_var, NULL, fixQ, n_fQ, q_by_T_samples);
  else
    sampleQs(nTests, q, -1, -1, n_rQ, restrictQ, fixQ, n_fQ, q_by_T_samples);


  Qmat = Calloc(q*q, double);
  Qinv = Calloc(q*q*nTests, double);

  /* ** DEPRECATED NON-TRIANGULAR COVARIANCE MATRIX MANIPULATION **
  for (i=0; i < nTests; i++) {
    Q = (int*) (q_by_T_samples+i*q);
    for (j=0; j < q; j++) {
      for (k=0; k < j; k++)
        Qmat[j + k*q] = Qmat[k + j*q] = S[Q[j] + Q[k] * n_var];
      Qmat[j + j*q] = S[Q[j] + Q[j] * n_var];
    }
    matinv((double*) (Qinv+i*q*q), Qmat, q);
  }
  */
  for (i=0; i < nTests; i++) {
    Q = (int*) (q_by_T_samples+i*q);
    for (j=0; j < q; j++) {
      for (k=0; k < j; k++)
        Qmat[j + k*q] = Qmat[k + j*q] = S[UTE2I(Q[j], Q[k])];
      Qmat[j + j*q] = S[UTE2I(Q[j], Q[j])];
    }
    matinv((double*) (Qinv+i*q*q), Qmat, q, 0);
  }
  Free(Qmat);
  
  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  n_adj = l_int * (l_jni + l_ini) + l_ini * l_jni + l_int * (l_int - 1) / 2;

  firstAdj = (myrank-1) * (n_adj / clsze);
  lastAdj  = myrank * (n_adj / clsze);

  if (myrank == clsze)
    lastAdj += n_adj - lastAdj;

  lastAdj--;

  n_adj_this_proc = lastAdj - firstAdj + 1;

  PROTECT(result = allocVector(VECSXP,2));
  SET_VECTOR_ELT(result, 0, nrrR = allocVector(REALSXP, lastAdj-firstAdj+1));
  SET_VECTOR_ELT(result, 1, idxR = allocVector(INTSXP, lastAdj-firstAdj+1));
  PROTECT(result_names = allocVector(STRSXP, 2));
  SET_STRING_ELT(result_names, 0, mkChar("nrr"));
  SET_STRING_ELT(result_names, 1, mkChar("idx"));
  setAttrib(result, R_NamesSymbol, result_names);
  nrr = REAL(VECTOR_ELT(result, 0));
  idx = INTEGER(VECTOR_ELT(result, 1));

  elapsedTime = 0.0;
  if (startTime > 0.0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    /* initialize 'idx' so that the R code copying the result works as
     * in a normal execution */
    for (k=0; k < n_adj_this_proc; k++)
      idx[k] = firstAdj + k + 1;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = procTime[2] - startTime;
    startTime = procTime[2];
    UNPROTECT(2); /* call procTimeR */
  }

  k = firstAdj;
  ppct = -1;

  if (k < l_int * (l_ini + l_jni)) {
    int j_first = k % (l_ini + l_jni);

    /* intersection variables against ij-exclusive variables */
    for (i=((int) (k/(l_ini + l_jni))); i < l_int && k <= lastAdj; i++) {
      int i2 = pairup_ij_int[i] - 1;

      for (j=j_first; j < l_ini + l_jni && k <= lastAdj; j++) {
        int j2 = pairup_ij_noint[j] - 1;

        nrr[k-firstAdj] = qp_edge_nrr_identicalQs(S, n_var, q_by_T_samples, Qinv,
                                                  N, i2, j2, q, nTests, alpha);
        idx[k-firstAdj] = UTE2I(i2, j2) + 1;
        k++;
        if (startTime > 0 && k-firstAdj == nAdjEtime)
          break;
        if (verbose && startTime == 0) {
          pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
          if (pct != ppct) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("sendData")); t=CDR(t);
            SETCAR(t, masterNode);
            SET_TAG(t, install("node")); t=CDR(t);
            INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
            SETCAR(t, progressReport);
            SET_TAG(t, install("data"));
            eval(s, env);
            UNPROTECT(1);
          }
          ppct = pct;
        }
      }
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      j_first = 0;
    }
  }

  if (l_ini + l_jni > 0)
    Free(pairup_ij_noint);

  if (k <= lastAdj && k < l_int * (l_ini + l_jni) + l_ini * l_jni &&
      (startTime == 0 || k-firstAdj < nAdjEtime)) {
    int i_first = ((int) ((k - l_int * (l_ini + l_jni)) / l_jni));
    int j_first = (k - l_int * (l_ini + l_jni)) % l_jni;

    /* i-exclusive variables against j-exclusive variables */
    for (i=i_first; i < l_ini && k <= lastAdj; i++) {
      int i2 = pairup_i_noint[i] - 1;

      for (j=j_first; j < l_jni && k <= lastAdj; j++) {
        int j2 = pairup_j_noint[j] - 1;

        nrr[k-firstAdj] = qp_edge_nrr_identicalQs(S, n_var, q_by_T_samples, Qinv,
                                                  N, i2, j2, q, nTests, alpha);
        idx[k-firstAdj] = UTE2I(i2, j2) + 1;
        k++;
        if (startTime > 0 && k-firstAdj == nAdjEtime)
          break;
        if (verbose && startTime == 0) {
          pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
          if (pct != ppct) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("sendData")); t=CDR(t);
            SETCAR(t, masterNode);
            SET_TAG(t, install("node")); t=CDR(t);
            INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
            SETCAR(t, progressReport);
            SET_TAG(t, install("data"));
            eval(s, env);
            UNPROTECT(1);
          }
          ppct = pct;
        }
      }
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      j_first = 0;
    }
  }

  if (k <= lastAdj && (startTime == 0 || k-firstAdj < nAdjEtime)) {
    int i_first = k - l_int * (l_ini + l_jni) - l_ini * l_jni;
    int l;

    /* intersection variables against themselves (avoiding pairing the same) */
    for (l = i_first; l < (l_int * (l_int - 1)) / 2 && k <= lastAdj; l++) {
      int i,j,i2,j2;
      i2e(l, &i, &j);

      i2 = pairup_ij_int[i] - 1;
      j2 = pairup_ij_int[j] - 1;

      nrr[k-firstAdj] = qp_edge_nrr_identicalQs(S, n_var, q_by_T_samples, Qinv,
                                                N, i2, j2, q, nTests, alpha);
      idx[k-firstAdj] = UTE2I(i2, j2) + 1;
      k++;
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      if (verbose && startTime == 0) {
        pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
        if (pct != ppct) {
          SEXP s, t;
          PROTECT(t = s = allocList(3));
          SET_TYPEOF(s, LANGSXP);
          SETCAR(t, install("sendData")); t=CDR(t);
          SETCAR(t, masterNode);
          SET_TAG(t, install("node")); t=CDR(t);
          INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
          SETCAR(t, progressReport);
          SET_TAG(t, install("data"));
          eval(s, env);
          UNPROTECT(1);
        }
        ppct = pct;
      }
    }
  }

  Free(S);
  Free(Qinv);

  if (restrictQR != R_NilValue)
    Free(restrictQ);

  if (fixQR != R_NilValue)
    Free(fixQ);

  if (startTime > 0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = elapsedTime + ((procTime[2] - startTime) / (double) (k-firstAdj)) * (double) n_adj_this_proc;
    UNPROTECT(2); /* call procTimeR */

    nrr[0] = elapsedTime; /* store in the first position of the nrr vector the estimated time */
  }

  UNPROTECT(4);   /* result result_names progressReport progressReport_names */

  return result;
}



/*
  FUNCTION: qp_fast_all_ci_tests
  PURPOSE: compute for each pair of vertices indexed by the columns of the
           matrix X one conditional independence test. Vertex pairs may be restricted
           by using the pairup_* arguments
  RETURNS: matrix of p-values of all the tests of conditional independence
*/ 

static SEXP
qp_fast_all_ci_tests(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP QR,
                     SEXP pairup_i_nointR, SEXP pairup_j_nointR, SEXP pairup_ij_intR,
                     SEXP exactTest, SEXP useR, SEXP tol, SEXP return_typeR, SEXP verboseR,
                     SEXP startTimeR, SEXP nAdj2estimateTimeR, SEXP env) {
  int     n, n_co;
  int     n_var;
  double* S = NULL;
  double* ssdMat = NULL;
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int*    ijQ = NULL;
  int*    Q = NULL;
  int     q = 0;
  int     n_upper_tri;
  int     work_with_margin = FALSE;
  int     l_ini = length(pairup_i_nointR);
  int     l_jni = length(pairup_j_nointR);
  int     l_int = length(pairup_ij_intR);
  int*    I = NULL;
  int*    Y = NULL;
  int*    mapX2ssd = NULL;
  int*    pairup_i_noint = INTEGER(pairup_i_nointR);
  int*    pairup_j_noint = INTEGER(pairup_j_nointR);
  int*    pairup_ij_int = INTEGER(pairup_ij_intR);
  int*    pairup_ij_noint = NULL;
  int     i,j,k;
  int     use, n_adj, pct, ppct;
  int     return_type;
  SEXP    result, result_names;
  int     result_n;
  SEXP    p_valuesR, statistic_valuesR, n_valuesR;
  double* p_values = NULL;
  double* statistic_values = NULL;
  int*    n_values = NULL;
  double  df, a, b, lambda;
  int     verbose;
  double  startTime, elapsedTime;
  int     nAdjEtime;
  SEXP    pb=NULL;

  n = n_co    = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  n_var       = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  verbose     = INTEGER(verboseR)[0];
  startTime   = REAL(startTimeR)[0];
  nAdjEtime   = INTEGER(nAdj2estimateTimeR)[0];
  use         = INTEGER(useR)[0];
  return_type = INTEGER(return_typeR)[0];

  if (n_I == 0) {
    if (!missing_obs(REAL(XR), n_var, n, NULL, n_var, NULL, n)) {
      S = ssdMat = Calloc((n_var*(n_var+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, n, NULL, n_var, NULL, n, TRUE, NULL, S);
    } else
      work_with_margin = TRUE;
  } else {
    I = Calloc(n_I, int);
    for (i=0; i < n_I; i++)
      I[i] = INTEGER(IR)[i]-1;

    Y = Calloc(n_Y, int);
    for (i=0; i < n_Y; i++)
      Y[i] = INTEGER(YR)[i]-1;

    if (!missing_obs(REAL(XR), n_var, n, NULL, n_var, NULL, n)) {
      mapX2ssd = Calloc(n_var, int);
      for (i=0; i < n_var; i++) {
        j = 0;
        while (j < n_Y && i != Y[j])
          j++;

        mapX2ssd[i] = j;
      }

      S = ssdMat = Calloc((n_Y*(n_Y+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, n, Y, n_Y, NULL, n, FALSE, NULL, ssdMat);
    } else
      work_with_margin = TRUE;
  }

  if (QR != R_NilValue) {
    q = length(QR);

    if (q > n_var-2)
      error("q=%d > p-2=%d", q, n_var-2);

    if (q < 0)
      error("q=%d < 0", q);

    if (q > n-3)
      error("q=%d > n-3=%d", q, n-3);

    if (work_with_margin)
      ijQ = Calloc(q+2, int); /* stores the indices of the variables in the {i, j, Q} margin */

    Q = Calloc(q, int);
    for (i=0; i < q; i++)
      if (work_with_margin) {
        ijQ[i+2] = INTEGER(QR)[i] - 1;
        Q[i] = n_I == 0 ? 2+i : INTEGER(QR)[i] - 1;
      } else
        Q[i] = INTEGER(QR)[i] - 1;
  } else
    ijQ = Calloc(2, int);

  n_upper_tri = ( (q+2) * ((q+2)+1) ) / 2; /* this upper triangle includes the diagonal */

  if (work_with_margin)
    S = Calloc(n_upper_tri, double);

  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  /* here we assume that RETURN_TYPE_PVALUE (1 elem), RETURN_TYPE_STATN (2 elem), RETURN_TYPE_ALL (3 elem) */
  result_n = 0;
  PROTECT(result = allocVector(VECSXP, return_type));
  PROTECT(result_names = allocVector(STRSXP, return_type));

  if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE) {
    SET_VECTOR_ELT(result, result_n, p_valuesR = allocVector(REALSXP, (n_var*(n_var+1))/2)); /* diagonal should be allocated */
    SET_STRING_ELT(result_names, result_n, mkChar("p.value"));
    p_values = REAL(VECTOR_ELT(result, result_n));
    result_n++;
  }

  if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
    SET_VECTOR_ELT(result, result_n, statistic_valuesR = allocVector(REALSXP, (n_var*(n_var+1))/2)); /* diagonal should be allocated */
    SET_STRING_ELT(result_names, result_n, mkChar("statistic"));
    statistic_values = REAL(VECTOR_ELT(result, result_n));
    result_n++;
    SET_VECTOR_ELT(result, result_n, n_valuesR = allocVector(INTSXP, (n_var*(n_var+1))/2));
    SET_STRING_ELT(result_names, result_n, mkChar("n"));
    n_values = INTEGER(VECTOR_ELT(result, result_n));
    result_n++;
  }

  setAttrib(result, R_NamesSymbol, result_names);

  for (i=0;i < (n_var*(n_var+1))/2;i++) {
    if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)  
      p_values[i] = NA_REAL;
    if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {  
      statistic_values[i] = NA_REAL;
      n_values[i] = NA_INTEGER;
    }
  }

  n_adj = l_int * (l_jni + l_ini) + l_ini * l_jni + l_int * (l_int - 1) / 2;

  elapsedTime = 0.0;
  if (startTime > 0.0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = procTime[2] - startTime;
    startTime = procTime[2];
    UNPROTECT(2); /* call procTimeR */
  }

  ppct = -1;
  k = 0;
  if (verbose && startTime == 0) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("txtProgressBar")); t=CDR(t);
    SETCAR(t, ScalarInteger(3));
    SET_TAG(t, install("style"));
    PROTECT(pb = eval(s, R_GlobalEnv));
    UNPROTECT(1); /* t s */
  }

  /* intersection variables against ij-exclusive variables */
  for (i=0; i < l_int; i++) {
    int i2 = pairup_ij_int[i] - 1;

    for (j=0; j < l_ini + l_jni; j++) {
      int j2 = pairup_ij_noint[j] - 1;

      if (n_I == 0) {
        if (work_with_margin) {
          ijQ[0] = i2;
          ijQ[1] = j2;
          memset(S, 0, sizeof(double) * n_upper_tri);
          n_co = ssd(REAL(XR), n_var, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
          lambda = qp_ci_test_std(S, q+2, n_co, 0, 1, Q, q, NULL);
        } else
          lambda = qp_ci_test_std(S, n_var, n, i2, j2, Q, q, NULL);
      } else
        lambda = qp_ci_test_hmgm(REAL(XR), n_var, n, I, n_I, INTEGER(n_levelsR),
                                 Y, n_Y, ssdMat, mapX2ssd, i2, j2,
                                 Q, q, use, REAL(tol)[0], &df, &a, &b, &n_co);

      if (n_I == 0) {
        if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
          p_values[UTE2I(i2, j2)] = 2.0 * (1.0 - pt(fabs(lambda), n_co-q-2, 1, 0));
        if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
          statistic_values[UTE2I(i2, j2)] = lambda;
          n_values[UTE2I(i2, j2)] = n_co;
        }
      } else {
        if (!ISNAN(lambda)) {
          if (INTEGER(exactTest)[0]) {
            lambda = exp(lambda / ((double) -(use == USE_EM ? n : n_co)));
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
              p_values[UTE2I(i2, j2)] = pbeta(lambda, a, b, TRUE, FALSE);
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
              statistic_values[UTE2I(i2, j2)] = lambda;
              n_values[UTE2I(i2, j2)] = use == USE_EM ? n : n_co;
            }
          } else {
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
              p_values[UTE2I(i2, j2)] = 1.0 - pchisq(lambda, df, TRUE, FALSE);
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
              statistic_values[UTE2I(i2, j2)] = lambda;
              n_values[UTE2I(i2, j2)] = use == USE_EM ? n : n_co;
            }
          }
        }
      }

      k++;
      if (startTime > 0 && k == nAdjEtime)
        break;
      pct = (int) ((k * 100) / n_adj);
      if (pct != ppct) {
        if (verbose && startTime == 0) {
          SEXP s, t;
          PROTECT(t = s = allocList(3));
          SET_TYPEOF(s, LANGSXP);
          SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
          SETCAR(t, pb);
          SET_TAG(t, install("pb")); t=CDR(t);
          SETCAR(t, ScalarReal(((double) pct) / 100.0));
          SET_TAG(t, install("value"));
          eval(s, R_GlobalEnv);
          UNPROTECT(1); /* t s */
        }
        R_CheckUserInterrupt();
#ifdef Win32
        R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
    if (startTime > 0 && k == nAdjEtime)
      break;
  }

  if (l_ini + l_jni > 0)
    Free(pairup_ij_noint);

  /* i-exclusive variables against j-exclusive variables */
  if (startTime == 0 || k < nAdjEtime) {
    for (i=0; i < l_ini; i++) {
      int i2 = pairup_i_noint[i] - 1;

      for (j=0; j < l_jni; j++) {
        int j2 = pairup_j_noint[j] - 1;

        if (n_I == 0) {
          if (work_with_margin) {
            ijQ[0] = i2;
            ijQ[1] = j2;
            memset(S, 0, sizeof(double) * n_upper_tri);
            n_co = ssd(REAL(XR), n_var, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
            lambda = qp_ci_test_std(S, q+2, n_co, 0, 1, Q, q, NULL);
          } else
            lambda = qp_ci_test_std(S, n_var, n, i2, j2, Q, q, NULL);
        } else
          lambda = qp_ci_test_hmgm(REAL(XR), n_var, n, I, n_I, INTEGER(n_levelsR),
                                   Y, n_Y, ssdMat, mapX2ssd, i2, j2,
                                   Q, q, use, REAL(tol)[0], &df, &a, &b, &n_co);

        if (n_I == 0) {
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
            p_values[UTE2I(i2, j2)] = 2.0 * (1.0 - pt(fabs(lambda), n_co-q-2, 1, 0));
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
            statistic_values[UTE2I(i2, j2)] = lambda;
            n_values[UTE2I(i2, j2)] = n_co;
          }
        } else {
          if (!ISNAN(lambda)) {
            if (INTEGER(exactTest)[0]) {
              lambda = exp(lambda / ((double) -(use == USE_EM ? n : n_co)));
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[UTE2I(i2, j2)] = pbeta(lambda, a, b, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[UTE2I(i2, j2)] = lambda;
                n_values[UTE2I(i2, j2)] = use == USE_EM ? n : n_co;
              }
            } else {
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[UTE2I(i2, j2)] = 1.0 - pchisq(lambda, df, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[UTE2I(i2, j2)] = lambda;
                n_values[UTE2I(i2, j2)] = use == USE_EM ? n : n_co;
              }
            }
          }
        }

        k++;
        if (startTime > 0 && k == nAdjEtime)
          break;
        pct = (int) ((k * 100) / n_adj);
        if (pct != ppct) {
          if (verbose && startTime == 0) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
            SETCAR(t, pb);
            SET_TAG(t, install("pb")); t=CDR(t);
            SETCAR(t, ScalarReal(((double) pct) / 100.0));
            SET_TAG(t, install("value"));
            eval(s, R_GlobalEnv);
            UNPROTECT(1); /* t s */
          }
          R_CheckUserInterrupt();
#ifdef Win32
          R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
          R_ProcessEvents();
#endif
          ppct = pct;
        }
      }
      if (startTime > 0 && k == nAdjEtime)
        break;
    }
  }

  /* intersection variables against themselves (avoiding pairing the same) */
  if (startTime == 0 || k < nAdjEtime) {
    for (i = 0; i < l_int-1; i++) {
      int i2 = pairup_ij_int[i] - 1;

      for (j = i+1; j < l_int; j++) {
        int j2 = pairup_ij_int[j] - 1;

        if (n_I == 0) {
          if (work_with_margin) {
            ijQ[0] = i2;
            ijQ[1] = j2;
            memset(S, 0, sizeof(double) * n_upper_tri);
            n_co = ssd(REAL(XR), n_var, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
            lambda = qp_ci_test_std(S, q+2, n_co, 0, 1, Q, q, NULL);
          } else
            lambda = qp_ci_test_std(S, n_var, n, i2, j2, Q, q, NULL);
        } else
          lambda = qp_ci_test_hmgm(REAL(XR), n_var, n, I, n_I, INTEGER(n_levelsR),
                                   Y, n_Y, ssdMat, mapX2ssd, i2, j2,
                                   Q, q, use, REAL(tol)[0], &df, &a, &b, &n_co);

        if (n_I == 0) {
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
            p_values[UTE2I(i2, j2)] = 2.0 * (1.0 - pt(fabs(lambda), n_co-q-2, 1, 0));
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
            statistic_values[UTE2I(i2, j2)] = lambda;
            n_values[UTE2I(i2, j2)] = n_co;
          }
        } else {
          if (!ISNAN(lambda)) {
            if (INTEGER(exactTest)[0]) {
              lambda = exp(lambda / ((double) -(use == USE_EM ? n : n_co)));
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[UTE2I(i2, j2)] = pbeta(lambda, a, b, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[UTE2I(i2, j2)] = lambda;
                n_values[UTE2I(i2, j2)] = use == USE_EM ? n : n_co;
              }
            } else {
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[UTE2I(i2, j2)] = 1.0 - pchisq(lambda, df, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[UTE2I(i2, j2)] = lambda;
                n_values[UTE2I(i2, j2)] = use == USE_EM ? n : n_co;
              }
            }
          }
        }

        k++;
        if (startTime > 0 && k == nAdjEtime)
          break;
        pct = (int) ((k * 100) / n_adj);
        if (pct != ppct) {
          if (verbose && startTime == 0) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
            SETCAR(t, pb);
            SET_TAG(t, install("pb")); t=CDR(t);
            SETCAR(t, ScalarReal(((double) pct) / 100.0));
            SET_TAG(t, install("value"));
            eval(s, R_GlobalEnv);
            UNPROTECT(1); /* t s */
          }
          R_CheckUserInterrupt();
#ifdef Win32
          R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
          R_ProcessEvents();
#endif
          ppct = pct;
        }
      }
      if (startTime > 0 && k == nAdjEtime)
        break;
    }
  }

  Free(S); /* = Free(ssdMat) */

  if (n_I > 0) {
    if (!work_with_margin)
      Free(mapX2ssd);
    Free(Y);
    Free(I);
  }

  if (work_with_margin)
    Free(ijQ);

  if (QR != R_NilValue)
    Free(Q);

  if (verbose && startTime == 0) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("close")); t=CDR(t);
    SETCAR(t, pb);
    eval(s, R_GlobalEnv);
    UNPROTECT(2); /* t s pb */
  }

  UNPROTECT(2);   /* result result_names */

  if (startTime > 0) {
    SEXP procTimeR;
    double* procTime;
    SEXP nm;
    int* estimatedTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = elapsedTime + ((procTime[2] - startTime) / (double) k) * (double) n_adj;
    UNPROTECT(2); /* call procTimeR */

    PROTECT(result = allocVector(INTSXP, 4));
    PROTECT(nm = allocVector(STRSXP, 4));
    estimatedTime = INTEGER(result);
    estimatedTime[0] = (int) (elapsedTime / (24.0*3600.0));
    estimatedTime[1] = (int) ((elapsedTime - estimatedTime[0]*24.0*3600.0) / 3600.0);
    estimatedTime[2] = (int) ((elapsedTime - estimatedTime[0]*24.0*3600.0 -
                               estimatedTime[1]*3600.0) / 60.0);
    estimatedTime[3] = (int) (elapsedTime - estimatedTime[0]*24.0*3600.0 -
                              estimatedTime[1]*3600.0 - estimatedTime[2]*60.0 + 1.0);
    SET_STRING_ELT(nm, 0, mkChar("days"));
    SET_STRING_ELT(nm, 1, mkChar("hours"));
    SET_STRING_ELT(nm, 2, mkChar("minutes"));
    SET_STRING_ELT(nm, 3, mkChar("seconds"));
    setAttrib(result, R_NamesSymbol, nm);

    UNPROTECT(2); /* result nm */
  }

  return result;
}



/*
  FUNCTION: qp_fast_all_ci_tests_par
  PURPOSE: compute for each pair of vertices indexed by the rows (columns)
           a conditional independence test. Vertex pairs may be restricted
           by using the pairup_* arguments. This function should be called only
           within a parallel environment running in a cluster where arguments
           myRankR and clSzeR tell how many nodes form the cluster (clSzeR) and
           which is the node running the function (myRankR)
  RETURNS: matrix of p-values of all the tests of conditional independence
*/

static SEXP
qp_fast_all_ci_tests_par(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP QR,
                         SEXP pairup_i_nointR, SEXP pairup_j_nointR,
                         SEXP pairup_ij_intR, SEXP exactTest, SEXP useR, SEXP tol,
                         SEXP return_typeR, SEXP verboseR, SEXP startTimeR,
                         SEXP nAdj2estimateTimeR, SEXP myRankR, SEXP clSzeR,
                         SEXP masterNode, SEXP env) {
  int     n, n_co;
  int     n_var;
  double* S = NULL;
  double* ssdMat = NULL;
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int*    ijQ = NULL;
  int*    Q = NULL;
  int     q = 0;
  int     n_upper_tri;
  int     work_with_margin = FALSE;
  int     l_ini = length(pairup_i_nointR);
  int     l_jni = length(pairup_j_nointR);
  int     l_int = length(pairup_ij_intR);
  int*    I = NULL;
  int*    Y = NULL;
  int*    mapX2ssd = NULL;
  int*    pairup_i_noint = INTEGER(pairup_i_nointR);
  int*    pairup_j_noint = INTEGER(pairup_j_nointR);
  int*    pairup_ij_int = INTEGER(pairup_ij_intR);
  int*    pairup_ij_noint = NULL;
  int     i,j,k;
  int     use, n_adj, n_adj_this_proc, pct, ppct;
  int     return_type;
  SEXP    result, result_names;
  int     result_n;
  SEXP    idxR, p_valuesR, statistic_valuesR, n_valuesR;
  int*    idx;
  double* p_values = NULL;
  double* statistic_values = NULL;
  int*    n_values = NULL;
  double  df, a, b, lambda;
  int     verbose;
  int     myrank;
  int     clsze;
  int     firstAdj, lastAdj;
  double  startTime, elapsedTime;
  int     nAdjEtime;
  SEXP    progressReport,progressReportType,
          progressReportValue,progressReportSuccess,
          progressReportTag,progressReport_names;

  PROTECT(progressReport = allocVector(VECSXP,4));
  SET_VECTOR_ELT(progressReport,0,progressReportType = allocVector(STRSXP,1));
  SET_VECTOR_ELT(progressReport,1,progressReportValue = allocVector(INTSXP,1));
  SET_VECTOR_ELT(progressReport,2,progressReportSuccess = allocVector(LGLSXP,1));
  SET_VECTOR_ELT(progressReport,3,progressReportTag = allocVector(STRSXP,1));
  PROTECT(progressReport_names = allocVector(STRSXP,4));
  SET_STRING_ELT(progressReport_names,0,mkChar("type"));
  SET_STRING_ELT(progressReport_names,1,mkChar("value"));
  SET_STRING_ELT(progressReport_names,2,mkChar("success"));
  SET_STRING_ELT(progressReport_names,3,mkChar("tag"));
  setAttrib(progressReport,R_NamesSymbol,progressReport_names);
  SET_STRING_ELT(VECTOR_ELT(progressReport,0), 0, mkChar("VALUE"));
  INTEGER(VECTOR_ELT(progressReport,1))[0] = 0;
  LOGICAL(VECTOR_ELT(progressReport,2))[0] = TRUE;
  SET_STRING_ELT(VECTOR_ELT(progressReport,3), 0, mkChar("UPDATE"));

  n = n_co    = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  n_var       = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  verbose     = INTEGER(verboseR)[0];
  startTime   = REAL(startTimeR)[0];
  nAdjEtime   = INTEGER(nAdj2estimateTimeR)[0];
  use         = INTEGER(useR)[0];
  return_type = INTEGER(return_typeR)[0];
  myrank      = INTEGER(myRankR)[0];
  clsze       = INTEGER(clSzeR)[0];

  if (n_I == 0) {
    if (!missing_obs(REAL(XR), n_var, n, NULL, n_var, NULL, n)) {
      S = ssdMat = Calloc((n_var*(n_var+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, n, NULL, n_var, NULL, n, TRUE, NULL, S);
    } else
      work_with_margin = TRUE;
  } else {
    I = Calloc(n_I, int);
    for (i=0; i < n_I; i++)
      I[i] = INTEGER(IR)[i]-1;

    Y = Calloc(n_Y, int);
    for (i=0; i < n_Y; i++)
      Y[i] = INTEGER(YR)[i]-1;

    if (!missing_obs(REAL(XR), n_var, n, NULL, n_var, NULL, n)) {
      mapX2ssd = Calloc(n_var, int);
      for (i=0; i < n_var; i++) {
        j = 0;
        while (j < n_Y && i != Y[j])
          j++;

        mapX2ssd[i] = j;
      }

      S = ssdMat = Calloc((n_Y*(n_Y+1))/2, double); /* if this doesn't do memset(0) there'll be trouble */
      ssd(REAL(XR), n_var, n, Y, n_Y, NULL, n, FALSE, NULL, ssdMat);
    } else
      work_with_margin = TRUE;
  }

  if (QR != R_NilValue) {
    q = length(QR);

    if (q > n_var-2)
      error("q=%d > p-2=%d", q, n_var-2);

    if (q < 0)
      error("q=%d < 0", q);

    if (q > n-3)
      error("q=%d > n-3=%d", q, n-3);

    if (work_with_margin)
      ijQ = Calloc(q+2, int);

    Q = Calloc(q, int);
    for (i=0; i < q; i++)
      if (work_with_margin) {
        ijQ[i+2] = INTEGER(QR)[i] - 1;
        Q[i] = n_I == 0 ? 2+i : INTEGER(QR)[i] - 1;
      } else
        Q[i] = INTEGER(QR)[i] - 1;
  } else
    ijQ = Calloc(2, int);

  n_upper_tri = ( (q+2) * ((q+2)+1) ) / 2; /* this upper triangle includes the diagonal */

  if (work_with_margin)
    S = Calloc(n_upper_tri, double);

  if (l_ini + l_jni > 0) {
    pairup_ij_noint = Calloc(l_ini + l_jni, int);
    Memcpy(pairup_ij_noint, pairup_i_noint, (size_t) l_ini);
    Memcpy(pairup_ij_noint + l_ini, pairup_j_noint, (size_t) l_jni);
  }

  n_adj = l_int * (l_jni + l_ini) + l_ini * l_jni + l_int * (l_int - 1) / 2;

  firstAdj = (myrank-1) * (n_adj / clsze);
  lastAdj  = myrank * (n_adj / clsze);

  if (myrank == clsze)
    lastAdj += n_adj - lastAdj;

  lastAdj--;

  n_adj_this_proc = lastAdj - firstAdj + 1;

  /* here we assume that RETURN_TYPE_PVALUE (1 elem), RETURN_TYPE_STATN (2 elem), RETURN_TYPE_ALL (3 elem) */
  result_n = 0;
  PROTECT(result = allocVector(VECSXP, return_type+1));
  PROTECT(result_names = allocVector(STRSXP, return_type+1));
  SET_VECTOR_ELT(result, result_n, idxR = allocVector(INTSXP, lastAdj-firstAdj+1));
  SET_STRING_ELT(result_names, result_n, mkChar("idx"));
  idx = INTEGER(VECTOR_ELT(result, result_n));
  result_n++;

  if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE || startTime > 0) {
    SET_VECTOR_ELT(result, result_n, p_valuesR = allocVector(REALSXP, lastAdj-firstAdj+1));
    SET_STRING_ELT(result_names, result_n, mkChar("p.value"));
    p_values = REAL(VECTOR_ELT(result, result_n));
    result_n++;
  }

  if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
    SET_VECTOR_ELT(result, result_n, statistic_valuesR = allocVector(REALSXP, lastAdj-firstAdj+1));
    SET_STRING_ELT(result_names, result_n, mkChar("statistic"));
    statistic_values = REAL(VECTOR_ELT(result, result_n));
    result_n++;
    SET_VECTOR_ELT(result, result_n, n_valuesR = allocVector(INTSXP, lastAdj-firstAdj+1));
    SET_STRING_ELT(result_names, result_n, mkChar("n"));
    n_values = INTEGER(VECTOR_ELT(result, result_n));
    result_n++;
  }

  setAttrib(result, R_NamesSymbol, result_names);

  elapsedTime = 0.0;
  if (startTime > 0.0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    /* initialize 'idx' so that the R code copying the result works as
     * in a normal execution */
    for (k=0; k < n_adj_this_proc; k++)
      idx[k] = firstAdj + k + 1;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = procTime[2] - startTime;
    startTime = procTime[2];
    UNPROTECT(2); /* call procTimeR */
  }

  k = firstAdj;
  ppct = -1;

  if (k < l_int * (l_ini + l_jni)) {
    int j_first = k % (l_ini + l_jni);

    /* intersection variables against ij-exclusive variables */
    for (i=((int) (k/(l_ini + l_jni))); i < l_int && k <= lastAdj; i++) {
      int i2 = pairup_ij_int[i] - 1;

      for (j=j_first; j < l_ini + l_jni && k <= lastAdj; j++) {
        int j2 = pairup_ij_noint[j] - 1;

        if (n_I == 0) {
          if (work_with_margin) {
            ijQ[0] = i2;
            ijQ[1] = j2;
            memset(S, 0, sizeof(double) * n_upper_tri);
            n_co = ssd(REAL(XR), n_var, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
            lambda = qp_ci_test_std(S, q+2, n_co, 0, 1, Q, q, NULL);
          } else
            lambda = qp_ci_test_std(S, n_var, n, i2, j2, Q, q, NULL);
        } else
          lambda = qp_ci_test_hmgm(REAL(XR), n_var, n, I, n_I, INTEGER(n_levelsR),
                                   Y, n_Y, ssdMat, mapX2ssd, i2, j2,
                                   Q, q, use, REAL(tol)[0], &df, &a, &b, &n_co);
        idx[k-firstAdj] = UTE2I(i2, j2) + 1;

        if (n_I == 0) {
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
            p_values[k-firstAdj] = 2.0 * (1.0 - pt(fabs(lambda), n_co-q-2, 1, 0));
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
            statistic_values[k-firstAdj] = lambda;
            n_values[k-firstAdj] = n_co;
          }
        } else {
          if (!ISNAN(lambda)) {
            if (INTEGER(exactTest)[0]) {
              lambda = exp(lambda / ((double) -(use == USE_EM ? n : n_co)));
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[k-firstAdj] = pbeta(lambda, a, b, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[k-firstAdj] = lambda;
                n_values[k-firstAdj] = use == USE_EM ? n : n_co;
              }
            } else {
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[k-firstAdj] = 1.0 - pchisq(lambda, df, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[k-firstAdj] = lambda;
                n_values[k-firstAdj] = use == USE_EM ? n : n_co;
              }
            }
          }
        }

        k++;
        if (startTime > 0 && k-firstAdj == nAdjEtime)
          break;
        if (verbose && startTime == 0) {
          pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
          if (pct != ppct) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("sendData")); t=CDR(t);
            SETCAR(t, masterNode);
            SET_TAG(t, install("node")); t=CDR(t);
            INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
            SETCAR(t, progressReport);
            SET_TAG(t, install("data"));
            eval(s, env);
            UNPROTECT(1); /* t s */
          }
          ppct = pct;
        }
      }
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      j_first = 0;
    }
  }

  if (l_ini + l_jni > 0)
    Free(pairup_ij_noint);

  if (k <= lastAdj && k < l_int * (l_ini + l_jni) + l_ini * l_jni &&
      (startTime == 0 || k-firstAdj < nAdjEtime)) {
    int i_first = ((int) ((k - l_int * (l_ini + l_jni)) / l_jni));
    int j_first = (k - l_int * (l_ini + l_jni)) % l_jni;

    /* i-exclusive variables against j-exclusive variables */
    for (i=i_first; i < l_ini && k <= lastAdj; i++) {
      int i2 = pairup_i_noint[i] - 1;

      for (j=j_first; j < l_jni && k <= lastAdj; j++) {
        int j2 = pairup_j_noint[j] - 1;

        if (n_I == 0) {
          if (work_with_margin) {
            ijQ[0] = i2;
            ijQ[1] = j2;
            memset(S, 0, sizeof(double) * n_upper_tri);
            n_co = ssd(REAL(XR), n_var, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
            lambda = qp_ci_test_std(S, q+2, n_co, 0, 1, Q, q, NULL);
          } else
            lambda = qp_ci_test_std(S, n_var, n, i2, j2, Q, q, NULL);
        } else
          lambda = qp_ci_test_hmgm(REAL(XR), n_var, n, I, n_I, INTEGER(n_levelsR),
                                   Y, n_Y, ssdMat, mapX2ssd, i2, j2,
                                   Q, q, use, REAL(tol)[0], &df, &a, &b, &n_co);
        idx[k-firstAdj] = UTE2I(i2, j2) + 1;

        if (n_I == 0) {
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
            p_values[k-firstAdj] = 2.0 * (1.0 - pt(fabs(lambda), n_co-q-2, 1, 0));
          if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
            statistic_values[k-firstAdj] = lambda;
            n_values[k-firstAdj] = n_co;
          }
        } else {
          if (!ISNAN(lambda)) {
            if (INTEGER(exactTest)[0]) {
              lambda = exp(lambda / ((double) -(use == USE_EM ? n : n_co)));
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[k-firstAdj] = pbeta(lambda, a, b, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[k-firstAdj] = lambda;
                n_values[k-firstAdj] = use == USE_EM ? n : n_co;
              }
            } else {
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
                p_values[k-firstAdj] = 1.0 - pchisq(lambda, df, TRUE, FALSE);
              if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
                statistic_values[k-firstAdj] = lambda;
                n_values[k-firstAdj] = use == USE_EM ? n : n_co;
              }
            }
          }
        }

        k++;
        if (startTime > 0 && k-firstAdj == nAdjEtime)
          break;
        if (verbose && startTime == 0) {
          pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
          if (pct != ppct) {
            SEXP s, t;
            PROTECT(t = s = allocList(3));
            SET_TYPEOF(s, LANGSXP);
            SETCAR(t, install("sendData")); t=CDR(t);
            SETCAR(t, masterNode);
            SET_TAG(t, install("node")); t=CDR(t);
            INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
            SETCAR(t, progressReport);
            SET_TAG(t, install("data"));
            eval(s, env);
            UNPROTECT(1); /* t s */
          }
          ppct = pct;
        }
      }
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      j_first = 0;
    }
  }

  if (k <= lastAdj && (startTime == 0 || k-firstAdj < nAdjEtime)) {
    int i_first = k - l_int * (l_ini + l_jni) - l_ini * l_jni;
    int l;

    /* intersection variables against themselves (avoiding pairing the same) */
    for (l = i_first; l < (l_int * (l_int - 1)) / 2 && k <= lastAdj; l++) {
      int i,j,i2,j2;
      i2e(l, &i, &j);

      i2 = pairup_ij_int[i] - 1;
      j2 = pairup_ij_int[j] - 1;

      if (n_I == 0) {
        if (work_with_margin) {
          ijQ[0] = i2;
          ijQ[1] = j2;
          memset(S, 0, sizeof(double) * n_upper_tri);
          n_co = ssd(REAL(XR), n_var, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
          lambda = qp_ci_test_std(S, q+2, n_co, 0, 1, Q, q, NULL);
        } else
          lambda = qp_ci_test_std(S, n_var, n, i2, j2, Q, q, NULL);
      } else
        lambda = qp_ci_test_hmgm(REAL(XR), n_var, n, I, n_I, INTEGER(n_levelsR),
                                 Y, n_Y, ssdMat, mapX2ssd, i2, j2,
                                 Q, q, use, REAL(tol)[0], &df, &a, &b, &n_co);
      idx[k-firstAdj] = UTE2I(i2, j2) + 1;

      if (n_I == 0) {
        if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
          p_values[k-firstAdj] = 2.0 * (1.0 - pt(fabs(lambda), n_co-q-2, 1, 0));
        if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
          statistic_values[k-firstAdj] = lambda;
          n_values[k-firstAdj] = n_co;
         }
      } else {
        if (!ISNAN(lambda)) {
          if (INTEGER(exactTest)[0]) {
            lambda = exp(lambda / ((double) -(use == USE_EM ? n : n_co)));
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
             p_values[k-firstAdj] = pbeta(lambda, a, b, TRUE, FALSE);
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
              statistic_values[k-firstAdj] = lambda;
              n_values[k-firstAdj] = use == USE_EM ? n : n_co;
            }
          } else {
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_PVALUE)
              p_values[k-firstAdj] = 1.0 - pchisq(lambda, df, TRUE, FALSE);
            if (return_type == RETURN_TYPE_ALL || return_type == RETURN_TYPE_STATN) {
              statistic_values[k-firstAdj] = lambda;
              n_values[k-firstAdj] = use == USE_EM ? n : n_co;
            }
          }
        }
      }

      k++;
      if (startTime > 0 && k-firstAdj == nAdjEtime)
        break;
      if (verbose && startTime == 0) {
        pct = (int) (((k-firstAdj) * 100) / n_adj_this_proc);
        if (pct != ppct) {
          SEXP s, t;
          PROTECT(t = s = allocList(3));
          SET_TYPEOF(s, LANGSXP);
          SETCAR(t, install("sendData")); t=CDR(t);
          SETCAR(t, masterNode);
          SET_TAG(t, install("node")); t=CDR(t);
          INTEGER(VECTOR_ELT(progressReport,1))[0] = k-firstAdj;
          SETCAR(t, progressReport);
          SET_TAG(t, install("data"));
          eval(s, env);
          UNPROTECT(1); /* t s */
        }
        ppct = pct;
      }
    }
  }

  Free(S); /* = Free(ssdMat) */

  if (n_I > 0) {
    if (!work_with_margin)
      Free(mapX2ssd);
    Free(Y);
    Free(I);
  }

  if (work_with_margin)
    Free(ijQ);

  if (QR != R_NilValue)
    Free(Q);

  if (startTime > 0) {
    SEXP procTimeR;
    double* procTime;
    SEXP call;

    PROTECT(call = lang1(install("proc.time")));
    PROTECT(procTimeR = eval(call, env));
    procTime = REAL(procTimeR);
    elapsedTime = elapsedTime + ((procTime[2] - startTime) / (double) (k-firstAdj)) * (double) n_adj_this_proc;
    UNPROTECT(2); /* call procTimeR */

    p_values[0] = elapsedTime; /* store in the first position of the p_values vector the estimated time */
  }

  UNPROTECT(4);   /* result result_names progressReport progressReport_names */

  return result;
}



/*
  FUNCTION: qp_fast_ci_test_std
  PURPOSE: wrapper of the R-C interface for calling the function that performs
           a test for conditional independence between variables i and j
           given de conditioning set C using standard calculations
  RETURNS: a list with two members, the t-statistic value and the p-value
           on rejecting the null hypothesis of independence
*/

static SEXP
qp_fast_ci_test_std(SEXP SR, SEXP pR, SEXP nR, SEXP iR, SEXP jR, SEXP QR) {
  int    n = INTEGER(nR)[0];
  int    p = INTEGER(pR)[0];
  int    q;
  int*   Q;
  int    i,j,k;
  double t_value;
  double df;
  double p_value;
  double beta;
  char   dataname[4096];
  SEXP   result;
  SEXP   result_names;
  SEXP   result_stat;
  SEXP   result_param;
  SEXP   result_p_val;
  SEXP   result_est;
  SEXP   class;
  SEXP   stat_name, param_name, pval_name, est_name, nullval_name;

  PROTECT_INDEX Spi,Qpi;

  PROTECT_WITH_INDEX(SR, &Spi);
  PROTECT_WITH_INDEX(QR, &Qpi);

  REPROTECT(SR = coerceVector(SR, REALSXP), Spi);
  REPROTECT(QR = coerceVector(QR, INTSXP), Qpi);

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;
  q = length(QR);

  sprintf(dataname, "%s and %s given {", CHAR(STRING_ELT(getAttrib(iR, R_NamesSymbol), 0)),
          CHAR(STRING_ELT(getAttrib(jR, R_NamesSymbol), 0)));
  Q = Calloc(q, int);
  for (k=0;k<q;k++) {
    Q[k] = INTEGER(QR)[k]-1;
    if (k > 0)
      strcat(dataname, ", ");
    strcat(dataname, CHAR(STRING_ELT(getAttrib(QR, R_NamesSymbol), k)));
  }
  strcat(dataname, "}");

  t_value = qp_ci_test_std(REAL(SR), p, n, i, j, Q, q, &beta);
  df = n - q - 2;
  p_value = 2.0 * (1.0 - pt(fabs(t_value), df, 1, 0));

  PROTECT(result = allocVector(VECSXP,8));
  SET_VECTOR_ELT(result,0,result_stat = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,1,result_param = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,2,result_p_val = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,3,result_est = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,4,allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,5,allocVector(STRSXP, 1));
  SET_VECTOR_ELT(result,6,allocVector(STRSXP, 1));
  SET_VECTOR_ELT(result,7,allocVector(STRSXP, 1));
  PROTECT(result_names = allocVector(STRSXP,8));
  SET_STRING_ELT(result_names,0,mkChar("statistic"));
  SET_STRING_ELT(result_names,1,mkChar("parameter"));
  SET_STRING_ELT(result_names,2,mkChar("p.value"));
  SET_STRING_ELT(result_names,3,mkChar("estimate"));
  SET_STRING_ELT(result_names,4,mkChar("null.value"));
  SET_STRING_ELT(result_names,5,mkChar("alternative"));
  SET_STRING_ELT(result_names,6,mkChar("method"));
  SET_STRING_ELT(result_names,7,mkChar("data.name"));
  setAttrib(result,R_NamesSymbol,result_names);

  PROTECT(stat_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,0))[0] = t_value;
  SET_STRING_ELT(stat_name,0,mkChar("t"));
  setAttrib(VECTOR_ELT(result,0),R_NamesSymbol,stat_name);

  PROTECT(param_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,1))[0] = df;
  SET_STRING_ELT(param_name,0,mkChar("df"));
  setAttrib(VECTOR_ELT(result,1),R_NamesSymbol,param_name);

  PROTECT(pval_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,2))[0] = p_value;
  SET_STRING_ELT(pval_name,0,mkChar("two.sided"));
  setAttrib(VECTOR_ELT(result,2),R_NamesSymbol,pval_name);

  PROTECT(est_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,3))[0] = beta;
  SET_STRING_ELT(est_name,0,mkChar("beta"));
  setAttrib(VECTOR_ELT(result,3),R_NamesSymbol,est_name);

  PROTECT(nullval_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,4))[0] = 0;
  SET_STRING_ELT(nullval_name,0,mkChar("partial regression coefficient"));
  setAttrib(VECTOR_ELT(result,4),R_NamesSymbol,nullval_name);

  SET_STRING_ELT(VECTOR_ELT(result,5), 0, mkChar("two.sided"));
  SET_STRING_ELT(VECTOR_ELT(result,6), 0, mkChar("Conditional independence test for continuous data using a t test for zero partial regression coefficient"));
  SET_STRING_ELT(VECTOR_ELT(result,7), 0, mkChar(dataname));

  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("htest"));
  installAttrib(result, R_ClassSymbol, class);

  UNPROTECT(10); /* S QR result result_names stat_name param_name pval_name est_name nullval_name class */

  Free(Q);

  return result;
}



/*
  FUNCTION: qp_fast_ci_test_opt
  PURPOSE: wrapper of the R-C interface for calling the function that performs
           a test for conditional independence between variables i and j
           given de conditioning set C using optimized calculations
  RETURNS: a list with two members, the t-statistic value and the p-value
           on rejecting the null hypothesis of independence
*/

static SEXP
qp_fast_ci_test_opt(SEXP SR, SEXP pR, SEXP nR, SEXP iR, SEXP jR, SEXP QR) {
  int    n = INTEGER(nR)[0];
  int    p = INTEGER(pR)[0];
  int    q;
  int*   Q;
  int    i,j,k;
  double t_value;
  double df;
  double p_value;
  double beta;
  char   dataname[4096];
  SEXP   result;
  SEXP   result_names;
  SEXP   result_stat;
  SEXP   result_param;
  SEXP   result_p_val;
  SEXP   result_est;
  SEXP   class;
  SEXP   stat_name, param_name, pval_name, est_name, nullval_name;

  PROTECT_INDEX Spi,Qpi;

  PROTECT_WITH_INDEX(SR, &Spi);
  PROTECT_WITH_INDEX(QR, &Qpi);

  REPROTECT(SR = coerceVector(SR, REALSXP), Spi);
  REPROTECT(QR = coerceVector(QR, INTSXP), Qpi);

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;
  q = length(QR);

  sprintf(dataname, "%s and %s given {", CHAR(STRING_ELT(getAttrib(iR, R_NamesSymbol), 0)),
          CHAR(STRING_ELT(getAttrib(jR, R_NamesSymbol), 0)));
  Q = Calloc(q, int);
  for (k=0;k<q;k++) {
    Q[k] = INTEGER(QR)[k]-1;
    if (k > 0)
      strcat(dataname, ", ");
    strcat(dataname, CHAR(STRING_ELT(getAttrib(QR, R_NamesSymbol), k)));
  }
  strcat(dataname, "}");

  t_value = qp_ci_test_opt(REAL(SR), p, n, i, j, Q, q, NULL, &beta);
  df = n - q - 2;
  p_value = 2.0 * (1.0 - pt(fabs(t_value), df, 1, 0));

  PROTECT(result = allocVector(VECSXP,8));
  SET_VECTOR_ELT(result,0,result_stat = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,1,result_param = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,2,result_p_val = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,3,result_est = allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,4,allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,5,allocVector(STRSXP, 1));
  SET_VECTOR_ELT(result,6,allocVector(STRSXP, 1));
  SET_VECTOR_ELT(result,7,allocVector(STRSXP, 1));
  PROTECT(result_names = allocVector(STRSXP,8));
  SET_STRING_ELT(result_names,0,mkChar("statistic"));
  SET_STRING_ELT(result_names,1,mkChar("parameter"));
  SET_STRING_ELT(result_names,2,mkChar("p.value"));
  SET_STRING_ELT(result_names,3,mkChar("estimate"));
  SET_STRING_ELT(result_names,4,mkChar("null.value"));
  SET_STRING_ELT(result_names,5,mkChar("alternative"));
  SET_STRING_ELT(result_names,6,mkChar("method"));
  SET_STRING_ELT(result_names,7,mkChar("data.name"));
  setAttrib(result,R_NamesSymbol,result_names);

  PROTECT(stat_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,0))[0] = t_value;
  SET_STRING_ELT(stat_name,0,mkChar("t"));
  setAttrib(VECTOR_ELT(result,0),R_NamesSymbol,stat_name);

  PROTECT(param_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,1))[0] = df;
  SET_STRING_ELT(param_name,0,mkChar("df"));
  setAttrib(VECTOR_ELT(result,1),R_NamesSymbol,param_name);

  PROTECT(pval_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,2))[0] = p_value;
  SET_STRING_ELT(pval_name,0,mkChar("two.sided"));
  setAttrib(VECTOR_ELT(result,2),R_NamesSymbol,pval_name);

  PROTECT(est_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,3))[0] = beta;
  SET_STRING_ELT(est_name,0,mkChar("beta"));
  setAttrib(VECTOR_ELT(result,3),R_NamesSymbol,est_name);

  PROTECT(nullval_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,4))[0] = 0;
  SET_STRING_ELT(nullval_name,0,mkChar("partial regression coefficient"));
  setAttrib(VECTOR_ELT(result,4),R_NamesSymbol,nullval_name);

  SET_STRING_ELT(VECTOR_ELT(result,5), 0, mkChar("two.sided"));
  SET_STRING_ELT(VECTOR_ELT(result,6), 0, mkChar("Conditional independence test for continuous data using a t test for zero partial regression coefficient"));
  SET_STRING_ELT(VECTOR_ELT(result,7), 0, mkChar(dataname));

  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("htest"));
  installAttrib(result, R_ClassSymbol, class);

  UNPROTECT(10); /* S QR result result_names stat_name param_name pval_name est_name nullval_name class */

  Free(Q);

  return result;
}



/*
  FUNCTION: qp_ci_test_std
  PURPOSE: perform a test for conditional independence between variables
           indexed by i and j given the conditioning set Q using standard calculations
  RETURNS: the t-statistic value on rejecting the null hypothesis of conditional independence
*/

static double
qp_ci_test_std(double* S, int n_var, int n, int i, int j, int* Q, int q, double* beta) {
  int*    subvars;
  int     subn = q + 2;
  int     k,l;
  double* Mmar;
  double  S11;
  double* S12;
  double* S21;
  double* S22;
  double* S22inv;
  double* S22inv1col;
  double* tmpmat;
  double  tmpval;
  double  betahat;
  double  sigma;
  double  se;
  double  t_value;

  subvars     = Calloc(subn,int);
  Mmar        = Calloc(subn*subn,double);
  S12         = Calloc(subn,double);
  S21         = Calloc(subn,double);
  S22         = Calloc((subn-1)*(subn-1),double);
  S22inv      = Calloc((subn-1)*(subn-1),double);
  S22inv1col  = Calloc(subn-1,double);

  subvars[0] = i; /* order here is important, first variable i */
  subvars[1] = j; /* then variable j then the conditioning set */
  for (k=2;k<subn;k++)
    subvars[k] = Q[k-2];

  /* Mmar <- S[c(i, j, sp), c(i, j, sp)] 
     S11     <- Mmar[1,1]
     S12     <- Mmar[1,-1]
     S21     <- Mmar[-1,1]
     S22     <- Mmar[-1,-1] */
/*
  for (k=0;k<subn;k++)
    for (l=0;l<subn;l++) {
      Mmar[k+l*subn] = S[subvars[k]+subvars[l]*n_var];
      if (k == 0 && l > 0)
        S12[l-1] = Mmar[k+l*subn];
      if (k > 0 && l == 0)
        S21[k-1] = Mmar[k+l*subn];
      if (k > 0 && l > 0)
        S22[k-1+(l-1)*(subn-1)] = Mmar[k+l*subn];
    }
  S11 = Mmar[0];
*/
  for (k=0;k<subn;k++)
    for (l=0;l<subn;l++) {
      Mmar[k+l*subn] = S[UTE2I(subvars[k],subvars[l])]; /* S is a vector storing the upper */
      if (k == 0 && l > 0)                              /* triangle of the sample covariance */
        S12[l-1] = Mmar[k+l*subn];                      /* matrix in column-major order */
      if (k > 0 && l == 0)
        S21[k-1] = Mmar[k+l*subn];
      if (k > 0 && l > 0)
        S22[k-1+(l-1)*(subn-1)] = Mmar[k+l*subn];
    }
  S11 = Mmar[0];

  /* S22inv  <- solve(S22) */
  matinv(S22inv,S22,subn-1, 0);

  /* betahat <- S12 %*% S22inv[,1] */
  Memcpy(S22inv1col,S22inv,(size_t) (subn-1));
  matprod(S12,1,subn-1,S22inv1col,subn-1,1,&betahat);

  /* sigma   <- sqrt((S11 - S12 %*% S22inv %*% S21) * (n - 1) / (n - q - 2)) */
  tmpmat = Calloc(subn-1,double);
  matprod(S22inv,subn-1,subn-1,S21,subn-1,1,tmpmat);
  matprod(S12,1,subn-1,tmpmat,subn-1,1,&tmpval);
  Free(tmpmat);
  sigma = sqrt( (S11 - tmpval) * (n - 1) / (n - subn) );
  /* se      <- sigma * sqrt(S22inv[1,1] / (n - 1)) */
  se = sigma * sqrt(S22inv[0] / (n - 1));
  /* t.value <- betahat / se */
  t_value = betahat / se;

  if (beta != NULL)
    *beta = betahat;

  Free(S22inv1col);
  Free(S22inv);
  Free(S22);
  Free(S21);
  Free(S12);
  Free(Mmar);
  Free(subvars);

  return t_value;
}



/*
  FUNCTION: qp_ci_test_opt
  PURPOSE: perform a test for conditional independence between variables
           indexed by i and j given the conditioning set Q using optimized
           calculations that allow one to use inverted sample covariance
           matrices of the conditioning sets pre-calculated beforehand, thus
           enabling using a common subset of conditioning sets for all pairs
           of variables which significantly decreases the overall computational
           cost of estimating non-rejection rates for a large number of
           pairs of variables
  RETURNS: the t-statistic value of the test on rejecting the null hypothesis
           of conditional independence
*/

static double
qp_ci_test_opt(double* S, int n_var, int N, int i, int j, int* Q, int q,
               double* Qinv, double* beta) {
  int*    subvars;
  int     subn = q + 2;
  int     k,l;
  double* Qmat;
  double* Sij;
  double* Sijbyq;
  double* Sqbyij;
  double* tmpmat1;
  double* tmpmat2;
  double* par_cov;
  double* par_cor;
  double  betahat;
  double  se;
  double  t_value;
  int     flagNoQinv=0;

  subvars     = Calloc(subn,int);
  Sij         = Calloc(4,double);
  Sijbyq      = Calloc(2*q, double);
  Sqbyij      = Calloc(q*2, double);
  par_cov     = Calloc(4,double);
  par_cor     = Calloc(4,double);

  subvars[0] = i; /* order here is important, first variable i */
  subvars[1] = j; /* then variable j then the conditioning set */
  for (k=2;k<subn;k++)
    subvars[k] = Q[k-2];
/*
  for (k=0;k<subn;k++)
    for (l=0;l<subn;l++) {
      if (k < 2 && l < 2)
        Sij[k+l*2] = S[subvars[k]+subvars[l]*n_var];
      if (k < 2 && l > 1) {
        Sijbyq[k+(l-2)*2] = S[subvars[k]+subvars[l]*n_var];
        Sqbyij[l-2+k*q] = S[subvars[l]+subvars[k]*n_var];
      }
    }
*/
  for (k=0;k<subn;k++)
    for (l=0;l<subn;l++) {
      double x = S[UTE2I(subvars[k], subvars[l])];

      if (k < 2 && l < 2)
        Sij[k+l*2] = x;
      if (k < 2 && l > 1) {
        Sijbyq[k+(l-2)*2] = x;
        Sqbyij[l-2+k*q] = x;
      }
    }

  if (Qinv == NULL) {
    Qmat = Calloc(q*q, double);
    Qinv = Calloc(q*q, double);
    /*
    for (i=0; i < q; i++) {
      for (j=0; j < i; j++)
        Qmat[i + j*q] = Qmat[j + i*q] = S[Q[i] + Q[j] * n_var];
      Qmat[i + i*q] = S[Q[i] + Q[i] * n_var];
    }
    */
    for (i=0; i < q; i++) {
      for (j=0; j < i; j++)
        Qmat[i + j*q] = Qmat[j + i*q] = S[UTE2I(Q[i], Q[j])];
      Qmat[i + i*q] = S[UTE2I(Q[i], Q[i])];
    }
    if (q > 1)
      matinv(Qinv,Qmat,q, 0);
    else
      Qinv[0] = 1.0 / Qmat[0];
    Free(Qmat);
    flagNoQinv=1;
  }

  tmpmat1 = Calloc(q*2,double);
  tmpmat2 = Calloc(4,double);
  matprod(Qinv,q,q,Sqbyij,q,2,tmpmat1);
  matprod(Sijbyq,2,q,tmpmat1,q,2,tmpmat2);
  Free(tmpmat1);
  matsumf(par_cov, 2, 2, Sij, tmpmat2, -1.0);
  Free(tmpmat2);
  Free(Sij);
  Free(Sijbyq);
  Free(Sqbyij);
  cov2cor(par_cor, par_cov, 2);
  Free(par_cov);
  betahat = sqrt(N - q -2) * par_cor[2];
  se = sqrt(1.0 - par_cor[2] * par_cor[2]);

  if (beta != NULL)
    *beta = betahat;

  /* t.value <- betahat / se */
  t_value = betahat / se;

  Free(par_cor);
  Free(subvars);

  if (flagNoQinv)
    Free(Qinv);

  return t_value;
}



/*
  FUNCTION: qp_fast_ci_test_hmgm
  PURPOSE: wrapper of the R-C interface for calling the function that performs
           a test for conditional independence between variables i and j
           given de conditioning set Q where variables indicated by i and Q can
           be discrete
  RETURNS: a list with two members, the likelihood ratio statistic value and the
           p-value on rejecting the null hypothesis of independence
*/

static SEXP
qp_fast_ci_test_hmgm(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP ssdR,
                     SEXP mapX2ssdR, SEXP iR, SEXP jR, SEXP QR, SEXP exactTestR,
                     SEXP use, SEXP tol) {
  int     n = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  int     p = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int     q = length(QR);
  int     i,j,k, n_co;
  int     exactTest = INTEGER(exactTestR)[0];
  int*    I;
  int*    Y;
  int*    Q;
  double* ssd = NULL;
  int*    mapX2ssd = NULL;
  double  lambda;
  double  df, a, b;
  double  p_value = R_NaReal;
  char    dataname[4096];
  SEXP    result;
  SEXP    result_names;
  SEXP    result_stat;
  SEXP    result_param;
  SEXP    result_p_val;
  SEXP    class;
  SEXP    stat_name, param_name, pval_name, nullval_name;

  PROTECT_INDEX ssd_pi;

  if (ssdR != R_NilValue) {
    PROTECT_WITH_INDEX(ssdR,&ssd_pi);

    REPROTECT(ssdR = coerceVector(ssdR,REALSXP), ssd_pi);
  }

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;

  I = Calloc(n_I, int);
  for (k=0; k < n_I; k++)
    I[k] = INTEGER(IR)[k]-1;

  Y = Calloc(n_Y, int);
  for (k=0; k < n_Y; k++)
    Y[k] = INTEGER(YR)[k]-1;

  sprintf(dataname, "%s and %s given {", CHAR(STRING_ELT(getAttrib(iR, R_NamesSymbol), 0)),
          CHAR(STRING_ELT(getAttrib(jR, R_NamesSymbol), 0)));
  Q = Calloc(q, int);
  for (k=0; k < q; k++) {
    Q[k] = INTEGER(QR)[k]-1;
    if (k > 0)
      strcat(dataname, ", ");
    strcat(dataname, CHAR(STRING_ELT(getAttrib(QR, R_NamesSymbol), k)));
  }
  strcat(dataname, "}");

  if (ssdR != R_NilValue) {
    mapX2ssd = Calloc(p, int);
    for (k=0; k < p; k++)
      mapX2ssd[k] = INTEGER(mapX2ssdR)[k]-1;
    ssd = REAL(ssdR);
  }

  lambda = qp_ci_test_hmgm(REAL(XR), p, n, I, n_I, INTEGER(n_levelsR), Y, n_Y,
                           ssd, mapX2ssd, i, j, Q, q, INTEGER(use)[0], REAL(tol)[0],
                           &df, &a, &b, &n_co);

  if (!ISNAN(lambda)) {
    if (exactTest) {
      lambda = exp(lambda / ((double) -(INTEGER(use)[0] == USE_EM ? n : n_co)));
      p_value = pbeta(lambda, a, b, TRUE, FALSE);
    } else
      p_value = 1.0 - pchisq(lambda, df, TRUE, FALSE);
  }

  PROTECT(result = allocVector(VECSXP,8));
  SET_VECTOR_ELT(result,0,result_stat = allocVector(REALSXP, 1));
  SET_VECTOR_ELT(result,1,result_param = allocVector(REALSXP, exactTest ? 3 : 2));
  SET_VECTOR_ELT(result,2,result_p_val = allocVector(REALSXP, 1));
  SET_VECTOR_ELT(result,3,R_NilValue); /* no estimated parameter, could be h = mu x Sigma^{-1} */
  SET_VECTOR_ELT(result,4,allocVector(REALSXP,1));
  SET_VECTOR_ELT(result,5,allocVector(STRSXP, 1));
  SET_VECTOR_ELT(result,6,allocVector(STRSXP, 1));
  SET_VECTOR_ELT(result,7,allocVector(STRSXP, 1));
  PROTECT(result_names = allocVector(STRSXP,8));
  SET_STRING_ELT(result_names,0,mkChar("statistic"));
  SET_STRING_ELT(result_names,1,mkChar("parameter"));
  SET_STRING_ELT(result_names,2,mkChar("p.value"));
  SET_STRING_ELT(result_names,3,mkChar("estimate"));
  SET_STRING_ELT(result_names,4,mkChar("null.value"));
  SET_STRING_ELT(result_names,5,mkChar("alternative"));
  SET_STRING_ELT(result_names,6,mkChar("method"));
  SET_STRING_ELT(result_names,7,mkChar("data.name"));
  setAttrib(result,R_NamesSymbol,result_names);

  PROTECT(stat_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,0))[0] = lambda;
  SET_STRING_ELT(stat_name,0,exactTest ? mkChar("Lambda") : mkChar("-n log Lambda"));
  setAttrib(VECTOR_ELT(result,0),R_NamesSymbol,stat_name);

  PROTECT(param_name = allocVector(STRSXP, exactTest ? 3 : 2));
  if (exactTest) {
    REAL(VECTOR_ELT(result,1))[0] = a;
    REAL(VECTOR_ELT(result,1))[1] = b;
    REAL(VECTOR_ELT(result,1))[2] = INTEGER(use)[0] == USE_EM ? n : n_co;
    SET_STRING_ELT(param_name,0,mkChar("a"));
    SET_STRING_ELT(param_name,1,mkChar("b"));
    SET_STRING_ELT(param_name,2,mkChar("n"));
    setAttrib(VECTOR_ELT(result,1),R_NamesSymbol,param_name);
  } else {
    REAL(VECTOR_ELT(result,1))[0] = df;
    REAL(VECTOR_ELT(result,1))[1] = INTEGER(use)[0] == USE_EM ? n : n_co;
    SET_STRING_ELT(param_name,0,mkChar("df"));
    SET_STRING_ELT(param_name,1,mkChar("n"));
    setAttrib(VECTOR_ELT(result,1),R_NamesSymbol,param_name);
  }

  PROTECT(pval_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,2))[0] = p_value;
  SET_STRING_ELT(pval_name,0,exactTest ? mkChar("less") : mkChar("greater"));
  setAttrib(VECTOR_ELT(result,2),R_NamesSymbol,pval_name);

  PROTECT(nullval_name = allocVector(STRSXP,1));
  REAL(VECTOR_ELT(result,4))[0] = exactTest ? 1 : 0;
  SET_STRING_ELT(nullval_name,0,exactTest ? mkChar("Lambda") : mkChar("-n log Lambda"));
  setAttrib(VECTOR_ELT(result,4),R_NamesSymbol,nullval_name);

  SET_STRING_ELT(VECTOR_ELT(result,5), 0, exactTest ? mkChar("less") : mkChar("greater"));
  SET_STRING_ELT(VECTOR_ELT(result,6), 0, exactTest ?
                                          mkChar("Conditional independence test for homogeneous mixed data using an exact likelihood ratio test") :
                                          mkChar("Conditional independence test for homogeneous mixed data using an asymptotic likelihood ratio test"));
  SET_STRING_ELT(VECTOR_ELT(result,7), 0, mkChar(dataname));

  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(class, 0, mkChar("htest"));
  installAttrib(result, R_ClassSymbol, class);

  UNPROTECT(7); /* result result_names stat_name param_name pval_name nullval_name class */

  Free(I);
  Free(Y);
  Free(Q);

  if (ssdR != R_NilValue) {
    UNPROTECT(1); /* ssdR */
    Free(mapX2ssd);
  }

  return result;
}



/*
  FUNCTION: qp_ci_test_hmgm
  PURPOSE: perform a test for conditional independence between variables
           indexed by i and j given the conditioning set Q.
  RETURNS: a list with two members, the t-statistic value and the p-value
           on rejecting the null hypothesis of independence
*/

static double
lr_complete_obs(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y,
                int n_Y, double* ucond_ssd, int* mapX2ucond_ssd, int total_n_Y,
                int i, int j, int* Q, int q, int* n_co) {
  int     k;
  int     n_I_i, n_Y_i, n_Y_j, n_Y_ij;
  int     sign;
  int     final_sign = 1;
  int     flag_zero = FALSE;
  int*    idx_misobs=NULL;
  double* ssd_mat;
  double  x;
  double  lr = 0.0;

  *n_co = n;

  /* in order to save memory we will use ssd_mat for calculating each ssd matrix */
  /* the call below to ssd_A() assumes that this Calloc() memsets ssd_mat to zeroes */
  ssd_mat = Calloc((n_Y * (n_Y + 1)) / 2, double);  /* upper triangle includes the diagonal */

  if (n_I > 0 || ucond_ssd == NULL) {
    idx_misobs = Calloc(n, int); /* assume Calloc() memsets everything to zeroes */

    /* when calculating the larger ssd, the number of complete observations and the
     * logical mask of missing observations is pulled out and employed in later calls
     * to ssd_A */
    *n_co = ssd_A(X, p, n, I, n_I, n_levels, Y, n_Y, NULL, idx_misobs, ssd_mat);
  } else {
    int* tmp;

    tmp = Calloc(n_Y, int);
    for (k=0; k < n_Y; k++)
      tmp[k] = mapX2ucond_ssd[Y[k]];

    symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y);

    Free(tmp);
  }
/*
  Rprintf("ssd (n_Y=%d):\n", n_Y);
  int m = 0;
  for (k=0; k < n_Y; k++) {
    int l;
    for (l=0; l <= k; l++) {
       Rprintf("%10.6f\t", ssd_mat[m++]);
    }
    Rprintf("\n");
  }
  Rprintf("n=%d\n", *n_co);
*/
  lr = x = symmatlogdet(ssd_mat, n_Y, &sign);
  if (x < -DBL_DIG)
    flag_zero = TRUE;
  final_sign *= sign;

/*
  Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y, &sign));
  Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/

  /* ssd_i = ssd_Gamma when i is discrete or ssd_{Gamma\i} when i is continuous */

  k = 0;
  while (i != I[k] && k < n_I)
    k++;

  n_I_i = n_I;
  n_Y_i = n_Y;
  if (k < n_I) {  /* i is discrete */
    int tmp;

    I[k] = I[n_I-1];
    I[n_I-1] = i;
    n_I_i = n_I - 1;
    tmp = n_levels[k];
    n_levels[k] = n_levels[n_I-1];
    n_levels[n_I-1] = tmp;
  } else {       /* i is continuous */
    k = 0;
    while (i != Y[k] && k < n_Y)
      k++;

    if (k < n_Y) {
      Y[k] = Y[n_Y-1];
      Y[n_Y-1] = i;
      n_Y_i = n_Y - 1;
    } else
      error("qp_ci_test_hmgm(): i does not form part of neither I nor Y\n");
  }

  if (n_I > 0 || ucond_ssd == NULL) {
    memset(ssd_mat, 0, ((n_Y_i * (n_Y_i + 1)) / 2) * sizeof(double));
    ssd_A(X, p, n, I, n_I_i, n_levels, Y, n_Y_i, idx_misobs, NULL, ssd_mat);
  } else {
    int* tmp;

    tmp = Calloc(n_Y_i, int);
    for (k=0; k < n_Y_i; k++)
      tmp[k] = mapX2ucond_ssd[Y[k]];

    symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y_i);

    Free(tmp);
  }

/*
  Rprintf("ssd_i(n_Y_i=%d):\n", n_Y_i);
  m = 0;
  for (k=0; k < n_Y_i; k++) {
    int l;
    Rprintf("%d", Y[k]+1);
    for (l=0; l <= k; l++) {
       Rprintf("\t%10.6f", ssd_mat[m++]);
    }
    Rprintf("\n");
  }
*/

  x = symmatlogdet(ssd_mat, n_Y_i, &sign);
  lr -= x;
  if (x < -DBL_DIG)
    flag_zero = TRUE;
  final_sign *= sign;

/*
  Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y_i, &sign));
  Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/

  if (n_Y > 1) {
    n_Y_j = n_Y;
    k = 0;
    while (j != Y[k] && k < n_Y)
      k++;

    if (k < n_Y) {
      Y[k] = Y[n_Y-1];
      Y[n_Y-1] = j;
      n_Y_j = n_Y - 1;
    }

    if (n_I > 0 || ucond_ssd == NULL) {
      memset(ssd_mat, 0, ((n_Y_j * (n_Y_j + 1)) / 2) * sizeof(double));
      ssd_A(X, p, n, I, n_I, n_levels, Y, n_Y_j, idx_misobs, NULL, ssd_mat);
    } else {
      int* tmp;

      tmp = Calloc(n_Y_j, int);
      for (k=0; k < n_Y_j; k++)
        tmp[k] = mapX2ucond_ssd[Y[k]];

      symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y_j);

      Free(tmp);
    }

/*
    Rprintf("ssd_j:\n");
    m = 0;
    for (k=0; k < n_Y_j; k++) {
      int l;
      Rprintf("%d", Y[k]+1);
      for (l=0; l <= k; l++) {
         Rprintf("\t%10.6f", ssd_mat[m++]);
      }
      Rprintf("\n");
    }
*/

    x = symmatlogdet(ssd_mat, n_Y_j, &sign);
    lr -= x;
    if (x < -DBL_DIG)
      flag_zero = TRUE;
    final_sign *= sign;

/* 
    Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y_j, &sign));
    Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/

    n_Y_ij = n_Y_j;
    k = 0;
    while (i != Y[k] && k < n_Y)
      k++;

    if (k < n_Y) {
      Y[k] = Y[n_Y-2];
      Y[n_Y-2] = i;
      n_Y_ij = n_Y_j - 1;
    }

    if (n_Y_ij > 0) {
      if (n_I > 0 || ucond_ssd == NULL) {
        memset(ssd_mat, 0, ((n_Y_j * (n_Y_j + 1)) / 2) * sizeof(double));
        ssd_A(X, p, n, I, n_I_i, n_levels, Y, n_Y_ij, idx_misobs, NULL, ssd_mat);
      } else {
        int* tmp;

        tmp = Calloc(n_Y_ij, int);
        for (k=0; k < n_Y_ij; k++)
          tmp[k] = mapX2ucond_ssd[Y[k]];

        symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y_ij);

        Free(tmp);
      }

/*
      Rprintf("ssd_ij:\n");
      m = 0;
      for (k=0; k < n_Y_ij; k++) {
        int l;
        Rprintf("%d", Y[k]+1);
        for (l=0; l <= k; l++) {
           Rprintf("\t%10.6f", ssd_mat[m++]);
        }
        Rprintf("\n");
      }
*/

      x = symmatlogdet(ssd_mat, n_Y_ij, &sign);
      lr += x;
      if (x < -DBL_DIG)
        flag_zero = TRUE;
      final_sign *= sign;

/*
      Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y_ij, &sign));
      Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/
    }

  }

  Free(ssd_mat);

  if (idx_misobs != NULL)
    Free(idx_misobs);

  if (flag_zero || final_sign == -1)
    lr = R_NaN;
  else
    lr = lr * ((double) -(*n_co));

  return lr;
}

/* complete sufficient statistics: calculate sufficient statistics (En, Es and Ess)
 * of those observations that do not have missing values */

com_stats_t
new_com_stats(int n_joint_levels, int n_Y) {
  com_stats_t cs = { NULL, NULL, NULL};

  if (n_joint_levels > 0 && n_Y > 0)
    cs.Es_com = Calloc(n_joint_levels * n_Y, double);

  if (n_Y > 0)
    cs.Ess_com = Calloc(n_Y * n_Y, double);

  if (n_joint_levels > 0)
    cs.n_com = Calloc(n_joint_levels, int); /* assuming this memsets to zeroes */

  return cs;
}

void
free_com_stats(com_stats_t cs) {
  if (cs.Es_com != NULL)
    Free(cs.Es_com);

  if (cs.Ess_com != NULL)
    Free(cs.Ess_com);

  if (cs.n_com != NULL)
    Free(cs.n_com);
}

com_stats_t
stat_com(double* X, int p, int n, int* missing_mask, int n_mis, int* Is, int n_Is,
         int *Y, int n_Y, int* n_levels, int n_joint_levels) {
  int         i,j,k,l,m;
  int         n_obs=0;
  com_stats_t cs;

  cs = new_com_stats(n_joint_levels, n_Y);

  if (n-n_mis > 0 && n_Y > 0) {
    int* obs_idx;

    obs_idx = Calloc(n, int);
    global_xtab = Calloc(n, int);
    for (i=0; i < n; i++) {
      global_xtab[i] = 1;
      if (missing_mask[i])
        global_xtab[i] = -1;
      else
        obs_idx[n_obs++] = i;
    }
    calculate_xtab(X, p, n, Is, n_Is, n_levels, global_xtab);

    /* group together observations from the same joint discrete level putting the *
     * observations from missing discrete levels at the beginning                 */
    qsort(obs_idx, n_obs, sizeof(int), indirect_int_cmp);

    /* skip missing discrete observations */
    i = 0;
    while (i < n_obs && global_xtab[obs_idx[i]] < 1)
      i++;

    l=0; /* counts distinct levels with complete observations */
    /* sum continuous variables through complete obs per level. figure out how to
     * assign 0 to the unobserved levels and whether that's necessary at all */
    while (i < n_obs) {
      j = i;
      while (j < n_obs && global_xtab[obs_idx[i]] == global_xtab[obs_idx[j]]) {
        for (k=0; k < n_Y; k++) {
          cs.Es_com[k*n_joint_levels+l] += X[Y[k]*n+obs_idx[j]];
          for (m=0; m < n_Y; m++) {
            cs.Ess_com[k*n_Y+m] += X[Y[k]*n+obs_idx[j]]*X[Y[m]*n+obs_idx[j]];
          }
        }
        j++;
      }

      cs.n_com[l] = j-i;
      l++;
      i = j;
    }

    Free(global_xtab);
    Free(obs_idx);
  }

  return cs;
}

/* sufficient statistics from missing data: perform the E step for the observations
 * with missing values */

 /* sufstat_i <- stat_mis(X, I_i, Y_i, levels_I_i, comStat_i, I, levels_I, mapX2I, Y, mapX2Y, p, mu, Sigma) */

suf_stats_t
new_suf_stats(int n_joint_levels, int n_Y) {
  suf_stats_t ss = { NULL, NULL, NULL, NULL};

  if (n_joint_levels > 0 && n_Y > 0) {
    ss.h = Calloc(n_joint_levels * n_Y, double);
    ss.bar_y = Calloc(n_joint_levels * n_Y, double);
  }

  if (n_Y > 0) {
    ss.ssd = Calloc(n_Y * n_Y, double);
    ss.K = Calloc(n_Y * n_Y, double);
  }

  if (n_joint_levels > 0)
    ss.m = Calloc(n_joint_levels, double);

  return ss;
}

void
free_suf_stats(suf_stats_t ss) {
  if (ss.h != NULL)
    Free(ss.h);

  if (ss.bar_y != NULL)
    Free(ss.bar_y);

  if (ss.ssd != NULL)
    Free(ss.ssd);

  if (ss.K != NULL)
    Free(ss.K);

  if (ss.m != NULL)
    Free(ss.m);
}

/* k(i) = y^T\Sigma^{-1}\mu(i) - 1/2 * [y^T\Sigma^{-1} y + \mu(i)^T \Sigma^{-1}\mu(i)] + log p(i) */
double
Ki(double* X, int p, int n, int i, int* Y, int n_Y, int j, double* Sigma, double* mu, int n_joint_levels, double* pr) {
  int     k,l;
  double  aux;
  double* xYs;
  double* t_xYs;
  double* sigmaYs;
  double* inv_sigmaYs;
  double* muYs;
  double* t_muYs;
  double* tmp;
  double  t1,t2,t3;

  aux = pr[j];
  if (n_Y > 0) {
    xYs = Calloc(n_Y, double);
    t_xYs = Calloc(n_Y, double);
    sigmaYs = Calloc(n_Y*n_Y, double);
    inv_sigmaYs = Calloc(n_Y*n_Y, double);
    muYs = Calloc(n_Y, double);
    t_muYs = Calloc(n_Y, double);
    tmp = Calloc(n_Y*n_Y, double);

    for (k=0; k < n_Y; k++)
      xYs[k] = X[n * Y[k] + i];
    mattran(t_xYs, xYs, 1, n_Y);
    for (k=0; k < n_Y; k++)
      for (l=0; l < n_Y; l++)
        sigmaYs[n_Y * k + l] = Sigma[n * Y[k] + l];
    matinv(inv_sigmaYs, sigmaYs, n_Y, n_Y);
    for (k=0; k < n_Y; k++)
      muYs[k] = mu[n_joint_levels * k + j];
    mattran(t_muYs, muYs, 1, n_Y);

    matprod(t_xYs, 1, n_Y, inv_sigmaYs, n_Y, n_Y, tmp);
    matprod(tmp, 1, n_Y, muYs, n_Y, 1, &t1);

    matprod(tmp, 1, n_Y, xYs, n_Y, 1, &t2);

    matprod(t_muYs, 1, n_Y, inv_sigmaYs, n_Y, n_Y, tmp);
    matprod(tmp, 1, n_Y, muYs, n_Y, 1, &t3);

    aux = exp(t1 - 0.5 * (t2 + t3) + log(pr[j]));

    Free(tmp);
    Free(t_muYs);
    Free(muYs);
    Free(inv_sigmaYs);
    Free(sigmaYs);
    Free(t_xYs);
    Free(xYs);

  }

  return aux;
}

/* calculate the probability of I=i for each observation with missing values
 * pr(I=i' | (i_{obs}, y}^{(\nu)}) exp k(i') / \sum_{s\in{\cal S}} exp k(s) */
double*
prob_i(double* X, int p, int n, int * missing_mask, int n_mis, int* I, int n_I,
       int* Is, int n_Is, int *Y, int n_Y, int* n_levels, int n_joint_levels,
       int k, double* pr, double* mu, double* Sigma) {
  int     i,j;
  int*    idx_I_obs;
  int     n_idx_I_obs=0;
  int*    index_S;
  int*    index_Si;
  int     n_index_S, n_index_Si;
  double  K_den;
  double* pr_i;

  pr_i = Calloc(n_mis, double); /* must be freed outside */
  idx_I_obs = Calloc(n_I, int);

  for (i=0; i < n; i++) {
    if (missing_mask[i]) {
      int base = 1;
      int joint_level_Is_obs=1; /* SHOULD THIS BE ZERO ?? */ 
      int joint_level_k=1; /* SHOULD THIS BE ZERO ?? */

      for (j=0; j < n_Is; j++) {
        int k2 = (int) (k % (base * n_levels[I[j]]));
        if (!ISNA(X[Is[j]*n + i])) {
          idx_I_obs[n_idx_I_obs++] = j;
          joint_level_Is_obs = joint_level_Is_obs + base * ((int) (X[Is[j]*n + i] - 1.0)); /* ASSUMING DISCR. LEVELS > 0 */
          joint_level_k = joint_level_k + base * k2;
        }
        base = base * n_levels[I[j]];
      }

      if (n_idx_I_obs > 0 && joint_level_Is_obs != joint_level_k) {
        pr_i[i] = 0.0;
      } else {
        n_index_S = 0;
        index_S = Calloc(n_joint_levels, int);
        if (n_idx_I_obs == 0) {
          for (j=0; j < n_joint_levels; j++)
            index_S[j] = j;
          n_index_S = n_joint_levels;
        } else {
          for (j=0; j < n_idx_I_obs; j++) {
            if (idx_I_obs[j] == joint_level_k) /* is this correct ???? */
              index_S[n_index_S++] = j;
          }
          error("implementation not finished yet\n");
        }
        n_index_Si = 0;
        index_Si = Calloc(n_joint_levels, int);
        for (j=0; j < n_index_S; j++)
          if (index_S[j] == joint_level_k)
            index_Si[n_index_Si++] = index_S[j];
        if (n_index_Si == n_index_S)
          pr_i[i] = 1.0;
        else {
          /* K_den <- sum(sapply(index_S, function(i) {Ki(x, Ys, i, mapX2Y, Sigma, mu, p)})) */
          K_den = 0.0;
          for (j=0; j < n_index_S; j++)
            K_den = K_den + Ki(X, p, n, i, Y, n_Y, index_S[j], Sigma, mu, n_joint_levels, pr);
          pr_i[i] = 0.0;
          if (K_den != 0) {
            for (j=0; j < n_index_Si; j++)
              pr_i[i] = pr_i[i] + Ki(X, p, n, i, Y, n_Y, index_Si[j], Sigma, mu, n_joint_levels, pr) / K_den;
          }
        }
        Free(index_S);
        Free(index_Si);
      }
    }
  }

  Free(idx_I_obs);

  return pr_i;
}

suf_stats_t
stat_mis(double* X, int p, int n, int* missing_mask, int n_mis, int* I, int n_I,
         int* Is, int n_Is, int *Y, int n_Y, int* n_levels, int n_joint_levels,
         com_stats_t com_stats, double* pr, double* mu, double* Sigma) {
  int         i,j,k;
  suf_stats_t ss;
  double*     m;
  double*     h;
  double*     bar_y;
  double*     K;

  m = Calloc(n_joint_levels, double);
  h = Calloc(n_joint_levels*n_Y, double);
  bar_y = Calloc(n_joint_levels*n_Y, double);
  K = Calloc((n_Y+(n_Y+1))/2, double);

  ss = new_suf_stats(n_joint_levels, n_Y);

  /* QUESTION: DO WE NEED TO MAP Y TO Sigma ?????? (AS IN THE R CODE?) */
  /* YES!! JUST AS IN qp_ci_test_hmgm() AND WE SHOULD USE symmatsubm() !!!! */
  /* HOW DO WE DO THE LEVEL MATCHING FROM THE i-th JOINT LEVEL TO WHATEVER OTHER ONE ?? */ 
  if (n_Is > 0) {
  }

  Free(K);
  Free(bar_y);
  Free(h);
  Free(m);

  return ss;
}

/* this function assumes that all discrete variables specified in I and all continuous
 * variables specified in Y are involved in the calculations, i.e., Y~I~c(i,j,Q) */
static double
lr_em(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
      int i, int j, int* Q, int q, double tol) {
  double*     pr0;
  double*     pr;
  double*     mu0;
  double*     mu;
  double*     Sigma0;
  double*     Sigma;
  com_stats_t comStats_i, comStats_j, comStats_ij;
  suf_stats_t sufStats_i, sufStats_j, sufStats_ij;
  int*        Y_i = NULL;
  int*        Y_j = NULL;
  int*        I_ij = NULL;
  int*        Y_ij = NULL;
  int*        missing_mask;
  double      mdiff;
  int         k, l, n_I_i, n_Y_i, n_Y_j, n_I_ij, n_Y_ij;
  int         n_mis, n_upper_tri, n_joint_levels=1;
  int         n_joint_levels_i, n_joint_levels_j, n_joint_levels_ij;

  error("EM not implemented yet in C code\n");

  if (missing_obs(X, p, n, Y, n_Y, NULL, n))
    error("EM not implemented yet for missing continuous values\n");

  comStats_i = comStats_j = comStats_ij = (com_stats_t) {NULL, NULL, NULL};
  missing_mask = Calloc(n, int);
  n_mis = find_missing_obs(X, p, n, I, n_I, NULL, n, missing_mask);

  l = -1;
  for (k=0; k < n_I; k++) {
    n_joint_levels = n_joint_levels * n_levels[I[k]];
    if (I[k] == i)
      l = k;
  }
  n_upper_tri = (n_Y * (n_Y+1)) / 2;

  /* I_i <- intersect(I, c(j, Q)) */
  /* I_j <- intersect(I, c(i, Q)) */
  /* we'll be assuming that j is never discrete, and therefore, I_j = I and Y_j = Y \ j */
  /* exchange i with the last in I, if i is in I */

  n_joint_levels_i = n_joint_levels_j = n_joint_levels_ij = n_joint_levels;
  n_I_i = n_I;
  if (l >= 0) {
    int tmp;

    n_joint_levels_i = n_joint_levels / n_levels[I[l]];
    n_joint_levels_ij = n_joint_levels / n_levels[I[l]];
    tmp = I[n_I-1];
    I[n_I-1] = I[l];
    I[l] = tmp;
    n_I_i--; /* we'll re-use n_I with n_I_i */
  }

  pr0 = Calloc(n_joint_levels, double); /* WATCH OUT, this grows exponentially in the # of discrete vars */
  pr = Calloc(n_joint_levels, double); /* WATCH OUT, this grows exponentially in the # of discrete vars */
  for (k=0; k < n_joint_levels; k++) /* pr0 initialized to the uniform distribution */
    pr[k] = pr0[k] = 1.0 / n_joint_levels;

  /* mu0 initialized to zeroes, assuming that memset initializes content to zeroes */
  mu0 = Calloc(n_Y * n_joint_levels, double);
  mu = Calloc(n_Y * n_joint_levels, double);
  Sigma0 = Calloc(n_upper_tri, double); /* storing the upper triangle, incl. diagonal, only */
  Sigma = Calloc(n_upper_tri, double); /* storing the upper triangle, incl. diagonal, only */

  /* Y_i <- intersect(Y, c(j, Q)) */
  Y_i = Calloc(n_Y, int);
  n_Y_i = 0;
  l = -1;
  /* Sigma0 initialized to the diagonal matrix, simultaneously, we initialize Y_i and find j in Y */
  for (k=0; k < n_Y; k++) {
    Sigma[UTE2I(k, k)] = Sigma0[UTE2I(k, k)] = 1.0;
    if (Y[k] == j)
      l = k;
    if (Y[k] != i)
      Y_i[n_Y_i++] = Y[k];
  }

  /* Y_j <- intersect(Y, c(i, Q)) */
  n_Y_j = n_Y - 1; /* re-use Y for Y_j by putting j at the end of Y */
  if (l >= 0) {
    int tmp;

    tmp = Y[n_Y-1];
    Y[n_Y-1] = Y[l];
    Y[l] = tmp;
  }

  I_ij = Calloc(q, int);
  Y_ij = Calloc(q, int);
  /* I_ij <- intersect(I, Q)
   * Y_ij <- intersect(Y, Q) */
  n_I_ij = n_Y_ij = 0;
  for (k=0; k < q; k++) {
    l = 0;
    while (l < n_I && I[l] != Q[k])
      l++;
    if (l < n_I)
      I_ij[n_I_ij++] = Q[k];
    l = 0;
    while (l < n_Y && Y[l] != Q[k])
      l++;
    if (l < n_Y)
      Y_ij[n_Y_ij++] = Q[k];
  }

  if (n_I_i > 0) {
    /* comStat_i <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
     *                       Is=I_i, Ys=Y_i, levels_Is=levels_I_i) */
    /* note here we're reusing I for I_i via n_I_i */
    comStats_i = stat_com(X, p, n, missing_mask, n_mis, I, n_I_i, Y_i, n_Y_i,
                          n_levels, n_joint_levels_i);
    Rprintf("Es_com_i=\n");
    for (k=0; k < n_joint_levels_i; k++) {
      for (l=0; l < n_Y; l++)
        Rprintf("%.5f\t", comStats_i.Es_com[l*n_joint_levels_j+k]);
      Rprintf("\n");
    }
    Rprintf("Ess_com_i=\n");
    for (k=0; k < n_Y_i; k++) {
      for (l=0; l < n_Y_i; l++)
        Rprintf("%.5f\t", comStats_i.Ess_com[l*n_Y+k]);
      Rprintf("\n");
    }
    Rprintf("n_com_i=\n");
    for (k=0; k < n_joint_levels_i; k++)
      Rprintf("%d\t", comStats_i.n_com[k]);
    Rprintf("\n");
  }

  if (n_I > 0) { /* assuming n_I_j = n_I because j should never be discrete */
    /* comStat_j <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
     *                       Is=I_j, Ys=Y_j, levels_Is=levels_I_j) */
    /* note here we're reusing Y for Y_j via n_Y_j */
    comStats_j = stat_com(X, p, n, missing_mask, n_mis, I, n_I, Y, n_Y_j,
                          n_levels, n_joint_levels_j);
    Rprintf("Es_com_j=\n");
    for (k=0; k < n_joint_levels_j; k++) {
      for (l=0; l < n_Y_j; l++)
        Rprintf("%.5f\t", comStats_j.Es_com[l*n_joint_levels_j+k]);
      Rprintf("\n");
    }
    Rprintf("Ess_com_j=\n");
    for (k=0; k < n_Y_j; k++) {
      for (l=0; l < n_Y_j; l++)
        Rprintf("%.5f\t", comStats_j.Ess_com[l*n_Y+k]);
      Rprintf("\n");
    }
    Rprintf("n_com_j=\n");
    for (k=0; k < n_joint_levels_j; k++)
      Rprintf("%d\t", comStats_j.n_com[k]);
    Rprintf("\n");
  }

  if (n_I_ij > 0) {
    /* comStat_ij <- stat_com(X, idxCompleteObs, idxMissingObs, mapAllObs2MissingObs,
     *                        Is=I_ij, Ys=Y_ij, levels_Is=levels_I_ij) */
    comStats_ij = stat_com(X, p, n, missing_mask, n_mis, I_ij, n_I_ij, Y_ij, n_Y_ij,
                           n_levels, n_joint_levels_ij);
    Rprintf("Es_com_ij=\n");
    for (k=0; k < n_joint_levels_ij; k++) {
      for (l=0; l < n_Y_ij; l++)
        Rprintf("%.5f\t", comStats_ij.Es_com[l*n_joint_levels_ij+k]);
      Rprintf("\n");
    }
    Rprintf("Ess_com_ij=\n");
    for (k=0; k < n_Y_ij; k++) {
      for (l=0; l < n_Y_ij; l++)
        Rprintf("%.5f\t", comStats_ij.Ess_com[l*n_Y+k]);
      Rprintf("\n");
    }
    Rprintf("n_com_ij=\n");
    for (k=0; k < n_joint_levels_ij; k++)
      Rprintf("%d\t", comStats_ij.n_com[k]);
    Rprintf("\n");
  }

  /* first round */
  mdiff = 1.0;
  while (mdiff > tol) {
    sufStats_i = stat_mis(X, p, n, missing_mask, n_mis, I, n_I, I, n_I_i, Y_i,
                          n_Y_i, n_levels, n_joint_levels, comStats_i,
                          pr, mu, Sigma);
    sufStats_j = stat_mis(X, p, n, missing_mask, n_mis, I, n_I, I, n_I, Y,
                          n_Y_j, n_levels, n_joint_levels, comStats_j,
                          pr, mu, Sigma);
    sufStats_ij = stat_mis(X, p, n, missing_mask, n_mis, I, n_I, I_ij, n_I_ij,
                           Y_ij, n_Y_ij, n_levels, n_joint_levels, comStats_ij,
                           pr, mu, Sigma);
    mdiff = 0.0;
  }

  free_com_stats(comStats_i);
  free_com_stats(comStats_j);
  free_com_stats(comStats_ij);
  Free(Y_ij);
  Free(I_ij);
  Free(Y_i);
  Free(Sigma);
  Free(Sigma0);
  Free(mu);
  Free(mu0);
  Free(pr);
  Free(pr0);
  Free(missing_mask);

  return R_NaN;
}

static double
qp_ci_test_hmgm(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y,
                int n_Y, double* ucond_ssd, int* mapX2ucond_ssd, int i, int j,
                int* Q, int q, int use, double tol, double* df, double* a,
                double* b, int* n_co) {
  int     k,l;
  int     n_I_int = 0;
  int     n_Y_int = 0;
  int     total_n_Y;
  int     mixed_edge = FALSE;
  int     n_joint_levels = 1;
  int     n_joint_levels_i = 1;
  int     n_levels_i = 0;
  double  lr = 0.0;

  *n_co = n;

  /* if any of (i, j) is discrete it should be i */
  k = 0;
  while (k < n_I && j != I[k])
    k++;

  if (k < n_I) {
    int tmp;

    tmp = i;
    i = j;
    j = tmp;
  }
/*
  Rprintf("\n%d ci %d\n", i+1, j+1);
*/

  /* I <- intersect(I, c(i, Q)) */

  k = 0;
  while (k < n_I && i != I[k])
    k++;

  if (k < n_I) {
    I[k] = I[0];
    I[0] = i;
    n_levels_i = n_levels[i];
    n_joint_levels = n_joint_levels * n_levels[i];

    n_I_int++;
    mixed_edge = TRUE;
  }

  for (k=0; k < q; k++) {
    l = 0;
    while (l < n_I && Q[k] != I[l])
      l++;

    if (l < n_I) {
      I[l] = I[n_I_int];
      I[n_I_int] = Q[k];

      n_joint_levels = n_joint_levels * n_levels[I[n_I_int]];
      n_joint_levels_i = n_joint_levels_i * n_levels[I[n_I_int]];

      n_I_int++;
    }
  }

  n_I = n_I_int;

  /* Y <- intersect(Y, c(i, j, Q)) */
  if (n_I > 0) {
    if (I[0] != i) {  /* i is continuous */
      k = 0;
      while (k < n_Y && i != Y[k])
        k++;

      if (k < n_Y) {
        Y[k] = Y[0];
        Y[0] = i;
        n_Y_int++;
      }
    }
  } else { /* no discrete variables, then i is continuous */
    k = 0;
    while (k < n_Y && i != Y[k])
      k++;

    if (k < n_Y) {
      Y[k] = Y[0];
      Y[0] = i;
      n_Y_int++;
    }
  }

  k = 0;
  while (k < n_Y && j != Y[k])
    k++;

  if (k < n_Y) {
    Y[k] = Y[n_Y_int];
    Y[n_Y_int] = j;
    n_Y_int++;
  }


  for (k=0; k < q; k++) {
    l = 0;
    while (l < n_Y && Q[k] != Y[l])
      l++;

    if (l < n_Y) {
      Y[l] = Y[n_Y_int];
      Y[n_Y_int] = Q[k];
      n_Y_int++;
    }
  }

  total_n_Y = n_Y;
  n_Y = n_Y_int;

  lr = R_NaN;
  switch(use) {
    case USE_COMPLETE_OBS:
      lr = lr_complete_obs(X, p, n, I, n_I, n_levels, Y, n_Y, ucond_ssd,
                           mapX2ucond_ssd, total_n_Y, i, j, Q, q, n_co);
      break;
    case USE_EM:
      lr = lr_em(X, p, n, I, n_I, n_levels, Y, n_Y, i, j, Q, q, tol);
      break;
    default:
      error("wrong 'use' argument in the internal call to qp_ci_test_hmgm()\n");
  }

  *df = 1.0;
  *a = ((double) ((use == USE_EM ? n : *n_co) - n_Y - n_joint_levels + 1)) / 2.0;
  if (mixed_edge) {
    *df = ((double) (n_joint_levels_i * (n_levels_i - 1)));
    *b = *df / 2.0;
  } else
    *b = 0.5;

  return lr;
}



/*
  FUNCTION: qp_ci_test_hmgm_sml
  PURPOSE: perform a test for conditional independence between variables
           indexed by i and j given the conditioning set Q. it differs from
           qp_ci_test_hmgm() in that it takes a list of snpMatrix objects as
           part of the input and tries to access efficiently that information
  RETURNS: a list with two members, the t-statistic value and the p-value
           on rejecting the null hypothesis of independence
*/

static double
qp_ci_test_hmgm_sml(SEXP Xsml, int* cumsum_sByChr, int s, int gLevels, double* XEP1q,
                    int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
                    double* ucond_ssd, int* mapX2ucond_ssd, int i, int j, int* Q, int q,
                    int use, double tol, double* df, double* a, double* b, int* n_co) {
  int     nChr = length(Xsml);
  int     k,l;
  int*    I_int;
  int     n_I_int = 0;
  int     n_Y_int = 0;
  int     n_I_i, n_Y_i, n_Y_j, n_Y_ij;
  int     total_n_Y;
  int*    n_levels_int;
  int     s1q=0;
  int     sign;
  int     final_sign = 1;
  int     flag_zero = 0;
  int     mixed_edge = FALSE;
  int     n_joint_levels = 1;
  int     n_joint_levels_i = 1;
  int     n_levels_i = 0;
  int*    idx_misobs = NULL;
  double* ssd_mat;
  double  x;
  double  lr = 0.0;

  *n_co = n;

  /* assume if any of (i, j) is discrete is always i */

  /* I <- intersect(I, c(i, Q)) */
  I_int = Calloc(q+1, int);
  /* n_levels_int = Calloc(q+1, int); */
  n_levels_int = Calloc(p+q+1, int);

  k = 0;
  if (i < p)
    while (k < n_I && i != I[k])
      k++;

  if (k < n_I || i >=p) { /* i is discrete */
    if (i < p) {
      I_int[n_I_int] = i;
      /* n_levels_int[n_I_int] = n_levels[k]; */
      n_levels_int[i] = n_levels[i];
    } else {
      int selChr=1;
      SEXP smR;
      const unsigned char* sm_i;

      while (cumsum_sByChr[selChr] < i-p && selChr < nChr)
        selChr++;

      if (selChr > nChr)
        error("chromosome not found for SNP in i\n");

      smR = VECTOR_ELT(Xsml, selChr-1);
      sm_i = RAW(smR) + (i-p-cumsum_sByChr[selChr])*n;
      for (l=0; l < n; l++) {
        int g = (int) sm_i[l];

        if (sm_i[l] < 4)
          XEP1q[p*n+l] = (double) g;
        else
          XEP1q[p*n+l] = NA_REAL;
      }
      /* n_levels_int[n_I_int] = gLevels; */
      n_levels_int[p+n_I_int] = gLevels;
      I_int[n_I_int] = p+s1q;
      i = p+s1q;
      s1q++;
    }

    /* n_levels_i = n_levels_int[n_I_int]; */
    n_levels_i = n_levels_int[p+n_I_int];
    /* n_joint_levels = n_joint_levels * n_levels_int[n_I_int]; */
    n_joint_levels = n_joint_levels * n_levels_int[p+n_I_int];
    n_I_int++;
    mixed_edge = TRUE;
  }

  for (k=0; k < q; k++) {
    l = 0;
    if (Q[k] < p)
      while (l < n_I && Q[k] != I[l])
        l++;

    if (l < n_I || Q[k] >= p) {
      if (Q[k] < p) {
        I_int[n_I_int] = Q[k];
        /* n_levels_int[n_I_int] = n_levels[l]; */
        n_levels_int[Q[k]] = n_levels[Q[k]];
      } else {
        int selChr=0;
        SEXP smR;
        const unsigned char* sm_Qk;

        while (cumsum_sByChr[selChr] < Q[k]-p && selChr < nChr)
          selChr++;

        if (selChr == nChr)
          error("chromosome not found for SNP in Q\n");

        smR = VECTOR_ELT(Xsml, selChr-1);
        sm_Qk = RAW(smR) + (Q[k]-p-cumsum_sByChr[selChr])*n;
        for (l=0; l < n; l++) {
          int g = (int) sm_Qk[l];

          if (sm_Qk[l] < 4)
            XEP1q[(p+s1q)*n+l] = (double) g;
          else
            XEP1q[(p+s1q)*n+l] = NA_REAL;
        }
        /* n_levels_int[n_I_int] = gLevels; */
        n_levels_int[p+n_I_int] = gLevels;
        I_int[n_I_int] = p+s1q;
        Q[k] = p+s1q;
        s1q++;
      }
      /* n_joint_levels = n_joint_levels * n_levels_int[n_I_int]; */
      n_joint_levels = n_joint_levels * n_levels_int[p+n_I_int];
      /* n_joint_levels_i = n_joint_levels_i * n_levels_int[n_I_int]; */
      n_joint_levels_i = n_joint_levels_i * n_levels_int[p+n_I_int];

      n_I_int++;
    }
  }

  /* Y <- intersect(Y, c(i, j, Q)) */
  if (n_I_int > 0) {
    if (I_int[0] != i) {  /* i is continuous */
      k = 0;
      while (k < n_Y && i != Y[k])
        k++;

      if (k < n_Y) {
        Y[k] = Y[0];
        Y[0] = i;
        n_Y_int++;
      }
    }
  } else { /* no discrete variables, then i is continuous */
    k = 0;
    while (k < n_Y && i != Y[k])
      k++;

    if (k < n_Y) {
      Y[k] = Y[0];
      Y[0] = i;
      n_Y_int++;
    }
  }

  k = 0;
  while (k < n_Y && j != Y[k]) /* assume j is always continuous */
    k++;

  if (k < n_Y) {
    Y[k] = Y[n_Y_int];
    Y[n_Y_int] = j;
    n_Y_int++;
  }


  for (k=0; k < q; k++) {
    l = 0;
    while (l < n_Y && Q[k] != Y[l])
      l++;

    if (l < n_Y) {
      Y[l] = Y[n_Y_int];
      Y[n_Y_int] = Q[k];
      n_Y_int++;
    }
  }

  total_n_Y = n_Y;
  n_Y = n_Y_int;

  /* the call below to ssd_A() assumes that this Calloc() memsets ssd_mat to zeroes */
  ssd_mat = Calloc((n_Y * (n_Y + 1)) / 2, double);  /* upper triangle includes the diagonal */

  if (n_I_int > 0 || ucond_ssd == NULL) {
    idx_misobs = Calloc(n, int);

    /* when calculating the larger ssd, the number of complete observations and the
     * logical mask of missing observations is pulled out and employed in later calls
     * to ssd_A */
    *n_co = ssd_A(XEP1q, p+s1q, n, I_int, n_I_int, n_levels_int, Y, n_Y, NULL, idx_misobs, ssd_mat);
  } else {
    int* tmp;

    tmp = Calloc(n_Y, int);
    for (k=0; k < n_Y; k++)
      tmp[k] = mapX2ucond_ssd[Y[k]];

    symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y);

    Free(tmp);
  }
/*
  Rprintf("ssd:\n");
  m = 0;
  for (k=0; k < n_Y; k++) {
    for (l=0; l <= k; l++) {
       Rprintf("%10.6f\t", ssd_mat[m++]);
    }
    Rprintf("\n");
  }
*/

  lr = x = symmatlogdet(ssd_mat, n_Y, &sign);
  if (x < -DBL_DIG)
    flag_zero = TRUE;
  final_sign *= sign;

/*
  Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y, &sign));
  Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/

  /* ssd_i = ssd_Gamma when i is discrete or ssd_{Gamma\i} when i is continuous */

  k = 0;
  while (i != I_int[k] && k < n_I_int)
    k++;

  n_I_i = n_I_int;
  n_Y_i = n_Y;
  if (k < n_I_int) {  /* i is discrete */
    I_int[k] = I_int[n_I_int-1];
    I_int[n_I_int-1] = i;
    n_I_i = n_I_int - 1;
  } else {       /* i is continuous */
    k = 0;
    while (i != Y[k] && k < n_Y)
      k++;

    if (k < n_Y) {
      Y[k] = Y[n_Y-1];
      Y[n_Y-1] = i;
      n_Y_i = n_Y - 1;
    } else
      error("qp_ci_test_hmgm_sml(): i does not form part of neither I nor Y\n");
  }

  if (n_I_int > 0 || ucond_ssd == NULL) {
    memset(ssd_mat, 0, ((n_Y_i * (n_Y_i + 1)) / 2) * sizeof(double));
    ssd_A(XEP1q, p+s1q, n, I_int, n_I_i, n_levels_int, Y, n_Y_i, idx_misobs, NULL, ssd_mat);
  } else {
    int* tmp;

    tmp = Calloc(n_Y_i, int);
    for (k=0; k < n_Y_i; k++)
      tmp[k] = mapX2ucond_ssd[Y[k]];

    symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y_i);

    Free(tmp);
  }

/*
  Rprintf("ssd_i:\n");
  m = 0;
  for (k=0; k < n_Y_i; k++) {
    Rprintf("%d", Y[k]+1);
    for (l=0; l <= k; l++) {
       Rprintf("\t%10.6f", ssd_mat[m++]);
    }
    Rprintf("\n");
  }
*/

  x = symmatlogdet(ssd_mat, n_Y_i, &sign);
  lr -= x;
  if (x < -DBL_DIG)
    flag_zero = TRUE;
  final_sign *= sign;
/*
  Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y_i, &sign));
  Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/

  if (n_Y > 1) {
    n_Y_j = n_Y;
    k = 0;
    while (j != Y[k] && k < n_Y)
      k++;

    if (k < n_Y) {
      Y[k] = Y[n_Y-1];
      Y[n_Y-1] = j;
      n_Y_j = n_Y - 1;
    }

    if (n_I_int > 0 || ucond_ssd == NULL) {
      memset(ssd_mat, 0, ((n_Y_j * (n_Y_j + 1)) / 2) * sizeof(double));
      ssd_A(XEP1q, p+s1q, n, I_int, n_I_int, n_levels_int, Y, n_Y_j, idx_misobs, NULL, ssd_mat);
    } else {
      int* tmp;

      tmp = Calloc(n_Y_j, int);
      for (k=0; k < n_Y_j; k++)
        tmp[k] = mapX2ucond_ssd[Y[k]];

      symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y_j);

      Free(tmp);
    }

/*
    Rprintf("ssd_j:\n");
    m = 0;
    for (k=0; k < n_Y_j; k++) {
      Rprintf("%d", Y[k]+1);
      for (l=0; l <= k; l++) {
         Rprintf("\t%10.6f", ssd_mat[m++]);
      }
      Rprintf("\n");
    }
*/

    x = symmatlogdet(ssd_mat, n_Y_j, &sign);
    lr -= x;
    if (x < -DBL_DIG)
      flag_zero = TRUE;
    final_sign *= sign;

    /*
    Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y_j, &sign));
    Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
    */

    n_Y_ij = n_Y_j;
    k = 0;
    while (i != Y[k] && k < n_Y)
      k++;

    if (k < n_Y) {
      Y[k] = Y[n_Y-2];
      Y[n_Y-2] = i;
      n_Y_ij = n_Y_j - 1;
    }

    if (n_Y_ij > 0) {
      if (n_I_int > 0 || ucond_ssd == NULL) {
        memset(ssd_mat, 0, ((n_Y_j * (n_Y_j + 1)) / 2) * sizeof(double));
        ssd_A(XEP1q, p+s1q, n, I_int, n_I_i, n_levels_int, Y, n_Y_ij, idx_misobs, NULL, ssd_mat);
      } else {
        int* tmp;

        tmp = Calloc(n_Y_ij, int);
        for (k=0; k < n_Y_ij; k++)
          tmp[k] = mapX2ucond_ssd[Y[k]];

        symmatsubm(ssd_mat, ucond_ssd, total_n_Y, tmp, n_Y_ij);

        Free(tmp);
      }

/*
      Rprintf("ssd_ij:\n");
      m = 0;
      for (k=0; k < n_Y_ij; k++) {
        Rprintf("%d", Y[k]+1);
        for (l=0; l <= k; l++) {
           Rprintf("\t%10.6f", ssd_mat[m++]);
        }
        Rprintf("\n");
      }
*/

      x = symmatlogdet(ssd_mat, n_Y_ij, &sign);
      lr += x;
      if (x < -DBL_DIG)
        flag_zero = TRUE;
      final_sign *= sign;

/*
      Rprintf("log(det(ssd_A))=%.5f\n", symmatlogdet(ssd_mat, n_Y_ij, &sign));
      Rprintf("sign(log(det(ssd_A)))=%d\n", sign);
*/
    }

  }

  Free(ssd_mat);
  Free(n_levels_int);
  Free(I_int);

  if (idx_misobs != NULL)
    Free(idx_misobs);

  if (flag_zero || final_sign == -1)
    lr = R_NaN;
  else
    lr = lr * ((double) -(*n_co));

  *df = 1.0;
  *a = ((double) (*n_co - n_Y - n_joint_levels + 1)) / 2.0; /* by now complete.obs w/ missing data */
  if (mixed_edge) {
    *df = ((double) (n_joint_levels_i * (n_levels_i - 1)));
    *b = *df / 2.0;
  } else
    *b = 0.5;

  return lr;
}



/*
  FUNCTION: qp_edge_nrr
  PURPOSE: estimate the non-rejection rate of the edge as the number of tests
           that accept the null hypothesis of independence given the
           q-order conditionals
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static double
qp_edge_nrr(double* X, double* S, int p, int n, int i, int j, int q,
            int* restrictQ, int n_rQ, int* fixQ, int n_fQ, int nTests, double alpha) {
  double thr;
  int*   q_by_T_samples;
  int*   ijQ = NULL;
  int*   Q = NULL;
  int    k, n_upper_tri;
  int    nAcceptedTests = 0;
  int    work_with_margin = FALSE;

  n_upper_tri = ( (q+2) * ((q+2)+1) ) / 2; /* this upper triangle includes the diagonal */
  if (S == NULL) {
    S = Calloc(n_upper_tri, double);
    ijQ = Calloc(q+2, int);
    Q = Calloc(q, int);
    ijQ[0] = i;
    ijQ[1] = j;
    for (k=0; k < q; k++)
      Q[k] = k+2;
    work_with_margin = TRUE;
  }

  q_by_T_samples = Calloc(q * nTests, int);

  if (n_rQ == 0)
    sampleQs(nTests, q, i, j, p, NULL, fixQ, n_fQ, q_by_T_samples);
  else
    sampleQs(nTests, q, i, j, n_rQ, restrictQ, fixQ, n_fQ, q_by_T_samples);

  thr = qt(1.0-(alpha/2.0), n-q-2, 1, 0);

  for (k = 0; k < nTests; k++) {
    double t_value;

    if (work_with_margin) {
      Memcpy((int*) (ijQ+2), (int*) (q_by_T_samples+k*q), (size_t) q);
      memset(S, 0, sizeof(double) * n_upper_tri);
      n = ssd(X, p, n, ijQ, q+2, NULL, n, TRUE, NULL, S);
      t_value = qp_ci_test_std(S, q+2, n, 0, 1, Q, q, NULL);
    } else
      t_value = qp_ci_test_std(S, p, n, i, j, (int*) (q_by_T_samples+k*q), q, NULL);

    if (fabs(t_value) < thr)
      nAcceptedTests++;
  }

  Free(q_by_T_samples);

  if (work_with_margin) {
    Free(S);
    Free(ijQ);
    Free(Q);
  }

  return (double) ( ((double) nAcceptedTests) / ((double) nTests) );
}



/*
  FUNCTION: qp_edge_nrr_identicalQs
  PURPOSE: estimate the non-rejection rate of the edge as the number of tests
           that accept the null hypothesis of independence given the
           q-order conditionals
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static double
qp_edge_nrr_identicalQs(double* S, int n_var, int* Qs, double* Qinvs, int N, int i,
                        int j, int q, int nTests, double alpha) {
  double thr;
  int    k;
  int    nAcceptedTests = 0;
  int    nActualTests = 0;
  double avgpr = 0;

  thr = qt(1.0-(alpha/2.0), N-q-2, 1, 0);

  for (k = 0; k < nTests; k++) {
    double t_value;
    int    l=0;

    while (l < q && Qs[k*q+l] != i && Qs[k*q+l] != j)
      l++;

    if (l >= q) {
      t_value = qp_ci_test_opt(S, n_var, N, i, j, (int*) (Qs+k*q), q,
                               (double*) (Qinvs+k*q*q), NULL);

      if (fabs(t_value) < thr)
        nAcceptedTests++;

      nActualTests++;
    } else
      avgpr++;
  }

  return (double) ( ((double) nAcceptedTests) / ((double) nActualTests) );
}



/*
  FUNCTION: qp_edge_nrr_hmgm
  PURPOSE: estimate the non-rejection rate of the edge as the number of tests
           that accept the null hypothesis of independence given the
           q-order conditionals for data where variables may be discrete or continuous
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static double
qp_edge_nrr_hmgm(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y,
                 int n_Y, double* ucond_ssd, int* mapX2ucond_ssd, int i, int j,
                 int q, int* restrictQ, int n_rQ, int* fixQ, int n_fQ, int nTests,
                 double alpha, int exactTest) {
  double thr=-1.0;
  double prev_thr, prev_df, prev_a, prev_b;
  int*   q_by_T_samples;
  int    k, n_co;
  int    nAcceptedTests = 0;
  int    nActualTests = 0;
  int*   problematicQ = NULL;

  q_by_T_samples = Calloc(q * nTests, int);

  if (n_rQ == 0)
    sampleQs(nTests, q, i, j, p, NULL, fixQ, n_fQ, q_by_T_samples);
  else
    sampleQs(nTests, q, i, j, n_rQ, restrictQ, fixQ, n_fQ, q_by_T_samples);

  if (n_I > 0) {
    k = 0;
    while (k < n_I && I[k] != j)
      k++;

    if (k < n_I) {  /* for a mixed edge i should be always the discrete variable */
      k = i;
      i = j;
      j = k;
    }
  }

  prev_thr = prev_df = prev_a = prev_b = -1;

  for (k = 0; k < nTests; k++) {
    double lambda, df, a, b;
/*
    int l;

    Rprintf("Q:");
    for (l=0; l < q; l++)
      Rprintf(" %d", *((int *) (q_by_T_samples+k*q+l)));
    Rprintf("\n");
*/
    lambda = qp_ci_test_hmgm(X, p, n, I, n_I, n_levels, Y, n_Y, ucond_ssd,
                             mapX2ucond_ssd, i, j, (int*) (q_by_T_samples+k*q),
                             q, USE_COMPLETE_OBS, 0.01, &df, &a, &b, &n_co);

    if (!ISNAN(lambda) && a > 0.0 && b > 0.0) {
      if (exactTest) {
        lambda = exp(lambda / ((double) -n));
        if (a == prev_a && b == prev_b)
          thr = prev_thr;
        else
          thr = qbeta(alpha, a, b, TRUE, FALSE);
        prev_a = a;
        prev_b = b;
        prev_thr = thr;
        if (lambda > thr)
          nAcceptedTests++;
      } else {
        if (df == prev_df)
          thr = prev_thr;
        else
          thr = qchisq(1.0-alpha, df, TRUE, FALSE);
        prev_df = df;
        prev_thr = thr;
        if (lambda < thr)
          nAcceptedTests++;
      }
      nActualTests++;
    } else
      problematicQ = (int*) (q_by_T_samples+k*q);
  }

  if (nActualTests < nTests) {
    char buf[4096];

    sprintf(buf, "Non-rejection rate estimation between i=%d and j=%d with q=%d was based on %d out of %d requested tests.\n"
                 "For instance, the CI test between i=%d and j=%d given Q={",
            i+1, j+1, q, nActualTests, nTests, i+1, j+1);
    for (k=0; k < q; k++) {
      char tmp[256];

      if (k == 0)
        sprintf(tmp, " %d", problematicQ[k]+1);
      else
        sprintf(tmp, ", %d", problematicQ[k]+1);
      strcat(buf, tmp);
    }
    strcat(buf, " }, could not be calculated. Try with smaller Q subsets or increase n if you can.\n");

    /* warning(buf); */
  }

  Free(q_by_T_samples);

  return (double) ( ((double) nAcceptedTests) / ((double) nActualTests) );
}



/*
  FUNCTION: qp_edge_nrr_hmgm_sml
  PURPOSE: estimate the non-rejection rate of the edge as the number of tests
           that accept the null hypothesis of independence given the
           q-order conditionals for data where variables may be discrete or continuous.
           The main difference with qp_edge_nrr_hmgm() is that this function works
           directly on the input snpMatrix object to avoid having to take the entire
           matrix of genotypes
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static double
qp_edge_nrr_hmgm_sml(SEXP X, int* cumsum_sByChr, int s, int gLevels, double* XEP1q,
                     int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
                     double* ucond_ssd, int* mapX2ucond_ssd, int i, int j, int q,
                     int* restrictQ, int n_rQ, int* fixQ, int n_fQ, int nTests,
                     double alpha, int exactTest) {
  double thr=-1.0;
  double prev_thr, prev_df, prev_a, prev_b;
  int*   q_by_T_samples;
  int    k, n_co;
  int    nAcceptedTests = 0;
  int    nActualTests = 0;
  int*   problematicQ = NULL;

  q_by_T_samples = Calloc(q * nTests, int);

  if (n_rQ == 0)
    sampleQs(nTests, q, i, j, p+s, NULL, fixQ, n_fQ, q_by_T_samples);
  else
    sampleQs(nTests, q, i, j, n_rQ, restrictQ, fixQ, n_fQ, q_by_T_samples);

  if (n_I > 0 || j >= p) {
    if (j < p) {
      k = 0;
      while (k < n_I && I[k] != j)
        k++;
    } else
      k = -1;

    if (k < n_I) {  /* for a mixed edge i should be always the discrete variable */
      k = i;
      i = j;
      j = k;
    }
  }

  prev_thr = prev_df = prev_a = prev_b = -1;

  for (k = 0; k < nTests; k++) {
    double lambda, df, a, b;
/*
    int l;

    Rprintf("Q:");
    for (l=0; l < q; l++)
      Rprintf(" %d", *((int *) (q_by_T_samples+k*q+l)));
    Rprintf("\n");
*/
    lambda = qp_ci_test_hmgm_sml(X, cumsum_sByChr, s, gLevels, XEP1q, p, n, I, n_I,
                                 n_levels, Y, n_Y, ucond_ssd, mapX2ucond_ssd, i, j,
                                 (int*) (q_by_T_samples+k*q), q, USE_COMPLETE_OBS,
                                 0.01, &df, &a, &b, &n_co);

    if (!ISNAN(lambda) && a > 0.0 && b > 0.0) {
      if (exactTest) {
        lambda = exp(lambda / ((double) -n));
        if (a == prev_a && b == prev_b)
          thr = prev_thr;
        else
          thr = qbeta(alpha, a, b, TRUE, FALSE);
        prev_a = a;
        prev_b = b;
        prev_thr = thr;
        if (lambda > thr)
          nAcceptedTests++;
      } else {
        if (df == prev_df)
          thr = prev_thr;
        else
          thr = qchisq(1.0-alpha, df, TRUE, FALSE);
        prev_df = df;
        prev_thr = thr;
        if (lambda < thr)
          nAcceptedTests++;
      }
      nActualTests++;
    } else
      problematicQ = (int*) (q_by_T_samples+k*q);
  }

  if (nActualTests < nTests) {
    char buf[4096];

    sprintf(buf, "Non-rejection rate estimation between i=%d and j=%d with q=%d was based on %d out of %d requested tests.\n"
                 "For instance, the CI test between i=%d and j=%d given Q={",
            i+1, j+1, q, nActualTests, nTests, i+1, j+1);
    for (k=0; k < q; k++) {
      char tmp[256];

      if (k == 0)
        sprintf(tmp, " %d", problematicQ[k]+1);
      else
        sprintf(tmp, ", %d", problematicQ[k]+1);
      strcat(buf, tmp);
    }
    strcat(buf, " }, could not be calculated. Try with smaller Q subsets or increase n if you can.\n");

    /* warning(buf); */
  }

  Free(q_by_T_samples);

  return (double) ( ((double) nAcceptedTests) / ((double) nActualTests) );
}


/*
  FUNCTION: qp_fast_edge_nrr
  PURPOSE: wrapper to the C function that estimates the non-rejection
           rate of the edge as the number of tests that accept the null
           hypothesis of independence given the q-order conditionals.
           It assumes that restrict.Q and fix.Q form a disjoint partition
           of all variables when fix.Q is set to something different from NULL
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static SEXP
qp_fast_edge_nrr(SEXP XR, SEXP SR, SEXP pR, SEXP nR, SEXP iR, SEXP jR, SEXP qR,
                 SEXP restrictQR, SEXP fixQR, SEXP nTestsR, SEXP alphaR) {
  int     i,j,k;
  int     p = INTEGER(pR)[0];
  int     n;
  int     q;
  int     nTests;
  double* X = NULL;
  double* S = NULL;
  double  alpha;
  int*    restrictQ=NULL;
  int     n_rQ=length(restrictQR);
  int*    fixQ=NULL;
  int     n_fQ=length(fixQR);
  SEXP    nrr;

  PROTECT_INDEX Xpi, Spi;

  if (XR != R_NilValue) {
    PROTECT_WITH_INDEX(XR, &Xpi);
    REPROTECT(XR = coerceVector(XR, REALSXP), Xpi);
    X = REAL(XR);
  }

  if (SR != R_NilValue) {
    PROTECT_WITH_INDEX(SR, &Spi);
    REPROTECT(SR = coerceVector(SR, REALSXP), Spi);
    S = REAL(SR);
  }

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;

  n = INTEGER(nR)[0];
  q = INTEGER(qR)[0];

  nTests = INTEGER(nTestsR)[0];

  alpha = REAL(alphaR)[0];

  if (i < 0 || i > p-1 || j < 0 || j > p-1)
    error("vertices of the selected edge (i=%d,j=%d) should lie in the range [1, p=%d]", i+1, j+1, p);

  if (q > p-2)
    error("q=%d > p-2=%d", q, p-2);

  if (q < 0)
    error("q=%d < 0",q);

  if (q > n-3)
    error("q=%d > n-3=%d", q, n-3);

  if (n_rQ > 0) {
    restrictQ = Calloc(n_rQ, int);
    for (k=0; k < n_rQ; k++)
      restrictQ[k] = INTEGER(restrictQR)[k]-1;
  }

  if (n_fQ > 0) {
    fixQ = Calloc(n_rQ, int);
    for (k=0; k < n_rQ; k++)
      fixQ[k] = INTEGER(fixQR)[k]-1;
  }

  PROTECT(nrr = allocVector(REALSXP,1));

  REAL(nrr)[0] = qp_edge_nrr(X, S, p, n, i, j, q, restrictQ, n_rQ,
                             fixQ, n_fQ, nTests, alpha);

  if (n_rQ > 0)
    Free(restrictQ);

  if (n_fQ > 0)
    Free(fixQ);

  UNPROTECT(1); /* nrr */

  if (XR != R_NilValue)
    UNPROTECT(1); /* XR */

  if (SR != R_NilValue)
    UNPROTECT(1); /* SR */

  return nrr;
}



/*
  FUNCTION: qp_fast_edge_nrr_hmgm
  PURPOSE: wrapper to the C function that estimates the non-rejection
           rate of the edge as the number of tests that accept the null
           hypothesis of independence given the q-order conditionals
           for data where variables may be discrete or continuous
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static SEXP
qp_fast_edge_nrr_hmgm(SEXP XR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP ssdR,
                      SEXP mapX2ssdR, SEXP iR, SEXP jR, SEXP qR, SEXP restrictQR,
                      SEXP fixQR, SEXP nTestsR, SEXP alphaR, SEXP exactTest) {
  int     n = INTEGER(getAttrib(XR, R_DimSymbol))[0];
  int     p = INTEGER(getAttrib(XR, R_DimSymbol))[1];
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int     q;
  int     i,j,k;
  int*    I;
  int*    Y;
  double* ssdMat = NULL;
  int*    mapX2ssd = NULL;
  int*    restrictQ = NULL;
  int     n_rQ = length(restrictQR);
  int*    fixQ=NULL;
  int     n_fQ=length(fixQR);
  int     nTests;
  double  alpha;
  SEXP    nrr;

  PROTECT_INDEX ssd_pi;

  if (ssdR != R_NilValue) {
    PROTECT_WITH_INDEX(ssdR, &ssd_pi);
    REPROTECT(ssdR = coerceVector(ssdR, REALSXP), ssd_pi);
    ssdMat = REAL(ssdR);
  }

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;

  q = INTEGER(qR)[0];

  nTests = INTEGER(nTestsR)[0];

  alpha = REAL(alphaR)[0];

  if (i < 0 || i > p-1 || j < 0 || j > p-1)
    error("vertices of the selected edge (i=%d,j=%d) should lie in the range [1, p=%d]", i+1, j+1, p);

  if (q > p-2)
    error("q=%d > p-2=%d", q, p-2);

  if (q < 0)
    error("q=%d < 0",q);

  if (q > n-3)
    error("q=%d > n-3=%d",q,n-3);

  I = Calloc(n_I, int);
  for (k=0; k < n_I; k++)
    I[k] = INTEGER(IR)[k]-1;

  Y = Calloc(n_Y, int);
  for (k=0; k < n_Y; k++)
    Y[k] = INTEGER(YR)[k]-1;

  if (ssdR != R_NilValue) {
    mapX2ssd = Calloc(p, int);
    for (k=0; k < p; k++)
      mapX2ssd[k] = INTEGER(mapX2ssdR)[k]-1;
  }

  if (n_rQ > 0) {
    restrictQ = Calloc(n_rQ, int);
    for (k=0; k < n_rQ; k++)
      restrictQ[k] = INTEGER(restrictQR)[k]-1;
  }

  if (n_fQ > 0) {
    fixQ = Calloc(n_rQ, int);
    for (k=0; k < n_rQ; k++)
      fixQ[k] = INTEGER(fixQR)[k]-1;
  }

  PROTECT(nrr = allocVector(REALSXP, 1));

  REAL(nrr)[0] = qp_edge_nrr_hmgm(REAL(XR), p, n, I, n_I, INTEGER(n_levelsR), Y,
                                  n_Y, ssdMat, mapX2ssd, i, j, q, restrictQ,
                                  n_rQ, fixQ, n_fQ, nTests, alpha, INTEGER(exactTest)[0]);

  UNPROTECT(1); /* nrr */

  if (ssdR != R_NilValue) {
    UNPROTECT(1); /* ssdR */
    Free(mapX2ssd);
  }

  if (n_rQ > 0)
    Free(restrictQ);

  if (n_fQ > 0)
    Free(fixQ);

  Free(I);
  Free(Y);

  return nrr;
}



/*
  FUNCTION: qp_fast_edge_nrr_hmgm_sml
  PURPOSE: wrapper to the C function that estimates the non-rejection
           rate of the edge as the number of tests that accept the null
           hypothesis of independence given the q-order conditionals
           for data where variables may be discrete or continuous. It
           differs from the previous function in that it takes the input
           data as an snpMatrix object (X) and a matrix of expression profiles
           and phenotypic data (XEP)
  RETURNS: the estimate of the non-rejection rate for the particular given edge
*/

static SEXP
qp_fast_edge_nrr_hmgm_sml(SEXP XR, SEXP cumsum_sByChrR, SEXP sR, SEXP gLevelsR,
                          SEXP XEPR, SEXP IR, SEXP n_levelsR, SEXP YR, SEXP ssdR,
                          SEXP mapX2ssdR, SEXP iR, SEXP jR, SEXP qR, SEXP restrictQR,
                          SEXP fixQR, SEXP nTestsR, SEXP alphaR, SEXP exactTest) {
  int     n = INTEGER(getAttrib(XEPR, R_DimSymbol))[0];
  int     p = INTEGER(getAttrib(XEPR, R_DimSymbol))[1];
  int     s = INTEGER(sR)[0];
  int     gLevels = INTEGER(gLevelsR)[0];
  int     n_I = length(IR);
  int     n_Y = length(YR);
  int     q;
  int     i,j,k;
  int*    I;
  int*    Y;
  double* XEP1q; /* store exp. profiles, phenotypic vars + up to (1+q) possible genotypes */
  double* ssdMat = NULL;
  int*    mapX2ssd = NULL;
  int*    restrictQ = NULL;
  int     n_rQ = length(restrictQR);
  int*    fixQ=NULL;
  int     n_fQ=length(fixQR);
  int     nTests;
  double  alpha;
  SEXP    nrr;

  PROTECT_INDEX ssd_pi;

  if (ssdR != R_NilValue) {
    PROTECT_WITH_INDEX(ssdR, &ssd_pi);
    REPROTECT(ssdR = coerceVector(ssdR, REALSXP), ssd_pi);
    ssdMat = REAL(ssdR);
  }

  i = INTEGER(iR)[0] - 1;
  j = INTEGER(jR)[0] - 1;

  q = INTEGER(qR)[0];

  nTests = INTEGER(nTestsR)[0];

  alpha = REAL(alphaR)[0];

  if (i < 0 || i > p+s-1 || j < 0 || j > p+s-1)
    error("vertices of the selected edge (i=%d,j=%d) should lie in the range [1, p=%d]", i+1, j+1, p+s);

  if (q > p+s-2)
    error("q=%d > p-2=%d", q, p+s-2);

  if (q < 0)
    error("q=%d < 0", q);

  if (q > n-3)
    error("q=%d > n-3=%d", q, n-3);

  I = Calloc(n_I, int);
  for (k=0; k < n_I; k++)
    I[k] = INTEGER(IR)[k]-1;

  Y = Calloc(n_Y, int);
  for (k=0; k < n_Y; k++)
    Y[k] = INTEGER(YR)[k]-1;

  if (ssdR != R_NilValue) {
    mapX2ssd = Calloc(p, int);
    for (k=0; k < p; k++)
      mapX2ssd[k] = INTEGER(mapX2ssdR)[k]-1;
  }

  if (n_rQ > 0) {
    restrictQ = Calloc(n_rQ, int);
    for (k=0; k < n_rQ; k++)
      restrictQ[k] = INTEGER(restrictQR)[k]-1;
  }

  if (n_fQ > 0) {
    fixQ = Calloc(n_rQ, int);
    for (k=0; k < n_rQ; k++)
      fixQ[k] = INTEGER(fixQR)[k]-1;
  }

  XEP1q = Calloc((p+q+1)*n, double); /* allocate extra space for all genotypes that may enter */
  Memcpy(XEP1q, REAL(XEPR), (size_t) (p * n));

  PROTECT(nrr = allocVector(REALSXP,1));

  REAL(nrr)[0] = qp_edge_nrr_hmgm_sml(XR, INTEGER(cumsum_sByChrR), s, gLevels,
                                      XEP1q, p, n, I, n_I, INTEGER(n_levelsR),
                                      Y, n_Y, ssdMat, mapX2ssd, i, j, q, restrictQ,
                                      n_rQ, fixQ, n_fQ, nTests, alpha, INTEGER(exactTest)[0]);

  UNPROTECT(1); /* nrr */

  if (ssdR != R_NilValue) {
    UNPROTECT(1); /* ssdR */
    Free(mapX2ssd);
  }

  if (n_rQ > 0)
    Free(restrictQ);

  if (n_fQ > 0)
    Free(fixQ);

  Free(XEP1q);
  Free(I);
  Free(Y);
  Free(mapX2ssd);

  return nrr;
}



/*
  FUNCTION: sampleQs
  PURPOSE: sample without replacement q elements from p, T times. this
           is a re-make of the SampleNoReplace function of random.c specifically
           tailored to sample in one shot everything we need
           if restrictQ != NULL then p is assumed to be its length and the
           sampling is restricted to the vertices in restrictQ taking care to
           remove first v_i and v_j if they occur in there.
  RETURN: a vector with of the T samples of q elements one after each other
*/

int
int_cmp(const void* a, const void* b);

static void
sampleQs(int T, int q, int v_i, int v_j, int p, int* restrictQ, int* fixQ,
         int n_fQ, int* y) {
  int  i;
  int  total_j = 0;
  int* x;
  int* z;

  if (restrictQ != NULL) {
    i = 0;
    while (i < p && v_i != restrictQ[i])
      i++;

    if (i < p) {
      restrictQ[i] = restrictQ[p-1];
      restrictQ[p-1] = v_i;
      p--;
    }

    i = 0;
    while (i < p && v_j != restrictQ[i])
      i++;

    if (i < p) {
      restrictQ[i] = restrictQ[p-1];
      restrictQ[p-1] = v_j;
      p--;
    }

    p = p + 2;
    v_i = v_j = -1;
  }

  x = Calloc(p,int);
  z = Calloc(p,int);

  for (i = 0; i < p; i++) {             /* x is a working-only vector */
    x[i] = i;
    z[i] = i;                           /* maps each vertex into a proper place */
  }

  if (v_i >=0 && v_j >=0) {               /* when sampling Qs outside v_i and v_j     */
    if (v_i < v_j) {                      /* we should take care that the mapping z   */
      z[v_i] = v_j != p-2 ? p-2 : p-1;    /* re-maps the v_i and v_j vertices to the  */
      z[v_j] = z[v_i] != p-1 ? p-1 : p-2; /* p-1 and p-2 properly when any of the two */
    } else {                              /* is smaller than p-2                      */
      z[v_j] = v_i != p-2 ? p-2 : p-1;
      z[v_i] = z[v_j] != p-1 ? p-1 : p-2;
    }
  }

  for (i = 0; i < T; i++) {
    int j;
    int m = p-2;                              /* we sample from p-2 elements */

    for (j = 0; j < q-n_fQ ; j++) {
      int r;

      r = (int) (((double) m) * unif_rand()); /* sample using R-builtin RNG */
      y[total_j + j] = x[r];
      x[r] = x[--m];                          /* sample without replacement */

    }

    /*Rprintf("Q:");*/
    for (j = total_j; j < total_j+q-n_fQ; j++) { /* replace again the sampled elements */
      x[y[j]] = y[j];                            /* for the next round of T       */
      y[j] = z[y[j]];                            /* use the mapping z to avoid choosing v_i or v_j */
      if (restrictQ != NULL)
        y[j] = restrictQ[y[j]];
      /*Rprintf(" %d", y[j]+1);*/
    }

    total_j += q-n_fQ;

    for (j = total_j ; j < total_j+n_fQ; j++) {
      y[j] = fixQ[j-total_j];
      /*Rprintf(" %d", y[j]+1);*/
    }
    /*Rprintf("\n");*/

    total_j += n_fQ;

  }

  Free(x);
  Free(z);
}



/*
  FUNCTION: qp_clique_number_lb
  PURPOSE: calculate a lower bound on the clique number from the input graph
  RETURNS: a lower bound of the clique number from the input graph
*/

typedef struct {
  int x;
  int ix;
} IntWithIdx;

int
int_cmp_desc_idx_decr(const void *a, const void *b) {
  const IntWithIdx* ia = (const IntWithIdx *) a;
  const IntWithIdx* ib = (const IntWithIdx *) b;

  return ib->x - ia->x;
}

static SEXP
qp_clique_number_lb(SEXP A, SEXP return_vertices, SEXP approx_iter, SEXP verbose) {
  int         n = INTEGER(getAttrib(A,R_DimSymbol))[0];
  IntWithIdx* deg;
  int*        pdeg;
  int*        ivec;
  int*        sset;
  int*        ssetelem;
  int         cliqueNumber=0;
  int*        cliqueVertices;
  int*        clq;
  int         i;
  int         ppct=-1;
  SEXP        return_value;

  PROTECT_INDEX Api;

  PROTECT_WITH_INDEX(A,&Api);

  REPROTECT(A = coerceVector(A,INTSXP),Api);

  deg = Calloc(n, IntWithIdx);
  pdeg = Calloc(n, int);
  ivec = Calloc(n, int);
  sset = Calloc(n, int);
  ssetelem = Calloc(n, int);
  cliqueVertices = Calloc(n, int);
  clq = Calloc(n, int);

  for (i=0; i < n; i++) {
    int j;

    deg[i].x = 0;
    for (j=0; j < n; j++)
      if (INTEGER(A)[j*n+i]) {
        deg[i].x++;
      }
    ivec[i] = i;
    deg[i].ix = i;
  } 
  qsort(deg, n, sizeof(IntWithIdx), int_cmp_desc_idx_decr);

  if (INTEGER(verbose)[0])
    Rprintf("calculating lower bound on the maximum clique size\n");

  for (i=0; i < INTEGER(approx_iter)[0]; i++) {
    int pct;
    int clqsze;
    int j;

    for (j=0; j < n; j++)
      pdeg[j] = deg[j].ix;

    if (i % n + 1 > 1) {
      int m;

      m = n;
      /* sample (i % n + 1) elements from n without replacement */
      for (j=0; j < i % n + 1; j++) {
        int r;

        r = (int) (((double) m) * unif_rand()); /* sample using R-builtin RNG */
        sset[j] = ivec[r];
        ivec[r] = ivec[--m];                    /* sample without replacement */
      }

      /* store the indices to the degree ranking using the sampled elements */
      for (j = 0; j < i % n + 1; j++) {  /* replace again the sampled elements */
        ssetelem[j] = pdeg[sset[j]];
        ivec[sset[j]] = sset[j];         /* so that ivec remains intact */
      }

      /* shuffle the stored indices using the Fisher-Yates algorithm */
      j = i % n;
      while (j >= 0) {
        int k;

        k = (int) (((double) j) * unif_rand());
        if (j != k) {
          int tmp;

          tmp = ssetelem[j];
          ssetelem[j] = ssetelem[k];
          ssetelem[k] = tmp;
        }
        j = j - 1;
      }

      /* shuffle the corresponding elements in the degree ranking */
      for (j=0; j < i % n + 1;j++)
        pdeg[sset[j]] = ssetelem[j];
    }

    /* go through the degree ranking building a clique */
    clq[0] = pdeg[0];
    clqsze = 1;
    for (j=1;j < n;j++) {
      int k;
      int isClique;

      clq[clqsze] = pdeg[j];

      isClique = 1;
      k = 0;
      while (k < (int) (((double) ((clqsze+1)*(clqsze)))/2.0) && isClique) {
        int u, v;

        i2e(k, &u, &v);
        if (!INTEGER(A)[clq[u] * n + clq[v]])
          isClique = 0;

        k++;
      }

      if (isClique)
        clqsze++; 
    }

    if (clqsze > cliqueNumber) {
      cliqueNumber = clqsze;
      Memcpy(cliqueVertices, clq, (size_t) clqsze);
    }

    pct = (int) ((i * 100) / INTEGER(approx_iter)[0]);
    if (pct != ppct) {
      if (INTEGER(verbose)[0]) {
        if (pct % 10 == 0)
          Rprintf("%d",pct);
        else
          Rprintf(".",pct);
        R_FlushConsole();
      }

      R_CheckUserInterrupt();
#ifdef Win32
      R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
      R_ProcessEvents();
#endif
      ppct = pct;
    }
  }

  UNPROTECT(1); /* A */

  if (INTEGER(verbose)[0])
    Rprintf("\n");

  if (INTEGER(return_vertices)[0]) {
    SEXP names;

    PROTECT(return_value = allocVector(VECSXP,2));

    SET_VECTOR_ELT(return_value,0,allocVector(INTSXP, 1));
    SET_VECTOR_ELT(return_value,1,allocVector(INTSXP, cliqueNumber));

    INTEGER(VECTOR_ELT(return_value,0))[0] = cliqueNumber;

    for (i=0;i < cliqueNumber;i++)
      INTEGER(VECTOR_ELT(return_value,1))[i] = cliqueVertices[i] + 1; /* in R vertices are 1-based */

    PROTECT(names = allocVector(VECSXP,2));
    SET_VECTOR_ELT(names, 0, mkChar("size"));
    SET_VECTOR_ELT(names, 1, mkChar("vertices"));
    setAttrib(return_value, R_NamesSymbol, names);
    UNPROTECT(1); /* names */

  } else {
    PROTECT(return_value = allocVector(INTSXP,1));

    INTEGER(return_value)[0] = cliqueNumber;
  }
  Free(deg);
  Free(pdeg);
  Free(ivec);
  Free(sset);
  Free(ssetelem);
  Free(cliqueVertices);
  Free(clq);

  UNPROTECT(1);  /* return_value */

  return return_value;
}



/*
  FUNCTION: qp_clique_number_os
  PURPOSE: wrapper of the function clique_unweighted_max_weight from the
           cliquer library that implements the algorithms by P. Ostergard
           to search for the maximal clique of maximum size
  RETURNS: the size of the largest clique
*/

static SEXP
qp_clique_number_os(SEXP I, SEXP return_vertices, SEXP verbose) {
  graph_t*       g;
  int            n = INTEGER(getAttrib(I,R_DimSymbol))[0];
  clique_options clq_opts;
  int            i;
  SEXP           return_value;

  if (INTEGER(verbose)[0])
    Rprintf("Niskanen and Ostergard algorithm for maximum clique running\n");

  PROTECT_INDEX Ipi;

  PROTECT_WITH_INDEX(I,&Ipi);

  REPROTECT(I = coerceVector(I,INTSXP),Ipi);

  g = graph_new(n);

  /* copy the incidence matrix 'I' into a graph_t structure
     from the cliquer library */

  for (i=0;i<n;i++) {
    int j;

    for (j=0;j<i;j++)
      if (INTEGER(I)[j*n+i])
        GRAPH_ADD_EDGE(g,i,j);
  }

  UNPROTECT(1); /* I */

  clq_opts.reorder_function   = reorder_by_default;
  clq_opts.reorder_map        = NULL;
  clq_opts.time_function      = INTEGER(verbose)[0] ? clique_print_time : NULL;
  clq_opts.output             = NULL;
  clq_opts.user_function      = NULL;
  clq_opts.user_data          = NULL;
  clq_opts.clique_list        = NULL;
  clq_opts.clique_list_length = 0;
  
  if (INTEGER(return_vertices)[0]) {
    int   i,j;
    set_t maxclq;
    SEXP  names;

    maxclq = clique_find_single(g,0,0,TRUE,&clq_opts);

    PROTECT(return_value = allocVector(VECSXP,2));

    SET_VECTOR_ELT(return_value,0,allocVector(INTSXP,1));
    SET_VECTOR_ELT(return_value,1,allocVector(INTSXP,set_size(maxclq)));

    INTEGER(VECTOR_ELT(return_value,0))[0] = set_size(maxclq);

    i=-1; j=0;
    while ((i=set_return_next(maxclq,i)) >= 0)
      INTEGER(VECTOR_ELT(return_value,1))[j++] = i + 1; /* in R vertices are 1-based */

    set_free(maxclq);

    PROTECT(names = allocVector(VECSXP,2));
    SET_VECTOR_ELT(names,0,mkChar("size"));
    SET_VECTOR_ELT(names,1,mkChar("vertices"));
    setAttrib(return_value,R_NamesSymbol,names);
    UNPROTECT(1); /* names */

  } else {
    PROTECT(return_value = allocVector(INTSXP,1));

    INTEGER(return_value)[0] = clique_unweighted_max_weight(g,&clq_opts);
  }

  UNPROTECT(1); /* return_value */

  graph_free(g);

  return return_value;
}



/*
  FUNCTION: cliquer_cb_add_clique_to_list
  PURPOSE: callback function for cliquer library function clique_unweighted_find_all
           it adds a clique to a linked list
  RETURNS: TRUE always to let continue the search till the end
*/

boolean
cliquer_cb_add_clique_to_list(set_t clique, graph_t* g, clique_options* opts) {
  clique_set_t* cset;
  clique_t*     c;

  cset = (clique_set_t *) opts->user_data;
  c = Calloc(1,clique_t);
  c->next = NULL;
                                                                                                
  if (cset->n == 0)
    cset->first = cset->last = c;
  else {
    clique_t* p;

    p = cset->last;
    p->next = cset->last = c;
  }

  c->u.vts = set_duplicate(clique);
  c->n     = set_size(clique);
  cset->n++;

  return TRUE;
}



/*
  FUNCTION: add_clique_vts
  PURPOSE: adds a clique to a linked list assuming that vertices
           are given and stored as cliquer sets
  RETURNS: nothing
*/

void
add_clique_vts(clique_set_t* cset, set_t clique) {
  clique_t* c;

  c = Calloc(1,clique_t);
  c->next = NULL;

  if (cset->n == 0)
    cset->first = cset->last = c;
  else {
    clique_t* p;

    p = cset->last;
    p->next = cset->last = c;
  }

  c->u.vts = set_duplicate(clique);
  c->n     = set_size(clique);
  cset->n++;
}



/*
  FUNCTION: add_clique_vta
  PURPOSE: adds a clique to a linked list assuming that vertices
           are given and stored as an array of integers
  RETURNS: nothing
*/

void
add_clique_vta(clique_set_t* cset, int* clique, int n) {
  clique_t* c;

  c = Calloc(1,clique_t);
  c->next = NULL;

  if (cset->n == 0)
    cset->first = cset->last = c;
  else {
    clique_t* p;

    p = cset->last;
    p->next = cset->last = c;
  }

  c->n     = n;
  c->u.vta = Calloc(n,int);
  Memcpy(c->u.vta,clique,(size_t) n);

  cset->n++;
}



/*
  FUNCTION: destroy_cliques_vts
  PURPOSE: destroys an object of the type clique_set_t which
           consists of going through a dynamically linked list
           and freeing the memory allocated for each of its elements
           this version assumes vertices were stored as cliquer sets
  RETURN: none
*/

void
destroy_cliques_vts(clique_set_t* cset) {
  clique_t* p;
                                                                                                
  if (cset->n == 0)
    return;
                                                                                                
  p = cset->first;
  while (p != NULL) {
    clique_t* tmp;
                                                                                                
    tmp = p->next;
    set_free(p->u.vts);
    Free(p);
    p = tmp;
  }

  cset->first = cset->last = NULL;
  cset->n = 0;
}



/*
  FUNCTION: destroy_cliques_vta
  PURPOSE: destroys an object of the type clique_set_t which
           consists of going through a dynamically linked list
           and freeing the memory allocated for each of its elements
           this version assumes vertices were stored as arrays of integers
  RETURN: none
*/

void
destroy_cliques_vta(clique_set_t* cset) {
  clique_t* p;
                                                                                                
  if (cset->n == 0)
    return;
                                                                                                
  p = cset->first;
  while (p != NULL) {
    clique_t* tmp;
                                                                                                
    tmp = p->next;
    Free(p->u.vta);
    Free(p);
    p = tmp;
  }

  cset->first = cset->last = NULL;
  cset->n = 0;
}



/*
  FUNCTION: init_cliques_list
  PURPOSE: initialize an object of the type clique_set_t
  RETURN: none
*/

void
init_cliques_list(clique_set_t* cset) {
  cset->first = cset->last = NULL;
  cset->n = 0;
}



int
int_cmp(const void *a, const void *b) {
  const int *ia = (const int *) a;
  const int *ib = (const int *) b;

  return *ia - *ib;
}



/*
  FUNCTION: qp_fast_cliquer_get_cliques
  PURPOSE: finds the (maximal) cliques of an undirected graph, it uses the
           library 'cliquer' from:

           Ostergard, PRJ. A fast algorithm for the maximum clique problem
           Discrete Appl. Math., 120:195-205, 2002
           http://users.tkk.fi/~pat/cliquer.html
  RETURNS: a list of (maximal) cliques
*/

static SEXP
qp_fast_cliquer_get_cliques(SEXP I, SEXP clqspervtx, SEXP verbose) {
  graph_t*       g;
  int            n = INTEGER(getAttrib(I,R_DimSymbol))[0];
  int            i;
  int            nclqs;
  clique_set_t   clqlst;
  SEXP           clqlstR;
  clique_options clq_opts;

  if (!isMatrix(I)) {
    error("qpGetCliques() expects an incidence matrix");
  }

  PROTECT_INDEX Ipi;

  PROTECT_WITH_INDEX(I,&Ipi);

  REPROTECT(I = coerceVector(I,INTSXP),Ipi);

  g = graph_new(n);

  /* copy the incidence matrix 'I' into a graph_t structure
     from the cliquer library */

  for (i=0;i<n;i++) {
    int j;

    for (j=0;j<i;j++)
      if (INTEGER(I)[j*n+i])
        GRAPH_ADD_EDGE(g,i,j);
  }

  UNPROTECT(1); /* I */

  init_cliques_list(&clqlst);
  clq_opts.reorder_function   = reorder_by_default;
  clq_opts.reorder_map        = NULL;
  clq_opts.time_function      = INTEGER(verbose)[0] ? clique_print_time : NULL;
  clq_opts.output             = NULL;
  clq_opts.user_function      = cliquer_cb_add_clique_to_list;
  clq_opts.user_data          = (void *) &clqlst;
  clq_opts.clique_list        = NULL;
  clq_opts.clique_list_length = 0;
  
  nclqs = clique_unweighted_find_all(g,1,n,TRUE,&clq_opts);

  graph_free(g);

  if (nclqs != clqlst.n)
    error("qpGetCliques: nclqs is different from clqlst.n!!!");

  /* allocate n components more where to put the cliques to which each vertex belongs to
     afterwards the intersection operation should quickly give the set of cliques including
     a given edge */

  if (INTEGER(clqspervtx)[0])
    PROTECT(clqlstR = allocVector(VECSXP,clqlst.n + n));
  else
    PROTECT(clqlstR = allocVector(VECSXP,clqlst.n));

  if (clqlst.n > 0) {
    clique_t* p;
    int       iclq;
    int**     idxclqs = NULL;
    int*      nidxclqs = NULL;
    int*      sidxclqs = NULL;

    if (INTEGER(clqspervtx)[0]) {
      int i;

      idxclqs = (int **) Calloc(n,int *);
      nidxclqs = (int *) Calloc(n,int);
      sidxclqs = (int *) Calloc(n,int);

      for (i=0;i<n;i++)
        nidxclqs[i]=0;
    }

    iclq = INTEGER(clqspervtx)[0] ? n : 0;
    p = clqlst.first;
    while (p != NULL) {
      clique_t* tmpp;
      int  v,i;
      SEXP clq;

      SET_VECTOR_ELT(clqlstR,iclq,clq = allocVector(INTSXP,p->n));

      v=-1; i=0;
      while ((v=set_return_next(p->u.vts,v)) >= 0) {
        INTEGER(VECTOR_ELT(clqlstR,iclq))[i] = v + 1; /* in R vertices are 1-based */

        if (INTEGER(clqspervtx)[0]) {
          if (nidxclqs[v] == 0) {
            sidxclqs[v] = 1;
            idxclqs[v] = (int *) Calloc(sidxclqs[v],int);
            idxclqs[v][nidxclqs[v]] = iclq + 1;
            nidxclqs[v]++;
          } else {
            if (sidxclqs[v] > nidxclqs[v]) {
              idxclqs[v][nidxclqs[v]] = iclq + 1;
              nidxclqs[v]++;
            } else {
              sidxclqs[v] = sidxclqs[v] * 2;
              idxclqs[v] = (int *) Realloc(idxclqs[v],sidxclqs[v],int);
              idxclqs[v][nidxclqs[v]] = iclq + 1;
              nidxclqs[v]++;
            }
          }
        }

        i++;
      }

      iclq++;

      /* free the elements of the linked list at the same time that the new R list
         structure is created to store the cliques in order to use as little memory as possible */

      tmpp = p->next;
      set_free(p->u.vts);
      Free(p);
      p = tmpp;
    }

    if (INTEGER(clqspervtx)[0]) {
      int i;

      for (i=0;i<n;i++) {
        qsort(idxclqs[i],nidxclqs[i],sizeof(int),int_cmp);
        SET_VECTOR_ELT(clqlstR,i,allocVector(INTSXP,nidxclqs[i]));
        Memcpy(INTEGER(VECTOR_ELT(clqlstR,i)),idxclqs[i],(size_t) nidxclqs[i]);
        Free(idxclqs[i]);
      }

      Free(sidxclqs);
      Free(nidxclqs);
      Free(idxclqs);
    }

  }

  UNPROTECT(1); /* clqlstR */

  return clqlstR;
}



/*
  FUNCTION: qp_fast_update_cliques_removing
  PURPOSE: modifies an input list of maximal cliques by removing one edge of the graph
  RETURNS: a list of maximal cliques
*/

static SEXP
qp_fast_update_cliques_removing(SEXP I, SEXP clqlstR, SEXP vR, SEXP wR, SEXP verbose) {
  SEXP new_clqlstR;
  int n = INTEGER(getAttrib(I,R_DimSymbol))[0];
  int v = INTEGER(coerceVector(vR,INTSXP))[0] - 1; /* internally we work with 0-based vertices */
  int w = INTEGER(coerceVector(wR,INTSXP))[0] - 1; /* internally we work with 0-based vertices */
  int nclqlstR = length(clqlstR);
  int nnew_clqlstR;
  int* clqs_v_w;
  int  nclqs_v,nclqs_w;
  int  nclqs_v_w;
  clique_set_t clqlst;
  int* new_clq_v;
  int* new_clq_w;
  set_t allvtc;
  clique_t* p;
  int removed_so_far;
  int i,inewclq;
  int** idxclqs;
  int*  nidxclqs;
  int*  sidxclqs;
  int ppct;

  if (!isMatrix(I)) {
    error("qpUpdateCliquesRemoving() expects an incidence matrix");
  }

  PROTECT_INDEX Ipi;

  PROTECT_WITH_INDEX(I,&Ipi);

  REPROTECT(I = coerceVector(I,INTSXP),Ipi);

  nclqs_v = length(VECTOR_ELT(clqlstR,v));
  nclqs_w = length(VECTOR_ELT(clqlstR,w));
  clqs_v_w = (int *) Calloc(nclqs_v+nclqs_w,int);

  if (INTEGER(verbose)[0]) {
    Rprintf("qpUpdateCliquesRemoving: initially there are %d maximal clique(s)\n",nclqlstR-n);
    Rprintf("qpUpdateCliquesRemoving: searching cliques to which the edge v-w belongs to (v belongs to %d, w belongs to %d)\n",nclqs_v,nclqs_w);
  }

  /* the cliques where the edge v-w is are the common ones to those where v and w belong to */
  nclqs_v_w = 0;
  for (i=0;i<nclqs_v;i++)
    if (bsearch(INTEGER(VECTOR_ELT(clqlstR,v))+i,INTEGER(VECTOR_ELT(clqlstR,w)),nclqs_w,sizeof(int),int_cmp))
        clqs_v_w[nclqs_v_w++] = INTEGER(VECTOR_ELT(clqlstR,v))[i] - 1; /* internally the clique index
                                                                          in the R list is 0-based */

  /* sort the cliques to be removed in ascending clique index order */
  qsort(clqs_v_w,nclqs_v_w,sizeof(int),int_cmp);

  INTEGER(I)[v*n+w] = INTEGER(I)[w*n+v] = 0; /* make sure the edge is removed */

  allvtc = set_new(n);
  for (i=0;i<n;i++)
    SET_ADD_ELEMENT(allvtc,i);

  new_clq_v = Calloc(n,int); /* store here each new clique that formerly contained vertex v */
  new_clq_w = Calloc(n,int); /* store here each new clique that formerly contained vertex w */

  if (INTEGER(verbose)[0])
    Rprintf("qpUpdateCliquesRemoving: going through the %d v-w-containing clique(s) and decide which one(s) to add\n",nclqs_v_w);

  /* go through all the cliques where the edge v-w forms part */

  init_cliques_list(&clqlst);
  for (i=0;i<nclqs_v_w;i++) {
    int clq_v_w_size = length(VECTOR_ELT(clqlstR,clqs_v_w[i]));
    set_t restvtc_v = set_duplicate(allvtc);
    set_t restvtc_w = set_duplicate(allvtc);
    int iclq_v = 0;
    int iclq_w = 0;
    int j;

    for (j=0;j<clq_v_w_size;j++) {
      int vtx = INTEGER(VECTOR_ELT(clqlstR,clqs_v_w[i]))[j] - 1; /* internally we work with 0-based vertices */

      if (vtx != v) {
        new_clq_v[iclq_v++] = vtx + 1;  /* in the R object vertices are 1-based */
        SET_DEL_ELEMENT(restvtc_v,vtx);
      }

      if (vtx != w) {
        new_clq_w[iclq_w++] = vtx + 1;  /* in the R object vertices are 1-based */
        SET_DEL_ELEMENT(restvtc_w,vtx);
      }
    }

    if (is_maximal_clique(INTEGER(I),n,new_clq_v,clq_v_w_size-1,restvtc_v))
      add_clique_vta(&clqlst,new_clq_v,clq_v_w_size-1);

    if (is_maximal_clique(INTEGER(I),n,new_clq_w,clq_v_w_size-1,restvtc_w))
      add_clique_vta(&clqlst,new_clq_w,clq_v_w_size-1);
  }

  Free(new_clq_v);
  Free(new_clq_w);

  UNPROTECT(1); /* I */

  /* remove the cliques that contain the edge v-w and add the new ones */

  nnew_clqlstR = nclqlstR + clqlst.n - nclqs_v_w;
  PROTECT(new_clqlstR=allocVector(VECSXP,nnew_clqlstR));

  if (INTEGER(verbose)[0])
    Rprintf("qpUpdateCliquesRemoving: going to remove %d clique(s) and add %d clique(s) ending with %d clique(s)\n",nclqs_v_w,clqlst.n,nnew_clqlstR-n);

  inewclq = n;
  removed_so_far = 0;
  for (i=n;i<nclqlstR;i++)
    if (removed_so_far >= nclqs_v_w || clqs_v_w[removed_so_far] != i)
      SET_VECTOR_ELT(new_clqlstR,inewclq++, Rf_duplicate(VECTOR_ELT(clqlstR,i)));
    else
      removed_so_far++;

  p = clqlst.first;
  while (p != NULL) {
    SEXP newclq;
    clique_t* tmpp;

    PROTECT(newclq = allocVector(INTSXP,p->n));
    Memcpy(INTEGER(newclq),p->u.vta,(size_t) p->n);
    SET_VECTOR_ELT(new_clqlstR,inewclq++,newclq);
    UNPROTECT(1); /* newclq */

    /* free the elements of the linked list at the same time that the new R list
       structure is created to store the cliques in order to use a little memory as possible */

    tmpp = p->next;
    Free(p->u.vta);
    Free(p);
    p = tmpp;
  }

  Free(clqs_v_w);

  if (INTEGER(verbose)[0])
    Rprintf("qpUpdateCliquesRemoving: rebuilding references to cliques\n");

  idxclqs = (int **) Calloc(n,int *);
  nidxclqs = (int *) Calloc(n,int);
  sidxclqs = (int *) Calloc(n,int);

  for (i=0;i<n;i++)
    nidxclqs[i]=0;

  ppct = -1;
  for (i=n;i<nnew_clqlstR;i++) {
    int j;

    for (j=0;j<length(VECTOR_ELT(new_clqlstR,i));j++) {
      v = INTEGER(VECTOR_ELT(new_clqlstR,i))[j] - 1;
      if (nidxclqs[v] == 0) {
        sidxclqs[v] = 1;
        idxclqs[v] = (int *) Calloc(sidxclqs[v],int);
        idxclqs[v][nidxclqs[v]] = i + 1;
        nidxclqs[v]++;
      } else {
        if (sidxclqs[v] > nidxclqs[v]) {
          idxclqs[v][nidxclqs[v]] = i + 1;
          nidxclqs[v]++;
        } else {
          sidxclqs[v] = sidxclqs[v] * 2;
          idxclqs[v] = (int *) Realloc(idxclqs[v],sidxclqs[v],int);
          idxclqs[v][nidxclqs[v]] = i + 1;
          nidxclqs[v]++;
        }
      }
    }

    if (INTEGER(verbose)[0]) {
      int pct = (int) ((i*100)/nnew_clqlstR);

      if (pct != ppct) {
        if (pct % 10 == 0)
          Rprintf("%d",pct);
        else
          Rprintf(".",pct);
        R_FlushConsole();
#ifdef Win32
        R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
  }
  if (INTEGER(verbose)[0])
    Rprintf("\n");

  for (i=0;i<n;i++) {
    SEXP clqrefs;

    qsort(idxclqs[i],nidxclqs[i],sizeof(int),int_cmp);
    PROTECT(clqrefs = allocVector(INTSXP,nidxclqs[i]));
    Memcpy(INTEGER(clqrefs),idxclqs[i],(size_t) nidxclqs[i]);
    Free(idxclqs[i]);

    SET_VECTOR_ELT(new_clqlstR,i,clqrefs);
    UNPROTECT(1); /* clqrefs */
  }

  Free(sidxclqs);
  Free(nidxclqs);
  Free(idxclqs);

  UNPROTECT(1); /* new_clqlstR */

  return new_clqlstR;
}



/*
  FUNCTION: is_maximal_clique
  PURPOSE: returns whether the clique in 'clq' is maximal. note that the
           vertices in clq come 1-based and vertices in noclq come 0-based
  RETURNS: TRUE if the clique is maximal; FALSE otherwise
*/
/* vertices in clq come 1-based vertices in noclq come 0-based */

Rboolean
is_maximal_clique(int* I, int n, int* clq, int cs, set_t noclq) {
  int i;
  Rboolean maximal = TRUE;

  i=-1;
  while ((i=set_return_next(noclq,i))>=0 && maximal) {
    int j=0;
    int allconnected = TRUE;

    while (j < cs && allconnected) {
      allconnected = allconnected && I[(clq[j]-1)*n+i];
      j++;
    }

    maximal = !allconnected;
  }

  return maximal;
}



/*
  FUNCTION: qp_fast_pac_se
  PURPOSE: calculate the standard errors for the edges of an undirected graphical
           Gaussian Markov model according to the method by:

           Roverato and Whittaker. Standard errors for the parameters of
           graphical Gaussian models, STATISTICS AND COMPUTING, 6:297-302, 1996

  PARAMETERS: S - estimate of the sample covariance matrix
              I - incidence matrix of the graph and thus it is assumed that the
                  diagonal is set to either 0s or FALSE truth values since there
                  should be no loops
  RETURNS: a matrix with the standard errors of the edges
  TODO: optimize!!!
*/

static SEXP
qp_fast_pac_se(SEXP Shat, SEXP I) {
  int  n_var = INTEGER(getAttrib(I,R_DimSymbol))[0];
  int  n_edges;
  int* r_nz;
  int* c_nz;
  int  i,j;
  int  rnz,cnz;
  double* tmp1;
  double* tmp2;
  double* tmp3;
  double* tmp4;
  double* Iss1;
  double* Iss2;
  double* Iss;
  double* Issinv;
  SEXP    SER;
  double* SE;
  PROTECT_INDEX Spi,Ipi;

  if (!isMatrix(Shat) || !isMatrix(I)) {
    error("qpPACSE: Shat or I is not a matrix!");
  }

  PROTECT_WITH_INDEX(Shat, &Spi);
  PROTECT_WITH_INDEX(I, &Ipi);

  REPROTECT(Shat = coerceVector(Shat,REALSXP), Spi);
  REPROTECT(I = coerceVector(I,INTSXP), Ipi);

  n_edges = 0;
  for (i=0;i<n_var;i++)
    for (j=0;j<=i;j++)
      if (INTEGER(I)[i+j*n_var] != 0)
        n_edges++;

  r_nz = Calloc(n_edges+n_var,int);
  c_nz = Calloc(n_edges+n_var,int);

  /* selection of row and column indices according to whether
     the cells in I have value 0 or not. indices of the diagonal
     are also selected */

  rnz = cnz = 0;
  for (i=0;i<n_var;i++) {
    for (j=0;j<=i;j++) {
      if (INTEGER(I)[i+j*n_var] != 0 || i==j) {
        r_nz[rnz++] = i;
        c_nz[cnz++] = j;
      }
    }
  }

  UNPROTECT(1); /* I */

  tmp1 = Calloc(cnz*cnz,double);
  tmp2 = Calloc(rnz*rnz,double);
  tmp3 = Calloc(cnz*rnz,double);
  tmp4 = Calloc(rnz*cnz,double);
  matsubm(tmp1,REAL(Shat),n_var,c_nz,cnz,c_nz,cnz);
  matsubm(tmp2,REAL(Shat),n_var,r_nz,rnz,r_nz,rnz);
  matsubm(tmp3,REAL(Shat),n_var,c_nz,cnz,r_nz,rnz);
  matsubm(tmp4,REAL(Shat),n_var,r_nz,rnz,c_nz,cnz);

  UNPROTECT(1); /* Shat */

  Iss1 = Calloc(cnz*cnz,double);
  Iss2 = Calloc(rnz*rnz,double);
  Iss  = Calloc(rnz*cnz,double);
  matscalarprod(Iss1,cnz,cnz,tmp1,tmp2);
  matscalarprod(Iss2,rnz,rnz,tmp3,tmp4);
  Free(tmp1); Free(tmp2); Free(tmp3); Free(tmp4);
  matsumf(Iss,rnz,cnz,Iss1,Iss2,1.0);
  Free(Iss1); Free(Iss2);

  Issinv = Calloc(rnz*cnz,double);
  matinv(Issinv, Iss, rnz, 0);
  Free(Iss);

  PROTECT(SER = allocMatrix(REALSXP,n_var,n_var));
  SE = REAL(SER);

  for (i=0;i<n_var;i++) {
    for (j=i;j<n_var;j++) {
      SE[i+n_var*j] = SE[j+n_var*i] = NA_REAL;
    }
  }

  for (i=0;i<rnz;i++)
    if (r_nz[i] != c_nz[i])
      SE[r_nz[i]+n_var*c_nz[i]] = SE[c_nz[i]+n_var*r_nz[i]] = Issinv[i*rnz+i];

  Free(Issinv);

  for (i=0;i<n_var;i++)
    SE[i+n_var*i] = NA_REAL;

  Free(r_nz); Free(c_nz);

  UNPROTECT(1); /* SER */

  return SER;
}



/*
  FUNCTION: qp_fast_ipf
  PURPOSE: implement the Iterative Proportional Fitting (IPF) algorithm. Part of the
           R code below has been borrowed from an implementation by Graham Wills in
           June of 1992
  PARAMETERS: vvR - input matrix (usually the sample variance-covariance matrix)
              clqlstR - list of (maximal) cliques
              tolR - tolerance under which the main loop stops
              verbose - when set to TRUE the algorithm shows progression
  RETURNS: the input matrix adjusted to the constraints of the list of cliques
*/

static SEXP
qp_fast_ipf(SEXP vvR, SEXP clqlstR, SEXP tolR, SEXP verbose) {
  int           n = INTEGER(getAttrib(vvR,R_DimSymbol))[0];
  int           nclqlst = length(clqlstR);
  double*       V;
  double*       Vold;
  double*       vv;
  double*       tmp;
  double        tol;
  double        precision;
  SEXP          VR;
  Rboolean      clqspervtx;
  int           fstclq = 0;
  PROTECT_INDEX vvpi,tolpi;
  int           pct,ppct;
  int           i,j;

  if (!isMatrix(vvR) || !isNewList(clqlstR)) {
    error("fast.ipf expects a matrix and a list of cliques");
  }

  clqspervtx = INTEGER(VECTOR_ELT(clqlstR,0))[0] > n ? TRUE : FALSE;
  if (clqspervtx && nclqlst <= n+1)
    error("qpIPF: the clique list seem to have vertex-clique indexes but then shows too few cliques\n");
  else if (clqspervtx)
    fstclq = n;

  PROTECT_WITH_INDEX(vvR,&vvpi);
  PROTECT_WITH_INDEX(tolR,&tolpi);

  REPROTECT(vvR  = coerceVector(vvR,REALSXP),vvpi);
  REPROTECT(tolR = coerceVector(tolR,REALSXP),tolpi);
  
  V    = Calloc(n*n,double);
  Vold = Calloc(n*n,double);
  vv   = Calloc(n*n,double);
  tmp  = Calloc(n*n,double);
  tol  = REAL(tolR)[0];

  UNPROTECT(1); /* tol */

  /* V <- diag(length(vv[, 1])) */
  for (i=0;i<n;i++)
    for (j=0;j<n;j++) {
      vv[i+j*n] = REAL(vvR)[i+j*n];
      V[i+j*n] = i == j ? 1.0 : 0.0;
    }
  UNPROTECT(1); /* vvR */

  if (INTEGER(verbose)[0])
    Rprintf("qpIPF: %d cliques\n",nclqlst-fstclq);

  precision = DBL_MAX;
  while (precision > tol) {
    /* Vold <- V */
    for (i=0;i<n;i++)
      for (j=0;j<n;j++)
        Vold[i+n*j] = V[i+n*j];

    if (INTEGER(verbose)[0])
      Rprintf("Iterating through cliques (%): ");

    ppct = -1;
    for (i=fstclq;i<nclqlst;i++) {
      SEXP cR;
      int  j;
      int* a;
      int  csze;

      PROTECT(cR = coerceVector(VECTOR_ELT(clqlstR,i),INTSXP));
      csze = length(cR);
      a = Calloc(csze,int);
      for (j=0;j<csze;j++)
        a[j] = INTEGER(cR)[j]-1;
      UNPROTECT(1); /* cR */

      fast_ipf_step(n,vv,V,a,csze);
      Free(a);

      pct = (int) (((i-fstclq)*100)/(nclqlst-fstclq));
      if (pct != ppct) {
        if (INTEGER(verbose)[0]) {
          if (pct % 10 == 0)
            Rprintf("%d",pct);
          else
            Rprintf(".",pct);
          R_FlushConsole();
        }

        R_CheckUserInterrupt();
#ifdef Win32
        R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
    if (INTEGER(verbose)[0])
      Rprintf("\n");

    matsumf(tmp,n,n,V,Vold,-1.0);
    precision = matmxab(tmp,n,n);
    if (INTEGER(verbose)[0])
      Rprintf("Precision: %.10f\n", precision);
  }

  PROTECT(VR = allocMatrix(REALSXP,n,n));
  Memcpy(REAL(VR),V,(size_t) (n * n));
  UNPROTECT(1); /* VR */

  Free(V);
  Free(Vold);
  Free(vv);
  Free(tmp);

  return VR;
}



/*
  FUNCTION: qp_fast_ipf_step
  PURPOSE: implement the Iterative Proportional Fitting (IPF) algorithm. Part of the
           R code below has been borrowed from an implementation by Graham Wills in
           June of 1992
  PARAMETERS: n - number of rows/columns of the input matrix
              Vf -
              Vn - (Vf and Vn have the same dimensions)
              a - clique
              csze - clique size
  RETURNS: the input matrix adjusted to the constraints of the clique
*/

static void
fast_ipf_step(int n, double* Vf, double* Vn, int* a, int csze) {
  int*    b;     /* vertices the clique a */
  double* Vfaa;
  double* Vnaa;
  double* Vni;
  double* Vnab;
  double* Vnba;
  double* Vnbb;
  double* Bnba;
  double* Vnbba;
  double* tmp1,*tmp2,*tmp3,*tmp4;
  int     i,j;

  b = Calloc(n-csze,int);

  /* b <- (1:length(Vf[, 1]))[ - a] */
  setdiff(n,csze,a,b);

  Vfaa  = Calloc(csze*csze,double);
  Vnaa  = Calloc(csze*csze,double);
  Vni   = Calloc(csze*csze,double);

  /* Vfaa <- Vf[a, a] */
  for (i=0;i<csze;i++)
    for (j=0;j<csze;j++) {
      Vfaa[i+j*csze] = Vf[a[i]+a[j]*n];
      Vni[i+j*csze] = i == j ? 1.0 : 0.0;
    }

  /* Vni <- solve(Vn[a, a]) */
  matsubm(Vnaa,Vn,n,a,csze,a,csze);
  matinv(Vni,Vnaa,csze,0);

  /* Bnba <- Vn[b, a] %*% Vni */
  Vnba = Calloc((n-csze)*csze,double);
  matsubm(Vnba,Vn,n,b,n-csze,a,csze);
  Bnba = Calloc((n-csze)*csze,double);
  matprod(Vnba,n-csze,csze,Vni,csze,csze,Bnba);

  /* Vnbba <- Vn[b, b] - Vn[b, a] %*% Vni %*% Vn[a, b] */
  Vnbb = Calloc((n-csze)*(n-csze),double);
  matsubm(Vnbb,Vn,n,b,n-csze,b,n-csze);
  Vnab = Calloc(csze*(n-csze),double);
  matsubm(Vnab,Vn,n,a,csze,b,n-csze);
  tmp1 = Calloc(csze*(n-csze),double);
  matprod(Vni,csze,csze,Vnab,csze,n-csze,tmp1);
  tmp2 = Calloc((n-csze)*(n-csze),double);
  matprod(Vnba,n-csze,csze,tmp1,csze,n-csze,tmp2);
  Vnbba = Calloc((n-csze)*(n-csze),double);
  matsumf(Vnbba,n-csze,n-csze,Vnbb,tmp2,-1.0);
  Free(tmp1);
  Free(tmp2);

  /* V <- Vf
  for (i=0;i<n;i++)
    for (j=0;j<ncol;j++)
      Vn[i+n*j] = Vf[i+n*j]; */
  Memcpy(Vn,Vf,(size_t) n*n);

  /* V[b, a] <- Bnba %*% Vfaa */
  tmp1 = Calloc((n-csze)*csze,double);
  matprod(Bnba,n-csze,csze,Vfaa,csze,csze,tmp1); 
  matrepm(Vn,n,b,n-csze,a,csze,tmp1);
  Free(tmp1);

  /* V[a, b] <- t(V[b, a]) */
  matsubm(Vnba,Vn,n,b,n-csze,a,csze);
  mattran(Vnab,Vnba,n-csze,csze);
  matrepm(Vn,n,a,csze,b,n-csze,Vnab);

  /* V[b, b] <- Vnbba + Bnba %*% Vfaa %*% t(Bnba) */
  tmp1 = Calloc(csze*(n-csze),double);
  mattran(tmp1,Bnba,n-csze,csze);
  tmp2 = Calloc(csze*(n-csze),double);
  matprod(Vfaa,csze,csze,tmp1,csze,n-csze,tmp2);
  Free(tmp1);
  tmp3 = Calloc((n-csze)*(n-csze),double);
  matprod(Bnba,n-csze,csze,tmp2,csze,n-csze,tmp3);
  Free(tmp2);
  tmp4 = Calloc((n-csze)*(n-csze),double);
  matsumf(tmp4,n-csze,n-csze,Vnbba,tmp3,1.0);
  Free(tmp3);
  matrepm(Vn,n,b,n-csze,b,n-csze,tmp4);
  Free(tmp4);

  Free(Vfaa);
  Free(Vnaa);
  Free(Vni);
  Free(Vnab);
  Free(Vnba);
  Free(Vnbb);
  Free(Bnba);
  Free(Vnbba);
  Free(b);
}



/*
  FUNCTION: qp_fast_htf
  PURPOSE: implement the algorithm of Hastie, Tibshirani and Friedman to perform
           maximum likelihood estimation of the sample covariance matrix given the
           independence constraints from an input undirected graph. Part of the
           R code below has been borrowed from an implementation by Giovanni Marchetti
           in the 'ggm' package (thanks Giovanni!!)
  PARAMETERS: SR - sample variance-covariance matrix
              AR - adjacency matrix of the undirected graph encoding independence constraints
              tolR - tolerance under which the main loop stops
              verbose - when set to TRUE the algorithm shows progression
  RETURNS: the input matrix adjusted to the constraints of the list of cliques
*/

static SEXP
qp_fast_htf(SEXP SR, SEXP AR, SEXP tolR, SEXP verbose) {
  int           n_var = INTEGER(getAttrib(SR,R_DimSymbol))[0];
  int*          vtc_wo_i;
  int*          A;
  int*          Ai;
  double*       W11;
  double*       W11Ai;
  double*       S;
  double*       W;
  double*       w;
  double*       s12Ai;
  double*       beta;
  double*       Wprev;
  double*       tmp;
  double        tol;
  double        precision;
  double        prev_precision;
  SEXP          WR;
  int           pct,ppct;
  int           i,j,k;

  S   = REAL(SR);
  A   = INTEGER(AR);
  tol = REAL(tolR)[0];

  vtc_wo_i = Calloc(n_var, int);
  Ai       = Calloc(n_var, int);
  beta     = Calloc(n_var-1, double);
  w        = Calloc(n_var-1, double);
  W11      = Calloc((n_var-1)*(n_var-1), double);
  s12Ai    = Calloc(n_var-1, double);
  Wprev    = Calloc(n_var*n_var, double);
  tmp      = Calloc(n_var*n_var, double);

  for (i=0; i < n_var; i++)
    vtc_wo_i[i] = i;

  PROTECT(WR = allocMatrix(REALSXP, n_var, n_var));
  W = REAL(WR);
  /* W <- Wprev <- S */
  Memcpy(W, REAL(SR), (size_t) (n_var * n_var));

  precision = DBL_MAX;
  while (precision > tol) {
    prev_precision = precision;
    Memcpy(Wprev, W, (size_t) (n_var * n_var));

    ppct = -1;
    for (i=0; i < n_var; i++) { /* i is the vertex currently adjusted */

      vtc_wo_i[i] = n_var-1;

      /* W11 <- W[-i, -i, drop=FALSE] */
      matsubm(W11, W, n_var, vtc_wo_i, n_var-1, vtc_wo_i, n_var-1);

      /* s12 <- S[-i, i, drop=FALSE] */
      /* Ai <- A[i, ] */
      /* Ai <- Ai[-i] */
      k = 0;
      for (j=0; j < n_var-1; j++) {
        if (A[i*n_var+vtc_wo_i[j]]) {
          Ai[k] = j;
          s12Ai[k] = S[i*n_var + vtc_wo_i[j]];
          k++;
        }
      }

      /* beta <- rep(0, n.var-1) */
      memset(beta, 0, (n_var-1) * sizeof(double));

      /* beta[Ai] <- solve(W11[Ai, Ai], s12[Ai, ]) */
      if (k > 0) {
        W11Ai = Calloc(k*k, double);
        matsubm(W11Ai, W11, n_var-1, Ai, k, Ai, k);
        matinv(s12Ai, W11Ai, k, 1);
        Free(W11Ai);
        for (j=0; j < k; j++)
          beta[Ai[j]] = s12Ai[j];
      }

      /* w <- W11 %*% beta */
      matprod(W11, n_var-1, n_var-1, beta, n_var-1, 1, w);

      /* W[-i, i] <- W[i, -i] <- w */
      for (j=0; j < n_var-1; j++)
        W[i*n_var + vtc_wo_i[j]] = W[vtc_wo_i[j]*n_var + i] = w[j];
      vtc_wo_i[i] = i;

      pct = (int) (((double) (i*100))/((double) n_var));
      if (pct != ppct) {
        if (INTEGER(verbose)[0]) {
          if (pct % 10 == 0)
            Rprintf("%d",pct);
          else
            Rprintf(".",pct);
          R_FlushConsole();
        }

        R_CheckUserInterrupt();
#ifdef Win32
        R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
        R_ProcessEvents();
#endif
        ppct = pct;
      }
    }
    if (INTEGER(verbose)[0])
      Rprintf("\n");

    matsumf(tmp, n_var, n_var, W, Wprev,-1.0);
    precision = matmxab(tmp, n_var, n_var);
    if (INTEGER(verbose)[0])
      Rprintf("Precision: %.10f\n", precision);
    if (precision > prev_precision)
      error("HTF is not converging, probably the input sample covariance matrix is calculated from too few observations\n");
  }

  Free(vtc_wo_i);
  Free(Ai);
  Free(beta);
  Free(w);
  Free(W11);
  Free(s12Ai);
  Free(Wprev);
  Free(tmp);

  UNPROTECT(1); /* WR */

  return WR;
}



/*
  FUNCTION: cov2cor
  PURPOSE: this is a C implementation of the cov2cor function from the stats package
  RETURNS: a scaled covariance matrix so that the diagonal equals unity
*/

static void
cov2cor(double* R, double* V, int n) {
  double* Is;
  int     i,j;

  Is = Calloc(n, double);
  for (i=0; i < n; i++)
    Is[i] = sqrt(1.0 / V[i + i*n]);

  for (i=0; i < n; i++) {
    for (j=0; j < i; j++)
      R[i + j*n] = R[j + i*n] = Is[i] * V[i + j*n] * Is[j];      
    R[i + i*n] = 1.0;
  }

  Free(Is);
}



/*
  FUNCTION: matprod
  PURPOSE: multiply two matrices by using the LaPACK library that
           comes along with the R distribution, this code is taken from

           R-2.2.0/src/main/array.c
  RETURNS: none
*/

static void
matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z) {
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0, sum;
    Rboolean have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
        /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
         * The test is only O(n) here
         */
        for (i = 0; i < nrx*ncx; i++)
            if (ISNAN(x[i])) {have_na = TRUE; break;}
        if (!have_na)
            for (i = 0; i < nry*ncy; i++)
                if (ISNAN(y[i])) {have_na = TRUE; break;}
        if (have_na) {
            for (i = 0; i < nrx; i++)
                for (k = 0; k < ncy; k++) {
                    sum = 0.0;
                    for (j = 0; j < ncx; j++)
                        sum += x[i + j * nrx] * y[j + k * nry];
                    z[i + k * nrx] = sum;
                }
        } else
            F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
                            x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
        for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}



/*
  FUNCTION: matinv
  PURPOSE: calculates de inverse of a matrix by using the LaPACK library
           that comes along with the R distribution, this code is taken from
           the function modLa_dgesv in file

           R-2.2.0/src/modules/lapack/Lapack.c
  RETURNS: none
*/

static void
matinv(double* inv, double* M, int n, int p) {
  int     i;
  int     info;
  int*    ipiv;
  double* avals;
  double* work;
  double  anorm;
  double  rcond;
  double  tol = DBL_MIN;

  if (p == 0) {
    memset(inv, 0, n * n * sizeof(double));
    for (i=0; i < n; i++)
      inv[i+i*n] = 1.0;
    p = n;
  }

  ipiv = (int *) Calloc(n, int);
  avals = (double *) Calloc(n*n, double);
  Memcpy(avals, M, (size_t) (n*n));

  F77_CALL(dgesv)(&n, &p, avals, &n, ipiv, inv, &n, &info);
  if (info < 0)
    error("argument %d of Lapack routine %s had invalid value",-info, "dgesv");
  if (info > 0)
    error("Lapack routine dgesv: system is exactly singular");

  anorm = F77_CALL(dlange)("1", &n, &n, M, &n, (double*) NULL);

  work = (double *) Calloc(4*n, double);

  F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info);
  if (rcond < tol)
    error("system is computationally singular: reciprocal condition number = %g",rcond);

  Free(ipiv);
  Free(avals);
  Free(work);
}



/*
  FUNCTION: symmatlogdet
  PURPOSE: calculates de determinant of a symmetric matrix by using the LaPACK
           library that comes along with the R distribution, this code is taken
           from the function moddet_ge_real in file

           R-2.13.0/src/modules/lapack/Lapack.c

           the input symmetric matrix should come as the upper triangle
  RETURNS: none
*/

static double
symmatlogdet(double* M, int n, int* sign) {
  double* A;
  double  modulus = 0.0;
  int     i,j;
  int*    jpvt;
  int     info;

  /* blowing up the memory footprint, it would be nice to
   * make these calculations directly on a symmetric matrix */
  A = Calloc(n*n, double);
  for (i=0; i < n; i++)
    for (j=0; j <= i; j++)
      A[i + j*n] = A[j + i*n] = M[UTE2I(i, j)];

  jpvt = Calloc(n, int);
  F77_CALL(dgetrf)(&n, &n, A, &n, jpvt, &info);

  *sign = 1;
  if (info < 0)
    error("error code %d from Lapack routine '%s'", info, "dgetrf");
  if (info > 0)
    warning("Lapack routine dgetrf: system is exactly singular");

  if (info == 0) {
    for (i=0; i < n; i++)
      if (jpvt[i] != (i + 1))
        *sign = -(*sign);

    for (i=0; i < n; i++) {
      double dii = A[i * (n + 1)];

      modulus += log(dii < 0 ? -dii : dii);
      if (dii < 0)
        *sign = -(*sign);
    }
  }

  Free(jpvt);
  Free(A);

  return modulus;
}



/*
  FUNCTION: matsumf
  PURPOSE: calculates the sum of two matrices multiplying the second one by a
           scalar factor
  RETURNS: none
*/

static void
matsumf(double* R, int nrow, int ncol, double* M, double* N, double factor) {
  int  i,j;

  for (i=0;i<nrow;i++)
    for (j=0;j<ncol;j++)
      R[i+nrow*j] = M[i+nrow*j] + factor * N[i+nrow*j];
}



/*
  FUNCTION: matscalarprod
  PURPOSE: calculates the scalar product of two matrices
  RETURNS: none
*/

static void
matscalarprod(double* R, int nrow, int ncol, double* M, double* N) {
  int  i,j;

  for (i=0;i<nrow;i++)
    for (j=0;j<ncol;j++)
      R[i+nrow*j] = M[i+nrow*j] * N[i+nrow*j];
}



/*
  FUNCTION: mattran
  PURPOSE: performs the transposition of a matrix
  RETURNS: none
*/

static void
mattran(double* T, double* M, int nrow, int ncol) {
  int  i,j;

  for (i=0;i<nrow;i++)
    for (j=0;j<ncol;j++)
      T[j+ncol*i] = M[i+nrow*j];
}



/*
  FUNCTION: matsubm
  PURPOSE: extracts a square submatrix from a square matrix
  RETURNS: none
*/

static void
matsubm(double* subM, double* M, int n, int* subrows, int nsubrows,
        int* subcols, int nsubcols) {
  int  i,j;

  for (i=0; i < nsubrows; i++)
    for (j=0; j < nsubcols; j++)
      subM[i + nsubrows*j] = M[subrows[i] + n*subcols[j]];
}



/*
  FUNCTION: symmatsubm
  PURPOSE: extracts a symmetric submatrix of a symmetric matrix
  RETURNS: none
*/

static void
symmatsubm(double* subM, double* M, int n, int* subrows, int nsubrows) {
  int  i,j;

  for (i=0; i < nsubrows; i++)
    for (j=0;j <= i; j++)
      subM[UTE2I(i,j)] = M[UTE2I(subrows[i], subrows[j])];
}



/*
  FUNCTION: matmxab
  PURPOSE: find the absolut maximum value in a matrix of doubles
  RETURNS: none
*/

static double
matmxab(double* M, int nrow, int ncol) {
  int    i,j;
  double maxabs = 0.0;

  for (i=0;i<nrow;i++)
    for (j=0;j<ncol;j++) {
      /* double abs = M[i+nrow*j] > 0 ? M[i+nrow*j] : -1.0*M[i+nrow*j]; */
      double abs = fabs(M[i+nrow*j]);
      maxabs = maxabs < abs ? abs : maxabs;
    }

  return maxabs;
}



/*
  FUNCTION: matrepm
  PURPOSE: replace a submatrix with another matrix
  RETURNS: none
*/

static void
matrepm(double* M, int n, int* subrows, int nsubrows,
        int* subcols, int nsubcols, double* N) {
  int  i,j;

  for (i=0;i<nsubrows;i++)
    for (j=0;j<nsubcols;j++)
      M[subrows[i]+n*subcols[j]] = N[i+nsubrows*j];
}



/*
  FUNCTION: matrepm
  PURPOSE: puts in b all integers from 0 till n-1 that are not
           in a where a is a vector of integers of size m
  RETURNS: none
*/

static void
setdiff(int n, int m, int* a, int* b) {
  int i;
  int k = 0;

  for (i=0;i<n;i++) {
    int j=0;
    while (j < m && a[j] != i)
      j++;

    if (j == m && k<n-m)
      b[k++] = i;
  }

}



/*
  FUNCTION: i2e
  PURPOSE: transform a non-negative integer representing an edge into two
           non-negative integers representing the vertices joined by this edge
  PARAMETERS: i - non-negative integer representing an edge
              e_i - vertex joined by the edge
              e_j - vertex joined by the edge, always strictly smaller than e_i
  RETURN: none
*/

void
i2e(int i, int* e_i, int* e_j) {
  *e_i = 1 + (unsigned int) (-0.5 + sqrt(0.25 + 2.0 * ((double) i)));
  *e_j = i - (unsigned int) ((double) ((*e_i)*((*e_i)-1)) / 2.0);
}



/*
  FUNCTION: e2i
  PURPOSE: transform two non-negative integers representing the vertices of an edge
           into a non-negative integer representing this same edge
  PARAMETERS: e_i - vertex joined by the edge
              e_j - vertex joined by the edge, always strictly smaller than e_i
              i - non-negative integer representing an edge
  RETURN: none
*/

int
e2i(int e_i, int e_j, int* i) {
  if (e_i < e_j) { /* e_j should always be smaller than e_i */
    e_i = e_i ^ e_j;
    e_j = e_i ^ e_j;
    e_i = e_i ^ e_j;
  }

  return(((int) (((double) (e_i * (e_i - 1))) / 2.0)) + e_j);
}



/*
  FUNCTION: missing_obs
  PURPOSE: check whether there are missing obserations
  PARAMETERS: X - vector containing the column-major stored matrix of values
              p - number of variables in X
              n - number of observations in X
              Y - vector containing the indices of the variables in X for which we
                  want to find missing values
              n_Y - number of elements in Y
              idx_obs - indices of the observations for which we want to check missingness
              n_idx_obs - number of observations for which we want to check missingness
  DESCRIPTION: this function serves the purpose of checking whether there are missing observations
               throughout X, or throughout a subset of variables in Y, or throughout a subset
               of observations in idx_obs, or both.
  RETURN: TRUE if there are missing observations; FALSE otherwise
*/

int
missing_obs(double* X, int p, int n, int* Y, int n_Y, int* idx_obs, int n_idx_obs) {
  int i,j,k,l;
  int missing_flag = 0;

  i = 0;
  while (i < n_idx_obs && !missing_flag) {
    k = n_idx_obs < n ? idx_obs[i] : i;
    j = 0;
    while (j < n_Y && !missing_flag) {
      l = n_Y < p ? Y[j] : j;
      if (ISNA(X[l * n + k]))
        missing_flag = 1;
      j++;
    }

    i++;
  }

  return missing_flag;
}



/*
  FUNCTION: find_missing_obs
  PURPOSE: create a logical mask indicating observations with at least one missing value
  PARAMETERS: X - vector containing the column-major stored matrix of values
              p - number of variables in X
              n - number of observations in X
              Y - vector containing the indices of the variables in X for which we
                  want to find missing values
              n_Y - number of elements in Y
              idx_obs - indices of the observations for which we want to check missingness
              n_idx_obs - number of observations for which we want to check missingness
              missing_mask - logical mask where this function sets to 1 if a value is missing
                             in an observation. it is assumed that initially is set to zeroes
                             on positions where observations are not missing.
  DESCRIPTION: this function serves the purpose of identifying missing observations either all
               throughout X, or throughout a subset of variables in Y, or throughout a subset
               of observations in idx_obs, or both.
  RETURN: the number of missing observations and, by reference, a logical mask (missing_mask)
           with positions corresponding to missing observations set to 1
*/

int
find_missing_obs(double* X, int p, int n, int* Y, int n_Y, int* idx_obs,
                 int n_idx_obs, int* missing_mask) {
  int i,j,k,l;
  int n_mis = 0;

  for (i=0; i < n_idx_obs; i++) {
    k = n_idx_obs < n ? idx_obs[i] : i;
    j = 0;
    while (!missing_mask[k] && j < n_Y) {
      l = n_Y < p ? Y[j] : j;
      if (ISNA(X[l * n + k]))
        missing_mask[k] = 1;
      j++;
    }

    if (missing_mask[k])
      n_mis++;
  }

  return n_mis;
}



/*
  FUNCTION: calculate_means
  PURPOSE: calculate the means of the values at the columns of the input matrix
           provided as a column-major vector
  PARAMETERS: X - vector containing the column-major stored matrix of values
              p - number of variables
              n - number of observations
              Y - vector containing the indices of the variables in X for which we
                  want to calculate the mean
              n_Y - number of elements in Y
              idx_obs - indices of the observations to employ in the calculation of the mean
              n_idx_obs - number of observations to employ in the calculation of the mean
              meanv - output vector of n_Y mean values
              missing_mask - logical mask indicating what observations contain at least one
                             missing value
              n_mis - number of missing observations
  RETURN: none
*/

void
calculate_means(double* X, int p, int n, int* Y, int n_Y, int* idx_obs,
                int n_idx_obs, int* missing_mask, int n_mis, double* meanv) {
  long double sum, tmp;
  double*     xx;
  int         i,j;

  for (i=0;i < n_Y;i++) {
    xx = n_Y < p ? &X[Y[i] * n] : &X[i * n];
    sum = 0.0;
    for (j=0;j < n_idx_obs;j++)
      if (!missing_mask[n_idx_obs < n ? idx_obs[j] : j])
        sum += n_idx_obs < n ? xx[idx_obs[j]] : xx[j];
    tmp = sum / (n_idx_obs-n_mis);
    if (R_FINITE((double) tmp)) {
      sum = 0.0;
      for (j=0;j < n_idx_obs;j++) {
        if (!missing_mask[n_idx_obs < n ? idx_obs[j] : j])
          sum += n_idx_obs < n ? (xx[idx_obs[j]] - tmp) : (xx[j] - tmp);
      }
      tmp = tmp + sum / (n_idx_obs-n_mis);
    }
    meanv[i] = tmp;
  }
}



/*
  FUNCTION: ssd
  PURPOSE: calculate the, corrected or not, sum of squares and deviations matrix
           returning only the upper triangle of the matrix in column-major order
           (for creating later a dspMatrix object). In the presence of missing data
           this function uses complete observations only and the returned value
           corresponds to its number.
  PARAMETERS: X - vector containing the column-major stored matrix of values
              p - number of variables
              n - number of observations
              Y - vector containing the indices of the variables in X for which we
                  want to calculate the ssd
              n_Y - number of elements in Y, it should equal p when Y is NULL
              idx_obs - indices of the observations to employ in the calculation of the ssd
              n_idx_obs - number of observations to employ in the calculation of the ssd, it
                          should be n when idx_obs is NULL
              corrected - flag indicating whether the sum be corrected (=covariance)
              missing_mask - logical mask indicating which observations are missing
              ssd_mat - pointer to the matrix where the result is returned
  DESCRIPTION: this function skips missing observations and can work cooperatively with
               ssd_A() to push through 'missing_mask' the detected missing observations
  RETURN: number of complete observations employed in the estimation
*/

int
ssd(double* X, int p, int n, int* Y, int n_Y, int* idx_obs, int n_idx_obs,
    int corrected, int* missing_mask, double* ssd_mat) {
  double* meanv;
  int     allocated_missing_mask = 0;
  int     n1, n_mis=0;
  int     i,j,k,l;

  meanv = Calloc(n_Y, double);
  if (missing_mask == NULL) {
    missing_mask = Calloc(n, int); /* assume Calloc() sets memory to zeroes */
    allocated_missing_mask = 1;
  }

  n_mis = find_missing_obs(X, p, n, Y, n_Y, idx_obs, n_idx_obs, missing_mask);

  calculate_means(X, p, n, Y, n_Y, idx_obs, n_idx_obs, missing_mask, n_mis, meanv);

  n1 = n_idx_obs - n_mis - 1;

  if (corrected && n1 < 1)
    error("not enough complete observations available to calculate a corrected SSD matrix (n-1=%d, n_obs=%d, n_mis=%d)\n", n1, n_idx_obs, n_mis);

  l = 0;
  for (i=0; i < n_Y; i++)
    for (j=0; j <= i; j++) {
      double*     xx;
      double*     yy;
      long double xxm, yym, sum;

      xx  = n_Y < p ? &X[Y[i] * n] : &X[i * n];
      xxm = meanv[i];
      yy  = n_Y < p ? &X[Y[j] * n] : &X[j * n];
      yym = meanv[j];
      sum = 0.0;
      for (k=0; k < n_idx_obs; k++)
        if (n_mis == 0 || (!missing_mask[n_idx_obs < n ? idx_obs[k] : k]))
          sum += n_idx_obs < n ? (xx[idx_obs[k]] - xxm) * (yy[idx_obs[k]] - yym) : (xx[k] - xxm) * (yy[k] - yym);

      /* here we assume that ssd_mat was properly initialized to zeroes before the call
         to this function, not because we could not memset() here but because of the
         convinience of summing through the output argument ssd_mat when calculating an
         ssd matrix from mixed data */
      ssd_mat[l] += corrected ? (double) (sum / ((long double) n1)) : (double) sum;
      l++;
    }

  if (allocated_missing_mask)
    Free(missing_mask);

  Free(meanv);

  return n_idx_obs - n_mis; /* return number of complete observations */
}



/*
  FUNCTION: qp_fast_cov_upper_triangular
  PURPOSE: calculate a covariance matrix returning only the upper triangle
           of the matrix in column-major order (for creating later a dspMatrix object)
  PARAMETERS: XR - matrix of multivariate observations
              corrected - flag set to TRUE when calculating the sample covariance matrix
                          and set to FALSE when calculating the uncorrected sum of
                          squares and deviations
  RETURN: an SsdMatrix object containing the (corrected) sum of squares and deviations and
          the number of complete observations employed in the calculations
*/

static SEXP
qp_fast_cov_upper_triangular(SEXP XR, SEXP corrected) {
  SEXP          ssdR, dimX, dimNamesX, dimNamesSsd;
  double*       X;
  double*       ssd_val;
  int*          dim_ssd;
  int           p, n, n_obs, n_upper_tri;
  PROTECT_INDEX Xpi;

  PROTECT_WITH_INDEX(XR,&Xpi);

  REPROTECT(XR  = coerceVector(XR,REALSXP),Xpi);
  X = REAL(XR);
  dimX = getAttrib(XR, R_DimSymbol);
  dimNamesX = getAttrib(XR, R_DimNamesSymbol);

  /* number of observations equals number of rows */
  n = INTEGER(dimX)[0];
  /* number of variables equals number of columns */
  p = INTEGER(dimX)[1];
  /* number of elements in the upper triangular matrix including diagonal */
  n_upper_tri = (p * (p + 1)) / 2;

  ssdR = PROTECT(NEW_OBJECT(MAKE_CLASS("SsdMatrix")));

  SET_SLOT(ssdR, SsdMatrix_ssdSym, NEW_OBJECT(MAKE_CLASS("dspMatrix")));

  dim_ssd = INTEGER(ALLOC_SLOT(GET_SLOT(ssdR, SsdMatrix_ssdSym), Matrix_DimSym, INTSXP, 2));
  dim_ssd[0] = dim_ssd[1] = INTEGER(dimX)[1];

  if (dimNamesX != R_NilValue) {
    dimNamesSsd = ALLOC_SLOT(GET_SLOT(ssdR, SsdMatrix_ssdSym), Matrix_DimNamesSym, VECSXP, 2);
    SET_VECTOR_ELT(dimNamesSsd, 0, duplicate(VECTOR_ELT(dimNamesX, 1)));
    SET_VECTOR_ELT(dimNamesSsd, 1, duplicate(VECTOR_ELT(dimNamesX, 1)));
  }

  ssd_val = REAL(ALLOC_SLOT(GET_SLOT(ssdR, SsdMatrix_ssdSym), Matrix_xSym, REALSXP, n_upper_tri));
  SET_SLOT(GET_SLOT(ssdR, SsdMatrix_ssdSym), Matrix_uploSym, mkString("U"));

  /* the following line is needed before each call to ssd() where the last argument points
   * to a memory chunk allocated with anything that does not set the memory to zeroes */
  memset(ssd_val, 0, sizeof(double) * n_upper_tri);

  n_obs = ssd(X, p, n, NULL, p, NULL, n, INTEGER(corrected)[0], NULL, ssd_val);

  SET_SLOT(ssdR, SsdMatrix_nSym, ScalarInteger(n_obs));

  UNPROTECT(2); /* XR ssdR */

  return(ssdR);
}



/*
  FUNCTION: calculate_xtab
  PURPOSE: calculate the cross-tabulation of the observations according to the
           discrete variables indicated in the arguments. The discrete values
           should be encoded as natural numbers 1, 2, ... just as the encoding
           of factor variables in R. Missing values encoded as NAs are treated
           by setting to -1 the cross-classified observation.
  IMPORTANT: it assumes the argument xtab has all its positions set to the value of 1 <- SHOULD BE 0 ??
  PARAMETERS: X - vector containing the column-major stored matrix of values
              p - number of variables
              n - number of values per variable
              xtab - pointer where the result is returned
  RETURN: none
*/

void
calculate_xtab(double* X, int p, int n, int* I, int n_I, int* n_levels, int* xtab) {
  int i,j;
  int base=1;

  /*
  for (i=0; i < n; i++)
    xtab[i] = 1;
  */

  for (i=0; i < n_I; i++) {
    for (j=0; j < n; j++) {
      if (xtab[j] > 0) {
        double level = X[j + I[i] * n];

        if (ISNA(level))
          xtab[j] = -1;
        else {
          if (level <= 0 || (level-((double) ((int) level))) > 0)
            error("observation %d contains discrete levels that are not positive integers\n", j+1);
          xtab[j] = xtab[j] + base * ((int) (level-1.0));
        }
      }
    }
    base = base * n_levels[I[i]]; /* WAS n_levels[i] */
  }

  return;
}



/*
  FUNCTION: ssd_A
  PURPOSE: calculate the uncorrected sum of squares and deviations matrix
           for a subset of variables A = I \cap Y returning only the upper
           triangle of the matrix in column-major order (for creating later a dspMatrix object)
  PARAMETERS: X - vector containing the column-major stored matrix of values
              p - number of variables
              n - number of values per variable
              I - vector containing the indices of the discrete variables in X which we
                  want to employ in the calculation of the ssd
              n_I - number of elements in I
              Y - vector containing the indices of the variables in X for which we
                  want to calculate the ssd
              n_Y - number of elements in Y
              excobs_mask - logical mask of observations that should be excluded from the calculations
              missing_mask - (output) logical mask of observations that contain missing values
              ssd_A - (output) pointer to the matrix where the result is returned
  RETURN: none
*/

int
ssd_A(double* X, int p, int n, int* I, int n_I, int* n_levels, int* Y, int n_Y,
      int* excobs_mask, int* missing_mask, double* ssd_A) {
  int*    obs_idx;
  int     n_obs, n_co;
  int     i,j,k,m;

  obs_idx = Calloc(n, int);
  global_xtab = Calloc(n, int);
  n_obs = 0;
  for (i=0; i < n; i++) {
    global_xtab[i] = 1;
    if (excobs_mask != NULL) {
      if (!excobs_mask[i])     /* if obs is not excluded,  use it */
        obs_idx[n_obs++] = i;
      else                     /* if obs is excluded, avoid using it later */
        global_xtab[i] = -1;
    } else
      obs_idx[n_obs++] = i;
  }

  if (n_I == 0) {
    n_co = ssd(X, p, n, Y, n_Y, obs_idx, n_obs, FALSE, missing_mask, ssd_A);

    Free(obs_idx);
    Free(global_xtab);

    return n_co;
  }

  calculate_xtab(X, p, n, I, n_I, n_levels, global_xtab);

  /* group together observations from the same joint discrete level putting the *
   * observations from missing discrete levels at the beginning                 */
  qsort(obs_idx, n_obs, sizeof(int), indirect_int_cmp);

  /* skip missing discrete observations */
  i = 0;
  while (i < n_obs && global_xtab[obs_idx[i]] == -1) {
    if (missing_mask != NULL)
      missing_mask[obs_idx[i]] = 1;
    i++;
  }

  n_co = 0;

  while (i < n_obs) {
    j = i;
    while (j < n_obs && global_xtab[obs_idx[i]] == global_xtab[obs_idx[j]])
      j++;

    n_co += ssd(X, p, n, Y, n_Y, obs_idx+i, j-i, FALSE, missing_mask, ssd_A);

    i = j;
  }

  Free(obs_idx);
  Free(global_xtab);

  return n_co;
}



/*
  FUNCTION: qp_fast_rnd_graph
  PURPOSE: samples a d-regular graph uniformly at random
  PARAMETERS: pR - number of vertices
              dR - vertex degree
              excludeR - vertex whose pairwise adjacencies should be excluded
  RETURN: upper triangle in column major of the adjacency matrix of the sampled graph
*/

typedef struct {
  double x;
  int ix, jx;
} DblWithIdx;

int
dbl_cmp_desc_idx_decr(const void *a, const void *b) {
  const DblWithIdx* ia = (const DblWithIdx *) a;
  const DblWithIdx* ib = (const DblWithIdx *) b;

  return ib->x - ia->x;
}

int
int_cmp_desc_idx_incr(const void *a, const void *b) {
  const IntWithIdx* ia = (const IntWithIdx *) a;
  const IntWithIdx* ib = (const IntWithIdx *) b;

  return ia->x - ib->x;
}

static SEXP
qp_fast_rnd_graph(SEXP pR, SEXP dR, SEXP excludeR, SEXP verboseR) {
  int         p=INTEGER(pR)[0];
  int         d=INTEGER(dR)[0];
  int         verbose=INTEGER(verboseR)[0];
  int*        exclude;
  SEXP        GR;
  SEXP        pb=R_NilValue;
  Rboolean*   G;
  double*     deg_diff;
  IntWithIdx* deg;
  IntWithIdx* working_deg;
  Rboolean    regular;
  int         n_upper_tri;
  int         i,j,k;

  /* number of elements in the upper triangular matrix including diagonal */
  n_upper_tri = (p * (p+1)) / 2;
  
  deg         = Calloc(p, IntWithIdx);
  working_deg = Calloc(p, IntWithIdx);
  deg_diff    = Calloc(p, double);
  exclude     = Calloc(p, int); /* assume Calloc() initializes memory to 0 */

  if (verbose) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("txtProgressBar")); t=CDR(t);
    SETCAR(t, ScalarInteger(3));
    SET_TAG(t, install("style"));
    SETCAR(t, ScalarInteger(0));
    SET_TAG(t, install("min"));
    SETCAR(t, ScalarInteger(p));
    SET_TAG(t, install("max"));
    PROTECT(pb = eval(s, R_GlobalEnv));
    UNPROTECT(1); /* t s */
  }

  PROTECT(GR = allocVector(LGLSXP, n_upper_tri)); /* upper triangle includes diagonal */
  G = LOGICAL(GR);

  if (excludeR != R_NilValue) {
    for (i=0; i < length(excludeR); i++)
      exclude[INTEGER(excludeR)[i]-1] = 1;
  }

  regular=FALSE;
  while (!regular) {
    int n_vtx_left=p;
    int  missing_edges=1;

    for (i=0; i < p; i++) {
      deg[i].x = 0;
      deg[i].ix = i;
    } 
    memset(G, FALSE, sizeof(Rboolean) * n_upper_tri);

    R_CheckUserInterrupt();
#ifdef Win32
    R_ProcessEvents();
#endif
#ifdef HAVE_AQUA
    R_ProcessEvents();
#endif

    while (n_vtx_left > 0 && missing_edges > 0) {
      Memcpy(working_deg, deg, (size_t) p);
      qsort(working_deg, p, sizeof(IntWithIdx), int_cmp_desc_idx_incr);
      n_vtx_left=0;
      while (n_vtx_left < p && working_deg[n_vtx_left].x < d) {
        deg_diff[n_vtx_left] = (double) (d - working_deg[n_vtx_left].x);
        n_vtx_left++;
      }

      if (verbose) {
        SEXP s, t;
        PROTECT(t = s = allocList(3));
        SET_TYPEOF(s, LANGSXP);
        SETCAR(t, install("setTxtProgressBar")); t=CDR(t);
        SETCAR(t, pb);
        SET_TAG(t, install("pb")); t=CDR(t);
        SETCAR(t, ScalarInteger(p-n_vtx_left));
        SET_TAG(t, install("value"));
        eval(s, R_GlobalEnv);
        UNPROTECT(1); /* t s */
      }

      if (n_vtx_left > 0) {
        double* S;
        DblWithIdx* cdf;
        int     n_cdf;
        double  sum_S;

        S = Calloc(n_vtx_left * n_vtx_left, double);
        cdf = Calloc((n_vtx_left * (n_vtx_left-1))/2, DblWithIdx);

        /* outer product */
        matprod(deg_diff, n_vtx_left, 1, deg_diff, 1, n_vtx_left, S);

        missing_edges = 0;
        sum_S = 0.0;
        for (i=0; i < n_vtx_left; i++)
          for (j=i+1; j < n_vtx_left; j++)
            if (G[UTE2I(working_deg[i].ix, working_deg[j].ix)])
              S[i+j*n_vtx_left] = S[j+i*n_vtx_left] = -1; /* exclude adjacent pairs of vertices */
            else {
              if (!exclude[working_deg[i].ix] || !exclude[working_deg[j].ix]) {
                sum_S = sum_S + S[i+j*n_vtx_left];
                missing_edges++;
              }
            }

        if (missing_edges > 0) {
          n_cdf = 0;
          for (i=0; i < n_vtx_left; i++)
            for (j=i+1; j < n_vtx_left; j++)
              if (!G[UTE2I(working_deg[i].ix, working_deg[j].ix)] &&
                  (!exclude[working_deg[i].ix] || !exclude[working_deg[j].ix])) {
                cdf[E2I(i, j)].x = S[i+j*n_vtx_left] / sum_S;
                cdf[E2I(i, j)].ix = i;
                cdf[E2I(i, j)].jx = j;
                n_cdf++;
              } else {
                cdf[E2I(i, j)].x = -1;
                cdf[E2I(i, j)].ix = -1;
                cdf[E2I(i, j)].jx = -1;
              }

          if (n_cdf > 0) {
            double r;
            double cumsum = 0.0;

            qsort(cdf, (n_vtx_left*(n_vtx_left-1))/2, sizeof(DblWithIdx), dbl_cmp_desc_idx_decr);
            r = unif_rand(); /* sample using R-builtin RNG */
            k = -1;
            while (k < n_cdf && r > cumsum) {
              k++;
              cumsum = cumsum + cdf[k].x;
            }

            i = cdf[k].ix;
            j = cdf[k].jx;
            G[UTE2I(working_deg[i].ix , working_deg[j].ix)] = TRUE;
            deg[working_deg[i].ix].x++;
            deg[working_deg[j].ix].x++;
          } /* else
            error("No cdf could be built\n"); */
    
          Free(cdf);
          Free(S);
        } else
          n_vtx_left = 0;
      }
    }

    regular=TRUE;
    i=0;
    while (i < p && regular) {
      if (deg[i].x != d)
        regular=FALSE;
      i++;
    }
  }

  Free(exclude);
  Free(deg_diff);
  Free(working_deg);
  Free(deg);

  if (verbose) {
    SEXP s, t;
    PROTECT(t = s = allocList(2));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, install("close")); t=CDR(t);
    SETCAR(t, pb);
    eval(s, R_GlobalEnv);
    UNPROTECT(2); /* t s pb */
  }

  UNPROTECT(1); /* GR */

  return(GR);
}
