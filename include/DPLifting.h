/**
 * @file DPLifting.h
 * @brief DPL / DPT hybrid lifting for knapsack sets.
 *
 * Preferred public entry points:
 *   - lifting()              (C++)
 *   - dplifting_lift_cover() (C ABI)
 * Other Lifting_* symbols are internal helpers exposed for advanced use / tests.
 *
 * Version: see DPLIFTING_VERSION.
 */
#ifndef __DPLIFTING_H__
#define __DPLIFTING_H__

#define DPLIFTING_VERSION_MAJOR 1
#define DPLIFTING_VERSION_MINOR 4
#define DPLIFTING_VERSION_PATCH 1
#define DPLIFTING_VERSION "1.4.1"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#ifndef DPLIFTING_API
#  if defined(_WIN32) && defined(DPLIFTING_BUILD_SHARED)
#    ifdef DPLIFTING_EXPORTS
#      define DPLIFTING_API __declspec(dllexport)
#    else
#      define DPLIFTING_API __declspec(dllimport)
#    endif
#  else
#    define DPLIFTING_API
#  endif
#endif

#define EPS_DPL 1e-6
#define INF_DPL 1e+20
/** lifting(..., isdpl_mode) — DPL/DPT axis (orthogonal to +R):
 *   MODE_DPL / MODE_DPT     — hard force, no mid-lift
 *   MODE_THRESHOLD        — manual τ-policy: initial + mid-lift by threshold vs bar b^k
 *   MODE_AUTO             — default feature selector on (ρ_w,β,ū); ignores threshold for mode
 */
#define DPLIFTING_MODE_AUTO      (-1)
#define DPLIFTING_MODE_THRESHOLD (-2)
/** Force DPT table (isDPL=0); no mid-lift switching */
#define DPLIFTING_MODE_DPT    0
/** Force DPL table (isDPL=1); no mid-lift switching */
#define DPLIFTING_MODE_DPL    1

/** Default empirical thresholds (factorab): rho_th = beta_th = 6 */
#define DPLIFTING_DEFAULT_RHO_TH   6.0
#define DPLIFTING_DEFAULT_BETA_TH  6.0
#define DPLIFTING_DEFAULT_UBAR_TH  3.0
/** +R AUTO: disable when λ = b / Σ w_i u_i is below this (tight rows). */
#define DPLIFTING_DEFAULT_LAMBDA_TH 0.1

/**
 * +R axis (orthogonal to DPL/DPT), lift->reduction_request:
 *   RED_ON / RED_OFF — manual
 *   RED_AUTO         — default: enable iff λ ≥ λ_th and bar b^0 > τ
 */
#define DPLIFTING_RED_AUTO 0
#define DPLIFTING_RED_ON   1
#define DPLIFTING_RED_OFF  2
#define MIN_DPL(a,b) (a<=b? a:b) 
#define MAX_DPL(a,b) (a>=b? a:b) 
#define FLOOR_DPL(a) ( floor( a + EPS_DPL ) )
#define CEIL_DPL(a) ( ceil( a - EPS_DPL ) )
#define FLOOR_INT(a) ( (long)floor( a + EPS_DPL ) )
#define CEIL_INT(a) ( (long)ceil( a - EPS_DPL ) )

inline double ABS_DPL(double a)
{
   return fabs(a);
}

inline int ISLE(double a, double b)
{
   return a <= b + EPS_DPL;
}

inline int ISLT(double a, double b)
{
   return a < b - EPS_DPL;
}

inline int ISGE(double a, double b)
{
   return a + EPS_DPL >= b;
}

inline int ISGT(double a, double b)
{
   return a - EPS_DPL > b;
}

inline int ISZERO(double a)
{
   return ABS_DPL(a) <= EPS_DPL;
}

inline int ISINF(double a)
{
   return a >= INF_DPL - 1;
}

inline int ISEQ(double a, double b)
{
   return ISZERO(a - b);
}

inline double Lifting_GetTime()
{
	return (double)clock()/(double)CLOCKS_PER_SEC;
}

typedef double DTptype;
typedef double DTrctype;
typedef double DTwtype;
typedef double DTutype;
typedef double DTctype;
   
typedef struct DPLifting
{
   int                isleq;          // 1: <= knapsack; 0: >= knapsack
   DTptype*           p;              // lifting coefficients 
   DTwtype*           w;              // knapsack weights 
   DTutype*           u;              // variable upper bounds 
   int*               seed;           // indices in the seed inequality
   int                n_seed;         // nembers of indices in the seed inequality
   int*               liftingorder;   // variables lifted sequence
   int                n_liftingorder; 
   int*               isuseub;        // 1: fix variable at upper bound; 0: fix variable at lower bound
   DTctype            subcap;         // residual knapsack capacity
   DTctype            cap;            // table size for DPT / DPL
   DTctype            maxcap;
   DTctype            minweight;      // original >= threshold (isleq = 0 only) 
   int                n;
   int                solvedsize;
   double*            x_old;
   double*            x;
   double             activity;
   double*            dptlist;         // dense DPT table when isDPL = 0 

   DTptype*           psum;           // DPL breakpoint profits 
   DTwtype*           wsum;           // DPL breakpoint weights
   int                n_soltable;
   int                maxsolsize;
   double             rhs;
   double             tableleft;     // REDUCTION: U^k = sum u_i*w_i on remaining bounded vars
   int                reduction_active; // 0 if any variable is unbounded (disable REDUCTION)
   int                reduction_usable; // 1 if U^k>0 and subcap>0 after init (REDUCTION enabled)
   double             reduction_U_init; // U^k at init (after Calsubcap)
   double             reduction_b_init; // b^k = subcap at init

   DTptype*           psum1;
   DTwtype*           wsum1;
   DTptype*           psum2;
   DTwtype*           wsum2;

   int                force_mode;       // DPLIFTING_MODE_*; only MODE_THRESHOLD mid-lifts
   int                isDPL;             // 1: DPL table is authoritative; 0: dptlist is */

   double             threshold;       // τ for MODE_THRESHOLD and +R AUTO; <=0 → beta_th*mean(w)
   double             duration;        // measured CPU seconds (set by lifting())
   double             time_limit;      // optional limit in seconds; <=0 means unlimited
   double             t_start;         // Lifting_GetTime() when lifting() starts

   /* Empirical MODE_AUTO thresholds (0 → DPLIFTING_DEFAULT_*) */
   double             rho_th;
   double             beta_th;
   double             u_bar_th;
   double             lambda_th;       // +R AUTO: disable if λ < lambda_th (0 → default 0.1)
   int                reduction_request; // RED_AUTO / RED_ON / RED_OFF
   double             switch_cap;      // resolved τ
   double             feat_rho;
   double             feat_beta;
   double             feat_ubar;
   double             feat_lambda;     // λ = cap / Σ w_i u_i (0 if undefined)
} Lifting;

/** Row features for backend / reduction selection. */
typedef struct DPLiftingFeatures {
   double rho_w;   /**< w_max / w_min (1 if degenerate) */
   double beta;    /**< cap / mean(w) */
   double u_bar;   /**< mean(u); INF treated as large */
   double lambda;  /**< cap / Σ w_i u_i; 0 if unbounded / empty */
   double w_mean;
   double w_min;
   double w_max;
} DPLiftingFeatures;

/** Tunable empirical map for MODE_AUTO; zeros mean defaults. */
typedef struct DPLiftingPolicy {
   double rho_th;
   double beta_th;
   double u_bar_th;
   double lambda_th; /**< +R AUTO off when λ < this; 0 → DPLIFTING_DEFAULT_LAMBDA_TH */
   int prefer_dpl_fuzzy; /**< 1: boundary band favors DPL (default) */
} DPLiftingPolicy;

/** Same type as @c Lifting (struct tag @c DPLifting). */
typedef Lifting DPLifting;

void Lifting_Printsum(DPLifting* lift);
void Lifting_Check(DPLifting* lift);

int Lifting_Alloc(DPLifting* lift, int len, int scale, double threshold);
int Lifting_Realloc(DPLifting* lift, int len);
int Lifting_Reset(DPLifting* lift, int len);
int Lifting_Free(DPLifting* lift);
int Lifting_Calsubcap(DPLifting* lift);
int Lifting_Calcap(DPLifting* lift);

// Initialise lifting state from knapsack data and lifting sequence
int Lifting_Init(
      DPLifting* lift, 
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub, 
      DTctype cap, int issubcap, 
      int* basis, int n_basis, 
      int* liftingorder, int n_liftingorder, 
      int isleq, double* x, DTctype maxRhs, int n);

int Lifting_Wiszero(DPLifting* lift, DTptype p, DTwtype w, DTutype u);
int Lifting_Piszero(DPLifting* lift, DTptype p, DTwtype w, DTutype u);

// Merge a bounded item (p, w) into the DPL table
int Lifting_Mergesort(DPLifting* lift, DTptype p, DTwtype w);

// Merge an item with effectively unbounded multiplicity
int Lifting_Mergesortinf(DPLifting* lift, DTptype p, DTwtype w);

// Add item (p, w) via binary splitting; respects force_mode / threshold
int Lifting_Multiply(DPLifting* lift, DTptype p, DTwtype w, DTutype u);

// Query DPL table: last breakpoint with wsum <= cap or min profit with wsum >= cap
int Lifting_Findind(DPLifting* lift, DTctype cap, int begin, int end, int isleq);
DTptype Lifting_Findsol(DPLifting* lift, DTctype cap, int begin, int end, int isleq);

// Build seed table and return initial rhs of the cover inequality
DTptype Lifting_Calinitrhs(DPLifting* lift);

int Lifting_Iter(DPLifting* lift, DTptype p, DTwtype w, DTutype u);
int Lifting_Lifting(DPLifting* lift, DTptype* rhs);

/**
 * Run full coefficient lifting for a cover inequality
 *
 * @param cap      Knapsack rhs 
 * @param isSubCap If 1, cap is the subproblem capacity; else global rhs
 * @param seed     Variables fixed in the seed inequality
 * @param liftingorder  Remaining variables in lifting order
 * @param isLeq    1 for <= knapsack; 0 for >= knapsack
 * @param threshold  Capacity τ for MODE_THRESHOLD (initial+mid) and +R AUTO.
 *                   Ignored for MODE_AUTO / MODE_DPL / MODE_DPT mode choice.
 *                   If <= 0 under MODE_THRESHOLD or +R AUTO: τ = beta_th*mean(w).
 * @param duration   Optional time limit in seconds; <= 0 means unlimited.
 *                   Measured runtime is written to lift->duration.
 * @param isdpl_mode  MODE_AUTO | MODE_THRESHOLD | MODE_DPL | MODE_DPT
 * @return 1 on success; 0 on failure. Writes lifted coefficients into p and rhs.
 */
DPLIFTING_API int lifting(
      DPLifting* lift,
      DTptype* p, DTwtype* w, DTutype* u, int* isuseub,
      DTctype cap, int isSubCap,
      int* seed, int n_seed,
      int* liftingorder, int n_liftingorder,
      double* rhs,
      int isLeq, double* x, int n, double threshold, double duration, int isdpl_mode);

/** Fill policy with library defaults (rho=beta=6, u_bar=3, lambda_th=0.1, fuzzy DPL). */
DPLIFTING_API void dplifting_policy_default(DPLiftingPolicy* pol);

/** Compute rho_w, beta, u_bar from a knapsack row. */
DPLIFTING_API void dplifting_compute_features(
      const DTwtype* w, const DTutype* u, int n, DTctype cap, DPLiftingFeatures* feat);

/**
 * Select DPLIFTING_MODE_DPL or _DPT from features + policy (MODE_AUTO).
 * @return DPLIFTING_MODE_DPL or DPLIFTING_MODE_DPT
 */
DPLIFTING_API int dplifting_select_backend(
      const DPLiftingFeatures* feat, const DPLiftingPolicy* pol);

/** Effective residual capacity \(\bar b^k = \min(b^k, U^k)\) (U^k only if reduction active). */
DPLIFTING_API double Lifting_bar_b(const DPLifting* lift);

int lifting_lifting(DPLifting* lift, DTptype* alpha, DTwtype* a, DTutype* u, int* isuseub, DTptype *rhs, int n, int isleq);
void Lifting_Printsoltable(DTptype* psum, DTwtype* wsum, int n);

// Up-lifting
int Lifting_Up(DPLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs);

// Down lifting
int Lifting_Down(DPLifting* lift, DTptype* alpha, DTwtype a, DTutype u, DTptype *rhs);

int Lifting_Compress(DPLifting* lift, int begin);
int Lifting_Expand(DPLifting* lift);

//One unbounded knapsack DPT
void Lifting_DPTiterInf(DPLifting* lift, int w, double p);

// One bounded knapsack DPT
void Lifting_DPTiter(DPLifting* lift, int w, double p);
void Lifting_DPTPrint(double* dp, int c);
void Lifting_DPTFree(double* dp);

#ifdef __cplusplus
extern "C" {
#endif

/** Return codes for dplifting_lift_cover */
#define DPLIFTING_OK           0
#define DPLIFTING_ERR_ALLOC   -1
#define DPLIFTING_ERR_ARGS    -2
#define DPLIFTING_ERR_INTERNAL -3

/**
 * Lift a cover inequality for a single knapsack row (C ABI).
 *
 * Modifies @p coef in place and writes the lifted inequality rhs to @p rhs.
 *
 * @param isdpl_mode  DPLIFTING_MODE_AUTO / _THRESHOLD / _DPT / _DPL
 */
DPLIFTING_API int dplifting_lift_cover(
      int n,
      double* coef,
      const double* weight,
      const double* ub,
      const int* use_ub,
      double cap,
      int is_subcap,
      const int* seed,
      int n_seed,
      const int* lifting_order,
      int n_order,
      double* rhs,
      int is_leq,
      double threshold,
      int isdpl_mode,
      const double* x_frac);

#ifdef __cplusplus
}
#endif

/* Integer benches still pass -DDLLIFTING_REDUCTION and use pre-rename names. */
#if defined(DLLIFTING_REDUCTION) && !defined(DPLIFTING_REDUCTION)
#  define DPLIFTING_REDUCTION
#endif
#ifndef DPLIFTING_NO_LEGACY_ALIASES
typedef DPLifting DLLifting;
#define EPS_DL EPS_DPL
#define INF_DL INF_DPL
#define DLLIFTING_MODE_AUTO      DPLIFTING_MODE_AUTO
#define DLLIFTING_MODE_THRESHOLD DPLIFTING_MODE_THRESHOLD
#define DLLIFTING_MODE_DP        DPLIFTING_MODE_DPT
#define DLLIFTING_MODE_DL        DPLIFTING_MODE_DPL
#endif

#endif
