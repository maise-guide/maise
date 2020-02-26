#ifndef NMIN
#define NMIN

#include <gsl/gsl_errno.h> 
#include <gsl/gsl_math.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_multimin.h>
#include <signal.h>
#include <omp.h>
#include "clib.h"
#include "cdef.h"
#include "ndef.h"
#include "nmlp.h"
#include "nutl.h"
#include "util.h"

void   sig_term_handler(int signum, siginfo_t *info, void *ptr);
void   catch_sigterm();
void   W2X_GSL(gsl_vector *x);
void   X2W_GSL(const gsl_vector *x);
double tot_err_gsl();
double func_gsl(const gsl_vector *x, void *params) ;
void   dfunc_num_gsl(gsl_vector *x, void *params, gsl_vector *d) ;
void   dfunc_gsl(const gsl_vector *x, void* params, gsl_vector *d) ;
void   fdfunc_gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df) ;
double MLP_MIN(ANN *R, LNK *L);
#endif
