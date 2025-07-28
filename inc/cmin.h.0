#ifndef CMIN
#define CMIN

#include <gsl/gsl_errno.h> 
#include <gsl/gsl_math.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_multimin.h>
#include "cdef.h"
#include "ndef.h"
#include "edef.h"
#include "cell.h"
#include "cmod.h"
#include "cutl.h"

void   STR_BUF(gsl_vector *x);
void   BUF_STR(const gsl_vector *x);
double cfunc_gsl(const gsl_vector *x, void *params);
void   cdfunc_gsl(const gsl_vector *x, void* params, gsl_vector *d);
void   cfdfunc_gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df);
double CELL_MIN(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
#endif
