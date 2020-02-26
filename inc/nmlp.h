#ifndef NMLP
#define NMLP

#include <omp.h>
#include "clib.h"
#include "cdef.h"
#include "ndef.h"
#include "nmin.h"
#include "nutl.h"
#include "util.h"

double GPc(int G, double x);
double GPp(int G, double x);
double GP(int G, double x);
double FRC_ANN_PARA(ANN *R, LNK *L, double ***e, double ***d);
double DF_ANN_PARA(ANN *R, LNK *L, double ***Bp, double ****Wp,  double ***e,  double ***d,  double ***c, double ***xn, double ***xm);
double DE_ANN_PARA(ANN *R, LNK *L, double ***Bp, double ****Wp,  double ***e);
double TOT_ERR(ANN *R, LNK *L);
double ENE_ANN(ANN *R, LNK *L);
void   frc_ann(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
double FRC_ANN(ANN *R, LNK *L);
void   INIT_MLP(ANN *R);
int    W4V(ANN *R, double *V, int *SPC);
void   W2V(ANN *R, double *V);
void   V2W(ANN *R, double *V);
void   TRAN_MLP(ANN *R, Cell *C, LNK *L);

#endif
