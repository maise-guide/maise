#ifndef NMOD
#define NMOD

#include "clib.h"
#include "cdef.h"
#include "ndef.h"
#include "cutl.h"
#include "cpot.h"
#include "nprs.h"
#include "nmlp.h"
#include "plot.h"

void   NNET_MAIN(ANN *R, PRS *P, Cell *C);
void   TRAN_ANN(ANN *R, Cell *C);
void   CHCK_ERR(ANN *R, LNK *L);
double CPU_TIME(double ti, char buf[200]);
void   LOAD_LNK(ANN *R, Cell *C, LNK *L);
void   ADJT_LNK(ANN *R, LNK *L);
void   OUT_ANN(ANN *R, LNK *L,double w_time, double c_time, double s_time, char *s1,char *s2, char *s3);
void   PLT_EVAL(char *path,char *s);
double EVAL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   EVAL_ANN(ANN *R, PRS *P, Cell *C, LNK *L);
#endif
