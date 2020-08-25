#ifndef CMOD
#define CMOD

#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "cpot.h"
#include "cmin.h"
#include "cutl.h"
#include "ndef.h"
#include "nprs.h"
#include "nmlp.h"
#include "nutl.h"
#include "util.h"

double CELL_PRS(Cell *C);
double CELL_KIN(Cell *C);
double Lindemann(Cell *C, int J);
double CELL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
double CELL_FRC(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   MDOUT(Cell *C, double PE,double KE,double LV,double IP,double DE,double LI,int PIPE, int DISK, int NS, int sam, double T,double E0,double li);
void   Maxwell(Cell *C, double T, int therm);
void   Dynamics(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L, double T, int PIPE, char *fname);
void   CELL_MD(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   CELL_OUT(Cell *C);
int    CELL_OK(Cell *C);
void   CELL_RELX(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   CELL_PHON(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   CELL_MAIN(ANN *R, PRS *P, Cell *C);
double CALL_MAISE(ANN *R, PRS *P, PRS *W, LNK *L, Cell *C, int CODE, int N, int NM, int ND, int NP, int XT, int *ATMN, double *LAT, double *X, double *FRC, double *STR);
#endif
