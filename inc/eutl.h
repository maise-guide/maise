#ifndef EUTL
#define EUTL
 
#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "cutl.h"
#include "util.h"

int    CHCK_Rm(Cell *C, double *Rm, double s);
int    ADJT_NP(Cell *C, double *Rm, double RM, double *b, int N);
int    ADJT_CL(Cell *C, double *Rm, int N);
void   MATE_LT(Cell *C1, Cell *C2, Cell *C, double s);
int    SHKE_CL(Cell *C, double dL, double dX);
void   RAND_VC(double *a);
void   RAND_LV(Cell *C);
void   RAND_CL(Tribe *T, Cell *C, Cell *D, int J);
double COMP_CL(Cell *C, Cell *D);
void   TEMP_CL(Tribe *T, Cell *C, int p);
#endif
