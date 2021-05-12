#ifndef CFNC
#define CFNC

#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "cutl.h"
#include "nutl.h"
#include "sutl.h"

extern const double Pi;
extern const int D3;

int    FIND_NDIM(Cell *C);
void   KILL_DBL (Cell *C, double tol);
void   APPL_SG  (Cell *C, double tol);
void   READ_CIF (Cell *C, char file[], double tol, int NM, char input[]);
int    FIND_MTY (Cell *C, double tol);
void   FIND_PRS (Cell *C, Cell *D, double tol);
void   FIND_CXC (Cell *C, Cell *D, int argc, char argv[20][200]);
void   COMP_STR (Cell *C, Cell *D, int argc, char argv[20][200]);
void   INIT_CELL(Cell *C, char filename[], int N, int NM, int J);
void   CELL_EXAM(Cell *C, Cell *D, int argc, char **argv);

#endif
