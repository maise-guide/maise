#ifndef CELL
#define CELL

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "clib.h"
#include "cdef.h"
#include "cutl.h"
#include "util.h"

extern const double Pi;

void   Build_Cell(Cell *C, int J);
void   LIST(Cell *C, int O);
double NDR(Cell *C, int i, int j);
double DX(Cell *C, int i, int n, int q);
double ndx(Cell *C, int i, int j, int k);
double NDX(Cell *C, int i, int n);
double NDD(Cell *C, int i, int n);
double Cos(Cell *C, int i, int j, int k);
int    NN(Cell *C, int i, int j);
void   SCALE_Cell(Cell *C, double a);
void   SCALE_LATTICE(Cell *C, double a);
void   MATCH_LATTICE(Cell *C, Cell *D);
void   Print_RDF(Cell *C, char *file);
void   Print_RDF_FULL(Cell *C, char *file);
void   RDF(Cell *C, int J);
double CxC(Cell *C, Cell *D);
void   Reciprocal(Cell *C);
void   Relative(Cell *C);
void   Real(Cell *C);
void   JAR(Cell *C);
int    SHRT_LV(Cell *C);
void   ORDER_Z(Cell *C);
void   ORDER(Cell *C);
void   ADD(Cell *C, Cell *D, double x, double y, double z);
#endif
