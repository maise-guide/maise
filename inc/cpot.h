#ifndef CPOT
#define CPOT

#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "util.h"

double dFC(double R1, double R2, double r);
double FC(double R1, double R2, double r);
double ENE_POT(Cell *C);
double FRC_POT(Cell *C);
void   READ_POT(Cell *C, char *file);
#endif
