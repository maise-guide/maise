#ifndef SUTL_H
#define SUTL_H

#include <gsl/gsl_math.h>
#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "cutl.h"
#include "ndef.h"
#include "nutl.h"
#include "spglib.h"

int   FIND_WYC(Cell *C, Cell *D,  double tol, int J);
void  READ_SG(Cell *C);
#endif
