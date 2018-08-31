#ifndef EFNC
#define EFNC  
//includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cdef.h"
#include "edef.h"
#include "ndef.h"
#include "nutl.h"
#include "cutl.h"
#include "cell.h"
#include "util.h"
#include "eutl.h"

extern const double Pi;

void   GEN_TYPE(Tribe *T, Cell *C, double *D, int *I);
void   NANO_ORDER_TYPE(Tribe *T, int J);
void   NANO_CHOP(Tribe *T, int J, Cell *C, int P);
void   NANO_BLOB(Tribe *T, int J, Cell *C, int P);
void   NANO_TETR(Tribe *T, int J, Cell *C, int P);
void   NANO_PLNT(Tribe *T, int J, int P);
void   NANO_PACK(Tribe *T, int J, Cell *C, int P);
void   NANO_SWAP(Tribe *T, int J);
void   NANO_MATE(Tribe *T, int J);
void   NANO_RUBE(Tribe *T, int J);
void   NANO_SYMM(Tribe *T, int J);
void   NANO_MUTE(Tribe *T, int J);
void   BULK_MATE(Tribe *T, int J);
void   BULK_MUTE(Tribe *T, int J);
#endif
