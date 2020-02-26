#ifndef EMOD
#define EMOD

#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "cfnc.h"
#include "cmod.h"
#include "cpot.h"
#include "cutl.h"
#include "edef.h"
#include "efnc.h"
#include "eutl.h"
#include "ndef.h"
#include "nmod.h"
#include "nutl.h"
#include "sutl.h"
#include "util.h"

void   INIT_TR(Tribe *T);
void   EXIT_TR(Tribe *T);
void   RANK_TR(Tribe *T);
void   SLCT_TR(Tribe *T);
void   PLOT_TR(Tribe *T);
void   QSUB_TR(Tribe *T, int p);
void   RELX_INT(Tribe *T);
void   RELX_TR(Tribe *T);
void   EVLV_TR(Tribe *T, int J);
void   ANA_EVOS(Tribe *T, Cell *C, Cell *D);
void   INIT_EVOS(Tribe *T, Cell *C);
void   EVOS_MAIN(Tribe *T, ANN *R, PRS *P, Cell *C);
#endif
