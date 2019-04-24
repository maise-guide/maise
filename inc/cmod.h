#ifndef CMOD
#define CMOD

double KINETIC(Cell *C);
double P_Lindemann(Cell *C, int J);
double Lindemann(Cell *C, int J);
double CELL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
double CELL_FRC(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   MDOUT(Cell *C, int PIPE, int DISK, int NS, int n, int k, double E0, double T,double li);
void   Maxwell(Cell *C, double T, int therm);
void   Dynamics(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L, double T, int PIPE);
void   Thermostat(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L, double T, int PIPE) ;
void   CELL_MD(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   CELL_OUT(Cell *C);
int    CELL_OK(Cell *C, double *Rm);
int    STOP_OK(Cell *C, double *Rm, double x);
void   CELL_RELX(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   CELL_PHON(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   CELL_MAIN(ANN *R, PRS *P, Cell *C);

#endif
