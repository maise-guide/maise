#ifndef NPRS
#define NPRS

void   GNcheck(char *s, int coor);
void   Build_PRS(PRS *P, PRS *W, int J);
double dfc(double Rc, double r);
double fc(double Rc, double r);
double dGNR(PRS *P, int n, double Rc, double r, double fcij, double dfcij);
double dgnr(PRS *P, int n, double Rc, double r);
double GNR(PRS *P, int n, double Rc, double r);
double GNA(PRS *P, int n, double R0, double n0, double r0, double R1, double n1, double r1, double R2, double n2, double r2, double a);
void   LNK_IN(LNK *L, int o, char *path);
void   LNK_OUT(LNK *L, int o, char *path, int EFS, int D);
void   PRS_STOP(Cell *C);
void   PRS_BP(PRS *P, PRS *W, Cell *C, LNK *L, int o, char *path);
double pl(PRS *P, int l, double x);
double gn(PRS *P, int n, double r);
void   PRS_PS(PRS *P, Cell *C, LNK *L, int o, char *path);
void   PRS_IJ(PRS *P, Cell *C, LNK *L, int o, char *path);
void   PARS_STR(PRS *P, PRS *W, Cell *C, LNK *L, int o, char *path);
int    CHCK_DAT(ANN *R, Cell *C, LNK *L);
void   SORT_FIT(int *NFIT, double *EFIT, int n, int N, int M, double EMAX, int TAG);
void   PARS_DAT(ANN *R, PRS *P, Cell *C, LNK *L);

#endif
