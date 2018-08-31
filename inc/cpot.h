#ifndef CPOT
#define CPOT

double dFc(double Rc, double r);
double Fc(double Rc, double r);
double dFC(double R1, double R2, double r);
double FC(double R1, double R2, double r);
double ENE_GP3D(Cell *C, int J);
double ENE_GP(Cell *C, int J);
double ENE_SC(Cell *C, int J);
void   frc_sc(Cell *C);
double FRC_GP3D(Cell *C, int J);
double FRC_GP(Cell *C, int J);
double FRC_SC(Cell *C, int J);
double ENE_GP_(ANN *R, LNK *L);
double ENE_SC_(ANN *R, LNK *L);
double DE_SC_(ANN *R, LNK *L);
double ENE_LJ(Cell *C);
void   frc_lj(Cell *C);
double FRC_LJ(Cell *C);
void   SAVE_SC(ANN *R, double time);
void   READ_MOD(Cell *C, char *file);
void   INIT_MOD(ANN *R, Cell *C);
double scfunc(ANN *R, LNK *L);
void   scfunc_num(ANN *R, LNK *L);
void   scdfunc(ANN *R, LNK *L);
void   TRAN_SC(ANN *R, Cell *C, LNK *L);
#endif
