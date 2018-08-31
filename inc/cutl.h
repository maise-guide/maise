#ifndef CUTL
#define CUTL

void   ROT_CELL(Cell *C, int k, int q0, int q1);
void   CENTER(Cell *C, double o);
void   ROTATE_CELL(Cell *C, double f);
void   NANO_ROT(Cell *C, int J);
double Cell_VOLUME(Cell *C);
void   Lat_Order(Cell *C) ;
void   Lat_Align(Cell *C);
void   Copy_C(Cell *C, Cell *C1);
void   Clone(Cell *C, Cell *C1, int N0, int N1, int N2) ;
void   abc(Cell *C);
void   ABC_LT(Cell *C);
void   KMESHOLD(Cell *C, int N, char name[]);
void   KMESH(Cell *C, double KM, char name[], int ND);
int    Check_OUTCAR(Cell *C, char file[200]);
int    Check_GULP_OUT(Cell *C, char file[200]);
void   atomSwap(Cell *C, int i, int j, int k);
#endif
