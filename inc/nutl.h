#ifndef NUTL
#define NUTL

void   READ_STRAT3(ANN *R);
void   READ_STRAT2(ANN *R);
void   ANA_STR(ANN *R);
void   RNDNM2(int *NM, int N,int sd);
double cpu_time( );
void   MARK_CL(Cell *C, int M, long seed);
void   Build_ANN(ANN *R);
void   Build_LNK(LNK *L, int N, int NM, int D, int EFS);
void   SAVE_ANN(ANN *R, double TIME);
void   READ_ANN(ANN *R);
int    symb_atom(char *s    );
void   atom_symb(int i,char *s );
double NP_VOL(Cell *C);

#endif
