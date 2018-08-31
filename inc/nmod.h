#ifndef NMOD
#define NMOD

void   NNET_MAIN(ANN *R, PRS *P, Cell *C);
void   TRAN_ANN(ANN *R, Cell *C);
void   CHCK_ERR(ANN *R, LNK *L);
double CPU_TIME(double ti, char buf[200]);
void   LOAD_LNK(ANN *R, Cell *C, LNK *L);
void   OUT_ANN(ANN *R, double time, char *s1,char *s2, char *s3);
void   PLT_EVAL(char *path,char *s);
double EVAL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L);
void   EVAL_ANN(ANN *R, PRS *P, Cell *C, LNK *L);
#endif
