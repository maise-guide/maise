#ifndef UTIL
#define UTIL

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "clib.h"
#include "cdef.h"
#include "cell.h"
#include "cutl.h"
#include "edef.h"
#include "ndef.h"
#include "nutl.h"
#include "util.h"

struct Node
{
  char line[200];
  struct Node *next;
};

extern const double Pi;
extern const int    NE;
extern const double sqrt3;
extern const double sqrt2;
extern const double au;
extern const double pi;
extern const double eV2kV;
extern const double eV2GPa;

void   PRNT_HEAD();
void   FPRNT_HEAD(FILE *out);
int    READ_CELL(Cell *C, char filename[]);
void   SAVE_CELL(Cell *C, char filename[], int SD) ;
void   Read_OUTCAR(Cell *C, char file[], int NC);
double Read_OSZI(char *file);
void   PRNT_LIST(Cell *C);
void   Print_LOG(char buf[]);
void   READ_MAIN(Tribe *T, ANN *R, PRS *P, Cell *C, int J, int ARGC);
void   EV (double A[3][3], double e[3][3], double b[3]);
void   EVN(double *A, double *B, double *e, double *b, int N);
void   EXIT(char *s);
double RANG();
void   VectorProd(double *a, double *b, double *c);
double VectorLen(double *c, int D);
void   VectorNorm(double *c) ;
double CrossProd(double *a, double *b, double *c);
double COS(double *x1, double *x2, int D);
double DiffLen(double *x1, double *x2, int D);
double DotProd(double *x1, double *x2, int D);
double sign(double a);
long   TIME(char file[200]);
int    TIME_DIF(int t1, int t2);
void   Sort(double *x, int *I, int N);
int    FIXED(int *vec);
void   inSwap(double D[1000], int I[1000][3], int i, int n);
void   iSwap(int *i, int *j);
void   dSwap(double *x, double *y);
void   Counter(int *ic, int *N, int i);
double dMin(double a, double b);
void   tprintf(struct Node** tmp_ref, char new_line[200]);
char   *tgets(char *buf,int N,struct Node **n);
void   tclose(struct Node **n);
double *make_d1D(int x);
void   free_d1D(double *data);
double **make_d2D(int x, int y);
void   free_d2D(double **data, int x);
double ***make_d3D(int x, int y, int z);
void   free_d3D(double ***data, int x, int y);
double ****make_d4D(int x, int y, int z,int w);
void   free_d4D(double ****data, int x, int y,int z);
int    *make_i1D(int x);
void   free_i1D(int *data);
int    **make_i2D(int x, int y);
void   free_i2D(int **data, int x);
int    ***make_i3D(int x, int y, int z);
void   free_i3D(int ***data, int x, int y);
int    ****make_i4D(int x, int y, int z,int w);
void   free_i4D(int ****data, int x, int y,int z);
double Random(void);
void   PlantSeeds(long x);
void   PutSeed(long x);
void   GetSeed(long *x);
void   SelectStream(int index);
void   TestRandom(void);
int    check_ver(char *fname);
void   READ_MODEL(ANN *R, PRS *P, Cell *C);

#endif
