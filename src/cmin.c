#include "cdef.h"
#include "ndef.h"
#include "edef.h"
#include "cell.h"
#include "cmod.h"
#include "cutl.h"
#include "cmin.h"

#define x(X)  gsl_vector_get(x, X)
#define d(X)  gsl_vector_get(d, X)
#define p(X)  gsl_vector_get(p, X)
/////////// get p?
//#define o(X)  gsl_vector_get(o, X)
#define set(V, X, Y) gsl_vector_set(V, X, Y)
long ECNT,FCNT;

double X32;

int   CLL;

ANN  *RRR;
PRS  *PPP;
PRS  *WWW;
Cell *CCC;
LNK  *LLL;

//======================================================================
void STR_BUF(gsl_vector *x)
{
  int i,q,j=0;
  
  Relative(CCC);
  if(CCC->RLXT!=7)
  for(i=0;i<CCC->N;i++)
    for(q=0;q<3;q++)
      if(CCC->FF[i][q]==1)
	set(x,j++,CCC->X[i][q]);
  if(CCC->RLXT==2)
    return;
  set(x,j++,CCC->L[0][0]);
  set(x,j++,CCC->L[1][1]);
  set(x,j++,CCC->L[2][2]);
  set(x,j++,CCC->L[0][1]);
  set(x,j++,CCC->L[1][2]);
  set(x,j++,CCC->L[2][0]);
  set(x,j++,CCC->L[1][0]);
  set(x,j++,CCC->L[2][1]);
  set(x,j++,CCC->L[0][2]);

  return;
}
//======================================================================
void BUF_STR(const gsl_vector *x)
{
  int i,q,j=0;

  if(CCC->RLXT!=7)
  for(i=0;i<CCC->N;i++)
    for(q=0;q<3;q++)
      if(CCC->FF[i][q]==1)
        CCC->X[i][q] = x(j++);
  if(CCC->RLXT==2)
  {
    Real(CCC);
    return;
  }
  CCC->L[0][0] = x(j++);
  CCC->L[1][1] = x(j++);
  CCC->L[2][2] = x(j++);
  CCC->L[0][1] = x(j++);
  CCC->L[1][2] = x(j++);
  CCC->L[2][0] = x(j++);  
  CCC->L[1][0] = x(j++);
  CCC->L[2][1] = x(j++);
  CCC->L[0][2] = x(j++);

  Real(CCC);
}
//======================================================================
double cfunc_gsl(const gsl_vector *x, void *params) 
{  
  double E;
  int EFS;

  EFS = PPP->EFS;
  PPP->EFS = 0;

  BUF_STR(x);
  E = CELL_ENE(RRR,PPP,WWW,CCC,LLL);  
  PPP->EFS = EFS;
  STR_BUF((gsl_vector*) x);

  ECNT++;
  return E;
}
//======================================================================
void cdfunc_gsl(const gsl_vector *x, void* params, gsl_vector *d) 
{
  int i,q,k,j=0;
  double V,s[3];

  FCNT++;
  BUF_STR(x);
  V = Cell_VOLUME(CCC);
  Reciprocal(CCC);
  CELL_FRC(RRR,PPP,WWW,CCC,LLL);
  
  if(CCC->RLXT!=7)
    for(i=0;i<CCC->N;i++)//,printf("\n"))
    for(q=0;q<3;q++)
      if(CCC->FF[i][q]==1)
	for(k=0,set(d,j++,0.0);k<3;k++)
	  set(d,j-1,d(j-1) - CCC->F[i][k]*CCC->L[q][k] );

  if(CCC->RLXT==2)
  {
    STR_BUF((gsl_vector*) x);
    return;
  }
  //  0  1  2  3  4  5
  // xx yy zz xy yz zx

  for(q=0;q<3;q++)
  {
    s[q] = CCC->U[q];
    CCC->U[q] -= CCC->p*V;
  }

  set(d,j++,-( CCC->U[0]*CCC->R[0][0] + CCC->U[3]*CCC->R[0][1] + CCC->U[5]*CCC->R[0][2] ) );
  set(d,j++,-( CCC->U[3]*CCC->R[1][0] + CCC->U[1]*CCC->R[1][1] + CCC->U[4]*CCC->R[1][2] ) );
  set(d,j++,-( CCC->U[5]*CCC->R[2][0] + CCC->U[4]*CCC->R[2][1] + CCC->U[2]*CCC->R[2][2] ) );

  set(d,j++,-( CCC->U[3]*CCC->R[0][0] + CCC->U[1]*CCC->R[0][1] + CCC->U[4]*CCC->R[0][2] ) );
  set(d,j++,-( CCC->U[5]*CCC->R[1][0] + CCC->U[4]*CCC->R[1][1] + CCC->U[2]*CCC->R[1][2] ) );
  set(d,j++,-( CCC->U[0]*CCC->R[2][0] + CCC->U[3]*CCC->R[2][1] + CCC->U[5]*CCC->R[2][2] ) );

  set(d,j++,-( CCC->U[0]*CCC->R[1][0] + CCC->U[3]*CCC->R[1][1] + CCC->U[5]*CCC->R[1][2] ) );
  set(d,j++,-( CCC->U[3]*CCC->R[2][0] + CCC->U[1]*CCC->R[2][1] + CCC->U[4]*CCC->R[2][2] ) );
  set(d,j++,-( CCC->U[5]*CCC->R[0][0] + CCC->U[4]*CCC->R[0][1] + CCC->U[2]*CCC->R[0][2] ) );

  for(q=0;q<3;q++)
    CCC->U[q] = s[q];

  for(i=0;i<j;i++)
    if(fabs(d(i))<1e-10)
      set(d,i,0.0);

  STR_BUF((gsl_vector *) x);

  return;
}
//======================================================================
// Compute both f and df together.
//======================================================================
void cfdfunc_gsl(const gsl_vector *x, void *params, double *f, gsl_vector *df) 
{
  *f = cfunc_gsl(x, params); 
  cdfunc_gsl(x, params, df);
}
//======================================================================
//======================================================================
double CELL_MIN(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  gsl_vector *p;
  void *params;
  double E,O;
  int i,q,iter,N,status;
  char min[200];

  RRR = R; PPP = P; CCC = C; LLL = L; WWW = W;

  CLL=1;
  ECNT = FCNT = 0;
  PPP->EFS  = 3;
  params = NULL;

  N = 0;
  if(CCC->RLXT!=7)
  for(i=0;i<CCC->N;i++)
    for(q=0;q<3;q++)
      if(CCC->FF[i][q]==1)
	N++;
  if(CCC->RLXT==3||CCC->RLXT==7)
    N += 9;

  if(RRR->MINT==0)
    sprintf(min,"%s","BFGS2");
  if(RRR->MINT==1)
    sprintf(min,"%s","CG-PR");
  if(RRR->MINT==2)
    sprintf(min,"%s","CG-FR");
  if(RRR->MINT==3)
    sprintf(min,"%s","STPDT");

  printf("\n  %s relaxation: %d DOF\n\n",min,N);

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_multimin_function_fdf gsl_fdf;

  gsl_fdf.n = N;
  gsl_fdf.f = cfunc_gsl;
  gsl_fdf.df = cdfunc_gsl;
  gsl_fdf.fdf = cfdfunc_gsl;
  gsl_fdf.params = params;

  p = gsl_vector_alloc (N);
  STR_BUF(p);
  O = cfunc_gsl(p, params);
  iter = 0;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  if(RRR->MINT==0)
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
  if(RRR->MINT==1)
    T = gsl_multimin_fdfminimizer_conjugate_pr;
  if(RRR->MINT==2)
    T = gsl_multimin_fdfminimizer_conjugate_fr;
  if(RRR->MINT==3)
    T = gsl_multimin_fdfminimizer_steepest_descent;

  s = gsl_multimin_fdfminimizer_alloc (T, N);
  printf("%5d %24.16lf\n",iter, C->H);

  E = O;

  gsl_multimin_fdfminimizer_set (s, &gsl_fdf, p, 0.01, 0.001);
  do
  {
    C->it++;
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    status = gsl_multimin_test_size( fabs(C->H-E), RRR->ETOL);

    printf("%5d % 24.16lf % 24.16lf\n",iter, C->H,C->H-E);
    E = C->H;

    if( C->OUT%10==2)
    {
      Real(C);
      CELL_OUT(C);    
      Relative(C);
    }
  }
  while (status == GSL_CONTINUE && iter < RRR->MITR );

  E = C->H;
  BUF_STR(s->x);

  printf("\n  ENERGY CALLS = %5ld\n"  ,ECNT);
  printf(  "  FORCE  CALLS = %5ld\n\n",FCNT);

  printf("  INI % 16.8lf \n",O);
  printf("  DIF % 16.8lf \n",E-O);
  printf("  FIN % 16.8lf \n",E);

  gsl_vector_free (p);
  gsl_multimin_fdfminimizer_free (s);

  return 0.0;
}
