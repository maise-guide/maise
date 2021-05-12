#include "cutl.h"

//======================================================================
// Rotate unit cell along q0 axis to put atom k in the q0-q1 plane
//======================================================================
void ROT_CELL(Cell *C, int k, int q0, int q1)
{
  int i,q2;
  double r,ang,ang0,w;

  q2 = 3 - q1 - q0;
  w  = 1e-12;

  r    = sqrt(pow(C->X[k][q1],2.0)+pow(C->X[k][q2],2.0));
  if(r>w)
    ang0  = acos(C->X[k][q1]/r);
  else
    ang0  = 0.0;
  if(C->X[k][q2]<0.0)
    ang0 *= -1.0;

  for(i=0;i<C->N;i++)
  {
    r   = sqrt(pow(C->X[i][q1],2.0)+pow(C->X[i][q2],2.0));
    if(r>w)
      ang = acos(C->X[i][q1]/r);
    else
      ang  = 0.0;
    if(C->X[i][q2]<0.0)
      ang *= -1.0;
    
    C->X[i][q1] = r*cos(ang-ang0);
    C->X[i][q2] = r*sin(ang-ang0);
  }
  
}
//===================================================================
void CENTER(Cell *C, double o)
{
  double r;
  int i,q;

  Relative(C);
  for(q=0;q<3;q++)
  {
    for(i=0,r=0.0;i<C->N;i++)
      r += C->X[i][q]/(double)C->N;
    for(i=0;i<C->N;i++)
      C->X[i][q] += -r + o;
  }
  Real(C);
}
//======================================================================================================
void ROTATE_CELL(Cell *C, double f)
{
  double t;
  int i;
  for(i=0;i<C->N;i++)
  {
    t = C->X[i][0];
    C->X[i][0] =  t*cos(f) - C->X[i][1]*sin(f);
    C->X[i][1] =  t*sin(f) + C->X[i][1]*cos(f);
  }
  for(i=0;i<3;i++)
  {
    t = C->L[i][0];
    C->L[i][0] =  t*cos(f) - C->L[i][1]*sin(f);
    C->L[i][1] =  t*sin(f) + C->L[i][1]*cos(f);
  }
}
//===================================================================
// Will align NP along moments of inertia, needs C->N+2 array
//===================================================================
void NANO_ROT(Cell *C, int J)
{
  double II[3][3],e[3][3],b[3],r,w,k,a;
  int    i,j,q;  

  w = 1e-6;

  for(j=0;j<D3;j++)
    for(q=0;q<D3;q++)
      for(i=0,II[j][q]=0.0;i<C->N;i++)
        if(j==q)
          II[j][q] += (double)(C->ATMN[i]+1) * (C->X[i][(q+1)%3]*C->X[i][(q+1)%3] + C->X[i][(q+2)%3]*C->X[i][(q+2)%3]);
        else
          II[j][q] -= (double)(C->ATMN[i]+1) * C->X[i][q]*C->X[i][j];

  if(0)
  for(j=0;j<D3;j++,printf("\n"))
    for(q=0;q<D3;q++)
      printf("% 12.6lf ",II[j][q]);

  EV(II,e,b);

  if(0)
  for(j=0;j<D3;j++,printf("\n"))
    for(q=0,printf("% 12.6lf  ",b[j]);q<D3;q++)
      printf("% 12.6lf ",e[j][q]);

  Relative(C);
  C->N++;

  //===== if all eigenvalues are degenerate place furthest atoms along x and z =====
  if( fabs(b[0]-b[1])/b[0] < w && fabs(b[1]-b[2])/b[1] < w )
  {
    r = C->X[0][0]*C->X[0][0] + C->X[0][1]*C->X[0][1] + C->X[0][2]*C->X[0][2];
    for(i=1,k=0;i<C->N-1;i++)
      if( r < (a = C->X[i][0]*C->X[i][0] + C->X[i][1]*C->X[i][1] + C->X[i][2]*C->X[i][2]) )
      {
        r = a;
        k = i;
      }
    ROT_CELL(C,k,0,1);
    ROT_CELL(C,k,2,0);
    r = C->X[0][1]*C->X[0][1] + C->X[0][2]*C->X[0][2];
    for(i=1,k=0;i<C->N-1;i++)
      if( r < (a = C->X[i][1]*C->X[i][1] + C->X[i][2]*C->X[i][2]) )
      {
        r = a;
        k = i;
      }
    ROT_CELL(C,k,0,2);
  }  
  //===== if the last two eigenvalues are degenerate place furthest y-z atom along z =====
  if( fabs(b[0]-b[1])/b[0] >= w && fabs(b[1]-b[2])/b[1] < w )
  {
    for(q=0;q<3;q++)
      C->X[C->N-1][q] = e[0][q];
    ROT_CELL(C,C->N-1,0,1);
    ROT_CELL(C,C->N-1,2,0);
    r = C->X[0][1]*C->X[0][1] + C->X[0][2]*C->X[0][2];
    for(i=1,k=0;i<C->N-1;i++)
      if( r < (a = C->X[i][1]*C->X[i][1] + C->X[i][2]*C->X[i][2]) )
      {
        r = a;
        k = i;
      }
    ROT_CELL(C,k,0,2);
  }
  //===== if the first two eigenvalues are degenerate place furthest x-y atom along x =====
  if( fabs(b[0]-b[1])/b[0] < w && fabs(b[1]-b[2])/b[1] >= w )
  {
    for(q=0;q<3;q++)
      C->X[C->N-1][q] = e[2][q];
    ROT_CELL(C,C->N-1,2,1);
    ROT_CELL(C,C->N-1,0,2);
    r = C->X[0][0]*C->X[0][0] + C->X[0][1]*C->X[0][1];
    for(i=1,k=0;i<C->N-1;i++)
      if( r < (a = C->X[i][0]*C->X[i][0] + C->X[i][1]*C->X[i][1]) )
      {
        r = a;
        k = i;
      }
    ROT_CELL(C,k,2,0);
  }
  //===== if all the eigenvalues are different rotate along the eigenvectors  =====
  if( fabs(b[0]-b[1])/b[0] >= w && fabs(b[1]-b[2])/b[1] >= w )
  {
    C->N++;
    for(j=0;j<2;j++)
      for(q=0;q<3;q++)
        C->X[C->N-2+j][q] = e[j][q];
    ROT_CELL(C,C->N-2,0,1);
    ROT_CELL(C,C->N-2,2,0);
    ROT_CELL(C,C->N-1,0,1);
    C->N--;
  }
  C->N--;
  Real(C);      
  if(J==1)
    SAVE_CELL(C,"CONTCAR",0);

}
//==================================================================
//
//==================================================================
double CELL_VOL(Cell *C)
{   
  if(C->ND==0)
    return NP_VOL(C);

  double V = C->L[0][0]*C->L[1][1]*C->L[2][2] + C->L[0][2]*C->L[1][0]*C->L[2][1] + C->L[0][1]*C->L[1][2]*C->L[2][0]
           - C->L[0][0]*C->L[1][2]*C->L[2][1] - C->L[0][1]*C->L[1][0]*C->L[2][2] - C->L[0][2]*C->L[1][1]*C->L[2][0];
  if(0&&V<0.0)
    printf("%s cell has negative volume %lf!\n",C->TAG,V);
  return fabs(V);
} 
//==================================================================
void Lat_Order(Cell *C) // chooses 1st and 2nd lattice vectors be closest in magnitude
{
  int q,i;
  double a[3],b[3];

  for(i=0;i<3;i++)
    a[i] = VectorLen(C->L[i],3);
  for(i=0;i<3;i++)
    b[i] = pow( a[(i+1)%3]-a[i],2.0 ) + pow( a[(i+2)%3]-a[i],2.0 );

  if( b[2]<b[0] || b[2]<b[1] )
  {
    i = 0;
    if(b[1]>b[0])
      i = 1;
    for(q=0;q<3;q++)
      dSwap(&C->L[i][q],&C->L[2][q]);
  }      
}
//=========================================================================
void Lat_Align(Cell *C) // aligns 1st along x, 2nd in x-y
{
  int i;
  double a,t,r[3],L[3],sina,cosa;

  for(i=0;i<3;i++)
    r[i] = sqrt( C->L[i][1]*C->L[i][1] + C->L[i][2]*C->L[i][2] );
  a = 0.0;
  if(r[0]>1e-10)
    a = acos( C->L[0][1]/r[0] );
  if( C->L[0][2]<0.0 )
    a *= -1.0;
  sina = sin(a);
  cosa = cos(a);
  for(i=0;i<3;i++)
  {
    t = C->L[i][1];
    C->L[i][1] = t*cosa + C->L[i][2]*sina;
    C->L[i][2] =-t*sina + C->L[i][2]*cosa; 
  }

  for(i=0;i<3;i++)
    L[i] = VectorLen(C->L[i],3);
  a = 0.0;
  if(r[0]>1e-10)
    a = acos( C->L[0][0]/L[0] );
  if( C->L[0][1]<0.0 )
    a *= -1.0;
  sina = sin(a);
  cosa = cos(a);
  for(i=0;i<3;i++)
  {
    t = C->L[i][0];
    C->L[i][0] = t*cosa + C->L[i][1]*sina;
    C->L[i][1] =-t*sina + C->L[i][1]*cosa;
  }

  for(i=1;i<3;i++)
    r[i] = sqrt( C->L[i][1]*C->L[i][1] + C->L[i][2]*C->L[i][2] );
  a = 0.0;
  if(r[0]>1e-10)
    a = acos( C->L[1][1]/r[1] );
  if( C->L[1][2]<0.0 )
    a *= -1.0;
  sina = sin(a);
  cosa = cos(a);
  for(i=1;i<3;i++)
  {
    t = C->L[i][1];
    C->L[i][1] = t*cosa + C->L[i][2]*sina;
    C->L[i][2] =-t*sina + C->L[i][2]*cosa;
  }
}
//======================================================================
void Copy_C(Cell *C, Cell *C1)
{
  int i,j,q;
  C1->N = C->N;
  C1->P = C->P;
  C1->XT = C->XT;
  for(i=0;i<C->N;i++)
  {
    for(q=0;q<D3;q++)
    {
      C1->X[i][q] = C->X[i][q];
      C1->FF[i][q] = C->FF[i][q];
    }
    C1->ATMN[i] = C->ATMN[i];
    C1->ATMZ[i] = C->ATMZ[i];
    C1->Nn[i] = C->Nn[i];
    for(j=0;j<C1->Nn[i];j++)
      C1->Ni[i][j] = C->Ni[i][j];
    for(j=0;j<C1->Nn[i];j++)
      for(q=0;q<D3;q++)
	C1->S[i][j][q] = C->S[i][j][q];
  }
  C1->A    = C->A;
  C1->MNT  = C->MNT;
  C1->nspc   = C->nspc;
  C1->NSPC = C->NSPC;
  for(i=0;i<C->NSPC;i++)
    C1->SPCZ[i] = C->SPCZ[i];
  for(i=0;i<C->NSPC;i++)
    C1->mass[C->SPCZ[i]] = C->mass[C->SPCZ[i]];

  for(i=0;i<C->nspc;i++)
    C1->spcz[i] = C->spcz[i];

  for(q=0;q<D3;q++)
  {
    C1->min[q] = C->min[q];
    C1->max[q] = C->max[q];
    C1->ANG[q] = C->ANG[q];
    C1->LAT[q] = C->LAT[q];
  }
  for(i=0;i<D3;i++)
    for(q=0;q<D3;q++)
      C1->L[i][q] = C->L[i][q];
  for(i=0;i<D3;i++)
    for(q=0;q<D3;q++)
      C1->R[i][q] = C->R[i][q];
  sprintf(C1->TAG,"%s",C->TAG);

  if(C1->EVOK==1&&C->EVOK==0)
    exit(0);

  if(C1->EVOK==1&&C->EVOK==1)
    for(i=0;i<C->N*D3;i++)
    {
      C1->ev[i] = C->ev[i];
      for(j=0;j<C->N*D3;j++)
	C1->EV[i][j] = C->EV[i][j];      
    }
}
//------------------------------------------------------------------------- 
void Clone(Cell *C, Cell *C1, int N0, int N1, int N2) 
{ 
  int i,q0,q1,q2; 


  for(q0=0;q0<N0;q0++) 
  for(q1=0;q1<N1;q1++)  
  for(q2=0;q2<N2;q2++)  
    if(!(q0==0&&q1==0&&q2==0))
    ADD(C,C1,C->L[0][0]*(double)q0+C->L[1][0]*(double)q1+C->L[2][0]*(double)q2,
	     C->L[0][1]*(double)q0+C->L[1][1]*(double)q1+C->L[2][1]*(double)q2,
             C->L[0][2]*(double)q0+C->L[1][2]*(double)q1+C->L[2][2]*(double)q2); 
  for(i=0;i<3;i++)
  {
    C->L[0][i] *= (double)N0;
    C->L[1][i] *= (double)N1; 
    C->L[2][i] *= (double)N2; 
  }
  ORDER(C);
  return;
}
//=======================================================================================
void abc(Cell *C)
{
  int i;

  for(i=0;i<3;i++)
    C->LAT[i] = VectorLen(C->L[i],3);
  for(i=0;i<3;i++)
    C->ANG[i] = acos( DotProd(C->L[(i+1)%3],C->L[(i+2)%3],3)/(C->LAT[(i+1)%3]*C->LAT[(i+2)%3]));
}
//-------------------------------------------------------------------------
void ABC_LT(Cell *C)
{
  double t = C->LAT[2]*(cos(C->ANG[0])-cos(C->ANG[1])*cos(C->ANG[2]))/sin(C->ANG[2]);

  C->L[0][0] = C->LAT[0];                  C->L[0][1] = 0.0;                        C->L[0][2] = 0.0;
  C->L[1][0] = C->LAT[1]*cos(C->ANG[2]);   C->L[1][1] = C->LAT[1]*sin(C->ANG[2]);   C->L[1][2] = 0.0;
  C->L[2][0] = C->LAT[2]*cos(C->ANG[1]);   C->L[2][1] = t;                          C->L[2][2] = sqrt( pow(C->LAT[2]*sin(C->ANG[1]),2.0)-t*t );
}
//-------------------------------------------------------------------------
void KMESHOLD(Cell *C, int N, char name[])
{
  int i,q,K[3];
  double a[3],b[3];
  FILE *out;

  if( N==1 )
    K[0] = K[1] = K[2] = 1;
  else
  {
    N = (int)( (double)N/(double)C->N );
    for(i=0;i<3;i++)
      a[i] = VectorLen(C->L[i],3);
    
    for(i=1,q=0;i<3;i++)
      if(a[i]<(a[q]-1e-10))
	q = i;
    
    for(q=0;q<3;q++)
    {
      b[q] = pow( (double)N*(a[(q+1)%3]*a[(q+2)%3])/(a[q]*a[q]),1.0/3.0);
      K[q] = floor(b[q]+0.5);
      if(K[q]<2 && a[q] < 25.0 )
	K[q] = 2;
    }
  }
  
  out = fopen(name,"w");
  fprintf(out,"KPOINTS file\n0\n");
  if( fabs(a[0]-a[1])<1e-10 && fabs(COS(C->L[0],C->L[1],3)+0.5)<1e-10 )
    fprintf(out,"G\n");
  else
    fprintf(out,"M\n");
  fprintf(out,"%2d %2d %2d\n 0  0  0\n",K[0],K[1],K[2]);
  fclose(out);
}
//-------------------------------------------------------------------------
void KMESH(Cell *C, double KM, char name[], int ND)
{
  int i,q,K[3];
  double a[3];
  FILE *out;

  Reciprocal(C);

  for(i=0;i<3;i++)
    a[i] = VectorLen(C->R[i],3);
  for(q=0;q<ND;q++)
    K[q] = ceil(a[q]/KM);
  for(q=ND;q<3;q++)
    K[q] = 1;
  out = fopen(name,"w");
  fprintf(out,"KPOINTS file\n0\n");
  if( fabs(a[0]-a[1])<1e-10 && fabs(COS(C->L[0],C->L[1],3)+0.5)<1e-10 )
    fprintf(out,"G\n");
  else
    fprintf(out,"M\n");
  fprintf(out,"%2d %2d %2d\n 0  0  0\n",K[0],K[1],K[2]);
  fclose(out);
}
//-------------------------------------------------------------------------
int Check_OUTCAR(Cell *C, char file[200])
{
  int  J,E;
  char buf[200];
  FILE *in;

  if((in=fopen(file,"r")))
    fclose(in);
  else
    return -1;
  sprintf(buf,"grep 'al CPU' %s | grep -c CPU",file);
  in = popen(buf,"r");
  fscanf(in,"%d",&J);
  pclose(in);
  sprintf(buf,"grep ERROR %s | grep -c ERROR",file);
  in = popen(buf,"r");
  fscanf(in,"%d",&E);
  pclose(in);
  if(E>0)
    E = 1;
  return J*(1-E);
}
//-------------------------------------------------------------------------
int Check_GULP_OUT(Cell *C, char file[200])
{
  int J;
  char buf[200];
  FILE *in;

  if((in=fopen(file,"r")))
    fclose(in);
  else
    return -1;
  sprintf(buf,"grep 'Job Finished' %s | grep -c Job ",file);
  in=popen(buf,"r");
  fscanf(in,"%d",&J);
  pclose(in);
  return J;
}
//-------------------------------------------------------------------------
void atomSwap(Cell *C, int i, int j, int k)
{
  int q;

  iSwap(&C->Ni[i][j],&C->Ni[i][k]);
  for(q=0;q<D3;q++)
    dSwap(&C->S[i][j][q],&C->S[i][k][q]);
}
