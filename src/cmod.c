#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "cdef.h"
#include "ndef.h"
#include "edef.h"
#include "cutl.h"
#include "cell.h"
#include "cpot.h"
#include "nprs.h"
#include "nmlp.h"
#include "nmod.h"
#include "nutl.h"
#include "cmin.h"
#include "cfnc.h"
#include "util.h"

extern const double Pi;
extern const double eV2GPa;
const double kb = 11604.0;

//========================================================================
double KINETIC(Cell *C)
{
  int i,q;

  for(q=0,C->K=0.0;q<D3;q++)
    for(i=0;i<C->N;i++)
      C->K += C->m[C->ATMZ[i]]*C->V[q][i]*C->V[q][i];
  return (C->K *= 0.5);
}

//=======================================================================
double P_Lindemann(Cell *C, int J)
{
  int i,j;
  double r,LI;

  if(J==0)
  {
    for(i=0,C->LI=0;i<C->N;i++)
      for(j=0;j<C->Nn[i];j++)
	C->R1[i][j] = C->R2[i][j] = 0.0;
    return 0.0;
  }
  if(J==1)
  {
    LIST(C);
    C->LI++;
    LI = 0.0;
    for(i=0;i<C->N;i++)
      for(j=0;j<C->Nn[i];j++)
      {
	C->R2[i][j] += C->NDX[i][j]*C->NDX[i][j];
	C->R1[i][j] += C->NDX[i][j];
	if( (r=C->R2[i][j]*(double)C->LI - C->R1[i][j]*C->R1[i][j]) > 1e-12 )
	  LI += sqrt( r ) / C->R1[i][j]/ (double)C->Nn[i];
      }
    return LI/(double)(C->N);
  }
  return 0.0;
}
//=======================================================================
double Lindemann(Cell *C, int J)
{
  int i,j,q;
  double x,r,LI;

  if(J==0)
  {
    for(i=0,C->LI=0;i<C->N;i++)
      for(j=i+1;j<C->N;j++)
	C->R1[i][j] = C->R2[i][j] = 0.0;
    return 0.0;
  }
  if(J==1)
  {
    C->LI++;
    LI = 0.0;
    for(i=0;i<C->N;i++)
      for(j=i+1;j<C->N;j++)
      {
	for(q=0,r=0.0;q<D3;q++)
	{
	  x  = C->X[q][i] - C->X[q][j];
	  r += x*x;
	}
	C->R2[i][j] += r;
	C->R1[i][j] += sqrt(r);
	if( (r=C->R2[i][j]*(double)C->LI - C->R1[i][j]*C->R1[i][j]) > 1e-12 )
	  LI += sqrt( r ) / C->R1[i][j];
      }
    return 2.0*LI/(double)(C->N*(C->N-1));
  }
  return 0.0;
}
//=========================================================================
double CELL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{   
  int i;

  if( R->MODT==1 )
  {
    P->EFS = 0;
    PARS_STR(P,W,C,L,0,"."); 
    C->E = ENE_ANN(R,L);

    for(i=0;i<C->N;i++)
      C->E+=R->E0[C->ATMZ[i]];

    for(i=0;i<C->N;i++)
      C->EA[i] = L->EA[i]+R->E0[C->ATMZ[i]];
  }
  if( R->MODT>1 )
    C->E = ENE_POT(C);

  C->H = C->E;
  if(C->ND==3)
    C->H += C->p*Cell_VOLUME(C);
  return C->H;
}
//=========================================================================
double CELL_FRC(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int i,q;

  if( R->MODT==1 )
  {
    P->EFS = 3;
    PARS_STR(P,W,C,L,0,".");
    C->E = FRC_ANN(R,L);

    for(i=0;i<C->N;i++)
      C->E+=R->E0[C->ATMZ[i]];

    for(i=0;i<C->N;i++)
      C->EA[i] = L->EA[i]+R->E0[C->ATMZ[i]];
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->F[i][q] = L->f[i][q];
    for(q=0;q<6;q++)
      C->U[q] = L->s[q];
  }
  if( R->MODT>1 )
    C->E = FRC_POT(C);

  //=====  remove noise to preserve symmetry in some cases =====
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      if( fabs(C->F[i][q])<1e-8 )
	C->F[i][q] = 0.0;
  for(q=0;q<6;q++)
    if( fabs(C->U[q])<1e-8 )
      C->U[q] = 0.0;

  C->H = C->E;

  if(C->ND==3)
    C->H += C->p*Cell_VOLUME(C);
  return C->H;
}
//=========================================================================
void  MDOUT(Cell *C, int PIPE, int DISK, int NS, int n, int k, double E0, double T,double li)
    // PUTS the center of mass at the origin!
{
  FILE   *out;
  char   file[100];

  if(PIPE)
    printf("%8d:  P = % 24.14lf       K = % 24.14lf        dE = % 24.14lf     LI = % 24.16lf\n",NS*n+k,C->H/(double)C->N,C->K*kb/1.5/(double)C->N,(C->H+C->K-E0)/(double)C->N,li);
  if(DISK==0)
    return;
  sprintf(file,"MD/T%04d.dat",(int)(T*kb));
  if(n==0&&k==0)
    out = fopen(file,"w");
  else
    out = fopen(file,"a");
  fprintf(out,"%8d  % 14.8lf  % 14.8lf  % 14.8lf % 14.8lf\n",NS*n+k,C->H/(double)C->N,C->K*kb/1.5/(double)C->N,(C->H+C->K-E0)/(double)C->N,li);
  fclose(out);
  JAR(C);

}
//=========================================================================
void Maxwell(Cell *C, double T, int therm)
{
  int    i,q;
  double x,v,f,Q,V[3];

  if(therm==2)
    return;

  v = sqrt(3.0*T);
  for(q=0;q<3;q++)
    V[q] = 0.0;
  for(i=0;i<C->N;i++)
  {
    f = Random()*2.0*Pi;
    Q = Random()    *Pi;
    C->V[i][0] = v * sin(Q)*cos(f)/sqrt(C->m[C->ATMZ[i]]);
    C->V[i][1] = v * sin(Q)*sin(f)/sqrt(C->m[C->ATMZ[i]]);
    C->V[i][2] = v * cos(Q)       /sqrt(C->m[C->ATMZ[i]]);
    for(q=0;q<3;q++)
      V[q] += C->V[i][q];
  }
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      C->V[i][q] -= V[q]/(double)C->N;

  if(therm==0)   //for dynamic thermostat
    {
      for(i=0,x=0.0;i<C->N;i++)
	for(q=0;q<3;q++)
	  x += C->V[i][q]*C->V[i][q]*(double)C->m[C->ATMZ[i]];
      x = v/sqrt(x/(double)C->N);
      
      for(i=0;i<C->N;i++)
	for(q=0;q<3;q++)
	  C->V[i][q] *= x;
    }  
}
//=========================================================================
void Dynamics(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L, double T, int PIPE)
{
  int    i,k,q,n,NS,NI;
  double E0;

  double dt,a;

  T    = T/kb;
  dt   = R->DELT;
  a    = 0.0;
  NS   = 1;
  NI   = 50;

  for(i=0,C->K=0.0;i<C->N;i++)
    for(q=0;q<D3;q++)
      C->K += 0.5*C->V[i][q]*C->V[i][q]*C->m[C->ATMZ[i]];

  printf("T = % lf K\n",C->K*kb/1.5/(double)C->N);

  C->H = CELL_FRC(R,P,W,C,L);
  E0 = C->H + C->K;

  for( n=0;n<NI;n++)
    for(k=0, MDOUT(C,PIPE,1,NS,n,k,E0,T,0.0) ;k<NS;k++)
    {
      for(i=0;i<C->N;i++)
	for(q=0;q<3;q++)
	{
	  C->X[i][q] += C->V[i][q]*dt + 0.5*C->F[i][q]*dt*dt/C->m[C->ATMZ[i]];
	  C->V[i][q] +=  0.5*C->F[i][q]*dt/C->m[C->ATMZ[i]];      
	}
      C->H = CELL_FRC(R,P,W,C,L);
      for(i=0,C->K=0.0;i<C->N;i++)
	for(q=0;q<D3;q++)
	  C->K += 0.5*C->V[i][q]*C->V[i][q]*C->m[C->ATMZ[i]];
      for(i=0;i<C->N;i++)
	for(q=0;q<D3;q++)
	{
	  C->V[i][q] = (C->V[i][q]+0.5*C->F[i][q]*dt/C->m[C->ATMZ[i]])/(1.0+a);
	  C->F[i][q] += -a*C->V[i][q];
	}

    }
  MDOUT(C,PIPE,1,NS,n-1,k,E0,T,0.0);
}
//========================================================================= 
//========================================================================= 
void Thermostat(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L, double T, int PIPE) 
{ 
  int    i,k,q,n,NS,NI; 
  double E0; 

  double dt,fi,xi,v2,fs,tau; 
  FILE   *out; 
  char   file[100]; 
  double KE,PE;
  int sam=-1;
  double r,LI=0.0;
  int x=0;

  KE=PE=0.0;
  dt   = R->DELT;
  T    = T/kb;
  tau  = R->COPL;
  xi   = 0.0; 
  NS   = 10;
  if(R->NSTP<10)
    NS = 1;
  NI   = R->NSTP/NS;

  if(C->ND==0)  Lindemann(C,0); else P_Lindemann(C,0);
  
  C->H = CELL_ENE(R,P,W,C,L);

  for(i=0,C->K=0.0;i<C->N;i++)
    for(q=0;q<D3;q++)
      C->K += 0.5*C->V[i][q]*C->V[i][q]*C->m[C->ATMZ[i]];
  E0 = C->H + C->K; 

  KE=-C->K;
  PE=-C->H;

  if(PIPE)        
    printf("%d %d %d %lf\n",NI,NS,NI*NS,T); 
  
  C->H = CELL_FRC(R,P,W,C,L);

  for( n=k=0;n<NI;n++)   
  {
    C->H = CELL_ENE(R,P,W,C,L);
    for(i=0,C->K=0.0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->K += 0.5*C->V[i][q]*C->V[i][q]*C->m[C->ATMZ[i]];

    if(C->ND==0)  r =  Lindemann(C,1); 
    else          r =P_Lindemann(C,1);    

    MDOUT(C,PIPE,1,NS,n,0,E0,T,r);
    KE+=C->K;
    PE+=C->H;
    sam++;
    LI+=r;

    for(k=0;k<NS;k++)   
    {   
      x++;
      C->H = CELL_FRC(R,P,W,C,L);
      for(i=0,v2=0.0;i<C->N;i++)   
	for(q=0;q<3;q++)   
	{   
	  fi = C->F[i][q]/C->m[C->ATMZ[i]]-xi*C->V[i][q];   
	  v2 += C->V[i][q]*C->V[i][q]*C->m[C->ATMZ[i]];   
	  C->X[i][q] += C->V[i][q]*dt + 0.5*fi*dt*dt;
	  C->V[i][q] += fi*dt;	  
	}   
      fs = 1.0/(tau*tau)*(v2/( (3.0*(double)C->N+0.0)*T)-1.0 );   
      xi += dt*fs;          
      JAR(C);
      if(R->MOVI>0&&x%R->MOVI==0)
      {
	sprintf(file,"mkdir -p %s/movie%04d",C->WDIR,(int)(T*kb));
	system(file);
	C->H=CELL_ENE(R,P,W,C,L);
	sprintf(C->TAG,"%lf",C->H);
	sprintf(file,"%s/movie%04d/POSCAR-%d",C->WDIR,(int)(T*kb),x);
	SAVE_CELL(C,file,0);
      }//for MD movie
    }
  }      
  C->H = CELL_ENE(R,P,W,C,L);
  for(i=0,C->K=0.0;i<C->N;i++)
    for(q=0;q<D3;q++)
      C->K += 0.5*C->V[i][q]*C->V[i][q]*C->m[C->ATMZ[i]];
  
  if(C->ND==0)  r=Lindemann(C,1); else r=P_Lindemann(C,1);    
  MDOUT(C,PIPE,1,NS,n-1,k,E0,T,r);KE+=C->K;PE+=C->H;sam++;
  LI+=r;
  
  sprintf(file,"%s/ave.dat",C->WDIR);
  out=fopen(file,"a");
  
  fprintf(out,"% lf % lf % lf % lf\n",T*kb,PE/(double)sam/(double)C->N,KE*kb/1.5/(double)C->N/(double)sam,LI/(double)sam);
  fclose(out);  
  printf("% lf % lf % lf % lf\n",T*kb,PE/(double)sam/(double)C->N,KE*kb/1.5/(double)C->N/(double)sam,LI/(double)sam);
  
  sprintf(file,"%s/CONTCAR_%04d",C->WDIR,(int)(kb*T));
  SAVE_CELL(C,file,0);
} 
//=========================================================================  
//=========================================================================
void CELL_MD(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int n,N;
  double Tmin, Tmax, dT, T;

  Tmin = R->TMIN;
  Tmax = R->TMAX;
  dT   = R->TSTP;
  
  system("mkdir -p MD");
  N = (int)fabs((Tmax-Tmin)/dT);
  Maxwell(C,Tmin/kb,R->THRM);
  C->POS = 1;
  P->EFS = 1;
  for(n=0;n<=N;n++)
  {
    T = Tmin + (double)n*dT;
    if(R->THRM==0) Dynamics(R,P,W,C,L,T,1);
    if(R->THRM> 0) Thermostat(R,P,W,C,L,T,1);
  }
}
//======================================================
void CELL_OUT(Cell *C)
{
  int    i;
  double V,a;
  char   s[200],b[200],t[200];
  FILE   *out;
  time_t rawtime;
  struct tm * timeinfo;

  a   = 180.0/Pi;

  sprintf(s,"OUTCAR");
  out = fopen(s,"a");
  
  abc(C);
  V = Cell_VOLUME(C);
  sprintf(b,"----------------------------------------------------------------------------------------------------------\n");

  if( C->it==0 )
  {
    time(&rawtime);
    timeinfo = localtime ( &rawtime );
    sprintf(t," %s",asctime(timeinfo));
    t[strlen(t)-1] = 0;
    fprintf(out,"%s    Cell optimization    with %s   model %s   on %s                     \n%s",b,C->VER,C->ID,t,b);
    sprintf(t," %d",C->MINT);
    fprintf(out,"   optimizer type         %s\n",t);
    sprintf(t," %d",C->RLXT);
    fprintf(out,"   relaxation type        %s\n",t);
    sprintf(t," %1.2lf GPa",C->p*eV2GPa);
    fprintf(out,"   target pressure        %s\n",t);

    sprintf(t," %d",C->ND);
    fprintf(out,"   dimensionality         %s\n",t);
    for(i=0,t[0]=0,sprintf(t+strlen(t)," ");i < C->NSPC;i++,sprintf(t+strlen(t),"    "))
      atom_symb(C->SPCZ[i],t+strlen(t));
    fprintf(out,"   species                %s\n",t);
    sprintf(  t,"   atoms of each species                     ");
    for(i=0;i < C->NSPC;i++)
      sprintf(t+27+6*i,"%d           ",C->SPCN[i]);
    fprintf(out,"%s\n",t);
    sprintf(t," %d",C->N);
    fprintf(out,"   total number of atoms  %s\n",t);
    fprintf(out,"%s\n\n\n",b);
  }  
  fprintf(out,"%s                   LATTICE CONSTANTS (Angst) AND ANGLES (degrees)                  VOLUME (Angst^3/atom) \n%s",b,b);
  fprintf(out," %3d % 12.6lf % 12.6lf % 12.6lf   % 12.6lf % 12.6lf % 12.6lf % 17.12lf vol\n",C->it,C->LAT[0],C->LAT[1],C->LAT[2],C->ANG[0]*a,C->ANG[1]*a,C->ANG[2]*a,V/(double)C->N);
  if(C->OUT/10>0)
  {
    fprintf(out,"%s                 POSITION (Angst)                     TOTAL-FORCE (eV/Angst)           ATOM ENERGY (eV)\n%s",b,b);
    for(i=0;i<C->N;i++)
      fprintf(out,"   % 12.6lf % 12.6lf % 12.6lf   % 12.6lf % 12.6lf % 12.6lf % 17.12lf\n",C->X[i][0],C->X[i][1],C->X[i][2],C->F[i][0],C->F[i][1],C->F[i][2],C->EA[i]);
    fprintf(out,"%s",b);
    fprintf(out,"  Total   ");
    for(i=0;i<6;i++)
      fprintf(out,"% 15.8lf ",C->U[i]/Cell_VOLUME(C));
    fprintf(out,"\n");
    fprintf(out,"  in kB   ");
    for(i=0;i<6;i++)
      fprintf(out,"% 15.8lf ",C->U[i]*eV2GPa*10.0/Cell_VOLUME(C));
    fprintf(out,"\n%s",b);
  }
  fprintf(out,"  iter %3d   total enthalpy= % 14.8lf   energy=  % 14.8lf   % 14.8lf  % 14.8lf\n%s\n\n\n",C->it,C->H,C->E,C->H/(double)C->N,C->E/(double)C->N,b);

  fclose(out);
  sprintf(s,"OSZICAR");
  out = fopen(s,"a");
  fprintf(out,"                            % 1.8lf\n",C->H);
  fclose(out);
}
//======================================================
int CELL_OK(Cell *C, double *Rm)
{
  int i;

  for(i=0;i<C->N;i++)
    if(NDX(C,i,0)<Rm[C->ATMZ[i]]+Rm[C->ATMZ[C->Ni[i][0]]])
      return 0;
  return 1;
}
//======================================================
int STOP_OK(Cell *C, double *Rm, double x)
{
  int i;

  for(i=0;i<C->N;i++)
    if(NDX(C,i,0)<x*Rm[C->ATMZ[i]]+Rm[C->ATMZ[C->Ni[i][0]]])
      return 0;
  return 1;
}
//======================================================
void CELL_RELX(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
    int    i;
    double H;
    FILE  *out;
    double time;

    if(C->RLXT == 0 && R->MITR > 0)
    {
      printf("Error: type of relaxation RLXT not specified in setup file.\n");
      exit(1);
    }
    
    system("rm -f OUTCAR OSZICAR");
    C->stop = 0;

    time=cpu_time();
    H = CELL_FRC(R,P,W,C,L);

    if(R->MITR==0)
      printf("%5d %24.16lf\n",0,H);
    if(C->OUT%10>0)
      CELL_OUT(C);
    if(R->MITR>0)
    {
      if(R->MINT>=0 && R->MINT<4)
      {
       	CELL_MIN(R,P,W,C,L);
      }
      else
      {
        fprintf(stderr,"ERROR: Enter valid value for MINT (0-3)\n");
        exit(1);
      }
       
      if( (C->OUT%10==1) )
	CELL_OUT(C);
    }
    if(C->OUT%10==0)
      CELL_OUT(C);
    LIST(C);
    if(!CELL_OK(C,R->Rm))
      system("echo >> OUTCAR;echo ERROR distances are too short >> OUTCAR");

    out=fopen("OUTCAR","a");
    fprintf(out,"  \n Total CPU time used (sec):  % 14.8lf \n",cpu_time()-time);
    fclose(out);

    return;
    for(i=0;i<C->N;i++)
      C->Nn[i] = 12;
    Print_List(C);
}
//======================================================
//
//======================================================
void CELL_TEST(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int    i,q,OK;
  double t,dx,Em,Ep,**F;

  dx = 0.00001;
  F  = make_d2D(C->N,3);

  OK = 1;
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      t = C->X[i][q];
      C->X[i][q] += dx;
      Ep = CELL_ENE(R,P,W,C,L);
      C->X[i][q] -= 2.0*dx;
      Em = CELL_ENE(R,P,W,C,L);
      F[i][q] = -(Ep-Em)/(dx*2.0);
      C->X[i][q] += dx;
    }
  CELL_FRC(R,P,W,C,L);

  for(i=0;i<C->N;i++,printf("\n"))
  {
    for(q=0,printf("NUM  ");q<3;q++)
      printf("% 24.12lf ",F[i][q]);
    printf("\n");
    for(q=0,printf("ANA  ");q<3;q++)
      printf("% 24.12lf ",C->F[i][q]);
    printf("\n");
    for(q=0,printf("DIF  ");q<3;q++)
      printf("% 24.12lf ",C->F[i][q]-F[i][q]);
    printf("\n");
    for(q=0;q<3;q++)
      if( fabs(C->F[i][q]-F[i][q])>1e-6 )
      	OK = 0;
  }
  if( OK==0 )
    printf("CHECK FAILED\n");
  else
    printf("CHECK PASSED\n");
  free_d2D(F,C->N);
}
//======================================================
//
//======================================================
void CELL_MAIN(ANN *R, PRS *P, Cell *C)
{
  int   i;
  PRS   W[9];
  LNK   L;

  if(R->NSPC==0)
  {
    fprintf(stderr,"Error: species are not defined in setup file!\n");
    exit(1);
  }
  Build_Cell(C,1);

  L.B = 0;
  Build_LNK(&L,C->N,C->NM,P->D,3);

  if(R->MODT==1)
  {
    Build_ANN(R);
    READ_ANN(R);
    Build_PRS(P,W,0);
  }
  else
    READ_POT(C,".");

  READ_CELL(C,"POSCAR");
  if( C->ND < 0 )
  {
    fprintf(stderr,"Error: please specify NDIM in setup file\n");
    exit(1);
  }

  P->IO = 0;
  for(i=0;i<C->N;i++)
    L.MRK[i]=1;
  P->EFS = 0;

  JAR(C);
  sprintf(C->ID,"%s",R->ID);
  C->MINT = R->MINT;

  if( R->JOBT==20 ) CELL_RELX(R,P,W,C,&L);
  if( R->JOBT==21 ) CELL_MD(R,P,W,C,&L);
  if( R->JOBT==22 ) CELL_TEST(R,P,W,C,&L);

  JAR(C);

  SAVE_CELL(C,"CONTCAR",0);
}
//======================================================
