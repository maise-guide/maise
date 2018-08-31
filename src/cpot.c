#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cdef.h"
#include "edef.h"
#include "ndef.h"
#include "nmlp.h"
#include "nmod.h"
#include "nutl.h"
#include "util.h"
#include "cell.h"
#include "cpot.h"
#include "nmin.h"

//==================================================
double dFc(double Rc, double r)
{
  if(r>Rc) return 0.0;
  return -exp(-r/(Rc-r))*Rc/((Rc-r)*(Rc-r));
}
//==================================================
double Fc(double Rc, double r)
{
  if(r>Rc) return 0.0;
  return exp(-r/(Rc-r));
}
//==================================================
double dFC(double R1, double R2, double r)
{
  if(r<R1||r>R2) return 0.0;
  return -0.5*sin( Pi*(r-R1)/(R2-R1))*Pi/(R2-R1);
}
//==================================================
double FC(double R1, double R2, double r)
{
  if(r<R1) return 1.0;
  if(r>R2) return 0.0;
  return 0.5*(1.0 + cos( Pi*(r-R1)/(R2-R1)));
}
//==================================================
double ENE_GP3D(Cell *C, int J)
{
  int i,j,p;
  double r,e,t,u,y;
  double fc(double Rc, double r);

  if(J==1)
    LIST(C);

  for(i=0,e=0.0;i<C->N;i++)
    {
      for(j=0,t=0.0;j<C->Nn[i];j++)
	{
	  p = C->ATMN[i]*C->ATMN[i]+C->ATMN[C->Ni[i][j]]*C->ATMN[C->Ni[i][j]];
	  r  = NDX(C,i,j);
	  u  =  FC(C->WW[p][5],C->WW[p][6],r);
	  y  = r/C->WW[p][4] - 1.0;
	  //printf("%lf %lf %lf % lf\n",r,u,y,C->WW[0][4]);
	  e += exp(-y*C->WW[p][1])*u*C->WW[p][0];
	  t += exp(-y*C->WW[p][3]*2.0)*u*C->WW[p][2]*C->WW[p][2];
	}
      e -= sqrt(t);
    }
  return e + 0.26*(double)C->N;
}
//==================================================    
double ENE_GP(Cell *C, int J)
{
  int i,j,p;
  double r,e,t,y;
  double fc(double Rc, double r);

  if(C->ND>0)
    return ENE_GP3D(C,1);

  for(i=0,e=0.0;i<C->N;i++)
  {
    for(j=0,t=0.0;j<C->N;j++)
    if(i!=j)
    {
      p = C->ATMN[i]*C->ATMN[i]+C->ATMN[j]*C->ATMN[j];
      r  = NDR(C,i,j);
      y  = r/C->WW[p][4] - 1.0;
      e += exp(-y*C->WW[p][1])*C->WW[p][0];
      t += exp(-y*C->WW[p][3]*2.0)*C->WW[p][2]*C->WW[p][2];
    }
    e -= sqrt(t);
  }
  return e;
}
//==================================================
// for NPs only and no cutoff
//==================================================
double ENE_SC(Cell *C, int J)
{
  int i,j;
  double r,e,t;
  double fc(double Rc, double r);

  for(i=0,e=0.0;i<C->N;i++)
  {
    for(j=0,t=0.0;j<C->N;j++)
    if(i!=j)
    {
      r  = NDR(C,i,j);
      e += pow(C->WW[0][0]/r,C->WW[0][1])*0.5;
      t += pow(C->WW[0][0]/r,C->WW[0][3]);
    }
    e -= sqrt(t)*C->WW[0][2];
  }
  return e*C->WW[0][4];
}
//==================================================        
void frc_sc(Cell *C)
{
  int i,q;
  double t,dx,Em,Ep;
  
  dx = 0.0001;
  for(i=0,Em=Ep=0.0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      t = C->X[i][q];
      C->X[i][q] += dx;
      if(C->MODT==3) Ep = ENE_SC(C,1);
      if(C->MODT==4) Ep = ENE_GP(C,1);
      C->X[i][q] = t - dx;
      if(C->MODT==3) Em = ENE_SC(C,1);
      if(C->MODT==4) Em = ENE_GP(C,1);
      C->F[i][q] = -(Ep-Em)/(dx*2.0);
      C->X[i][q] = t;
    }
}
//==================================================
double FRC_GP3D(Cell *C, int J)
{
  int i,j,q,p;
  double r,y,t,e,*h,u,z;
  double  fc(double Rc, double r);
  double dfc(double Rc, double r);

  h = make_d1D(C->N);
  if(J==1)
    LIST(C);

  for(i=0,e=0.0;i<C->N;i++)
  {
    for(j=0,t=0.0;j<C->Nn[i];j++)
    {
      p  = C->ATMN[i]*C->ATMN[i]+C->ATMN[C->Ni[i][j]]*C->ATMN[C->Ni[i][j]];
      r  = NDX(C,i,j);
      u  = FC(C->WW[p][5],C->WW[p][6],r);
      y  = r/C->WW[p][4] - 1.0;
      e += exp(-y*C->WW[p][1])*u*C->WW[p][0];      
      t += exp(-y*C->WW[p][3]*2.0)*u*C->WW[p][2]*C->WW[p][2];
    }
    h[i] = sqrt(t);
    e -= h[i];
  }

  if(J!=0)
    for(q=0;q<6;q++)
      C->U[q] = 0.0;
  for(i=0;i<C->N;i++)
  {
    if(J!=0)
      for(q=0;q<3;q++)
        C->F[i][q] = 0.0;
    for(j=0;j<C->Nn[i];j++)
    {
      p = C->ATMN[i]*C->ATMN[i]+C->ATMN[C->Ni[i][j]]*C->ATMN[C->Ni[i][j]];
      r = NDX(C,i,j);
      u = FC(C->WW[p][5],C->WW[p][6],r);
      z = dFC(C->WW[p][5],C->WW[p][6],r);
      y = r/C->WW[p][4] - 1.0;
      
      t = ( 2.0*  C->WW[p][0]*            (     u*C->WW[p][1]/C->WW[p][4]-z)*exp(-y*C->WW[p][1])
           +0.5*( C->WW[p][2]*C->WW[p][2]*(-2.0*u*C->WW[p][3]/C->WW[p][4]+z)*exp(-y*C->WW[p][3]*2.0)*(1.0/h[i]+1.0/h[C->Ni[i][j]]) ))/r;

      for(q=0;q<3;q++)
        C->F[i][q] -= t*DX(C,i,j,q);
      for(q=0;q<3;q++)
        C->U[q]   += 0.5*t*DX(C,i,j,q)*DX(C,i,j,q);
      for(q=0;q<3;q++)
        C->U[q+3] += 0.5*t*DX(C,i,j,q)*DX(C,i,j,(q+1)%3);
    }
  }

  free_d1D(h);
  return e;
}
//==================================================
double FRC_GP(Cell *C, int J)
{
  int i,j,q,p;
  double r,y,t,e,*h;
  double  fc(double Rc, double r);
  double dfc(double Rc, double r);

  if(C->ND>0)
    return FRC_GP3D(C,J);

  h = make_d1D(C->N);

  for(i=0,e=0.0;i<C->N;i++)
  {
    for(j=0,t=0.0;j<C->N;j++)
    if(i!=j)
    {
      p = C->ATMN[i]*C->ATMN[i]+C->ATMN[j]*C->ATMN[j];
      r  = NDR(C,i,j);
      y  = r/C->WW[p][4] - 1.0;
      e += exp(-y*C->WW[p][1])*C->WW[p][0];      
      t += exp(-y*C->WW[p][3]*2.0)*C->WW[p][2]*C->WW[p][2];
    }
    h[i] = sqrt(t);
    e -= h[i];
  }

  for(i=0;i<C->N;i++)
  {
    for(q=0;q<3;q++)
      C->F[i][q] = 0.0;
    for(j=0;j<C->N;j++)
    if(i!=j)
    {
      p = C->ATMN[i]*C->ATMN[i]+C->ATMN[j]*C->ATMN[j];
      r = NDR(C,i,j);
      y = r/C->WW[p][4] - 1.0;
      t = ( 2.0*  C->WW[p][0]*            (     C->WW[p][1]/C->WW[p][4])*exp(-y*C->WW[p][1])
	   +0.5*( C->WW[p][2]*C->WW[p][2]*(-2.0*C->WW[p][3]/C->WW[p][4])*exp(-y*C->WW[p][3]*2.0)*(1.0/h[i]+1.0/h[j]) ))/r;
      for(q=0;q<3;q++)
	C->F[i][q] -= t*(C->X[j][q]-C->X[i][q]);
    }
  }

  free_d1D(h);
  return e;
}
//==================================================        
// for NPs only and no cutoff
//==================================================
double FRC_SC(Cell *C, int J)
{
  int i,j,q;
  double r,t,e,*h;
  double  fc(double Rc, double r);
  double dfc(double Rc, double r);

  h = make_d1D(C->N);

  for(i=0,e=0.0;i<C->N;i++)
  {
    for(j=0,t=0.0;j<C->N;j++)
    if(i!=j)
    {
      r  = NDR(C,i,j);
      e += pow(C->WW[0][0]/r,C->WW[0][1])*0.5;
      t += pow(C->WW[0][0]/r,C->WW[0][3]);
    }
    h[i] = sqrt(t);
    e -= h[i]*C->WW[0][2];
  }
  
  for(i=0;i<C->N;i++)
  {
    for(q=0;q<3;q++)
      C->F[i][q] = 0.0;
    for(j=0;j<C->N;j++)
    if(i!=j)
    {
      r =   NDR(C,i,j);
      t = C->WW[0][4]*( (C->WW[0][1]/r)*pow(C->WW[0][0]/r,C->WW[0][1]) - 0.5*C->WW[0][2]*(C->WW[0][3]/r)*pow(C->WW[0][0]/r,C->WW[0][3])*(1.0/h[i]+1.0/h[j]) )/r;
      for(q=0;q<3;q++)
	C->F[i][q] -= t*(C->X[j][q]-C->X[i][q]);
    }
  }
  free_d1D(h);
  return e*C->WW[0][4];
}
//==================================================
double ENE_GP_(ANN *R, LNK *L)
{
  int i,j,p;
  double r,e,t,u,y;

  p = 0;
  for(i=0,e=0.0;i<L->N;i++)
  {
    for(j=0,t=0.0;j<R->D;j++)
      if(L->Cn[i][j]>1e-14)
      {
	r  = L->Cn[i][j];
	u  =  FC(R->WW[p][5],R->WW[p][6],r);
	y  = r/R->WW[p][4] - 1.0;
	e += exp(-y*R->WW[p][1])*u*R->WW[p][0];
	t += exp(-y*R->WW[p][3]*2.0)*u;
      }
    e -= sqrt(t)*R->WW[0][2]; // for single species only
  }
  return e;
}
//==================================================
double ENE_SC_(ANN *R, LNK *L)
{
  int i,j,k,p,n,N;
  double r,e,t,u;
  double fc(double Rc, double r);

  N = R->D/R->NSPC;

  for(i=0,e=0.0;i<L->N;i++)
  {
    for(k=0,t=0.0;k<R->NSPC;k++)
      for(n=0;n<N;n++)
      {
	j  = k*N + n;
	if(L->Cn[i][j]>1e-14)
        {
	  p  = 2*L->ATMN[i]*(1-k)+k;
	  r  = L->Cn[i][j];
	  u  = fc(R->WW[p][4],r);
	  e += pow(R->WW[p][0]/r,R->WW[p][1])*u*0.5;
	  t += pow(R->WW[p][2]/r,R->WW[p][3])*u;
	}
      }
    e -= sqrt(t);
  }
  return e;
}
//==================================================
//  Elemental or stratified training
//==================================================
double DE_SC_(ANN *R, LNK *L)
{
  int i,j,k,n,p,N;
  double x,y,E0,e,c0,u,w,v;
  double fc(double Rc, double r);
  double Wp[4];

  for(i=0;i<R->NW;i++)
    Wp[i] = 0.0;

  N = R->D/R->NSPC;
  for(i=0,E0=0.0;i<L->N;i++)
  {
    for(k=0,e=w=v=0.0;k<R->NSPC;k++)
      for(n=0;n<N;n++)
      {
        j  = k*N + n;
	if(L->Cn[i][j]>1e-14)
        {
	  p  = 2*L->ATMN[i]*(1-k)+k;
	  u = fc(R->WW[p][4],L->Cn[i][j]);
	  x = pow(R->WW[p][0]/L->Cn[i][j],R->WW[p][1]);
	  y = pow(R->WW[p][2]/L->Cn[i][j],R->WW[p][3]);

	  if(p==R->NSPC-1)
	  {
	    Wp[0] += 0.5*x*u*R->WW[p][1]/R->WW[p][0];
	    Wp[1] += 0.5*x*u*log(R->WW[p][0]/L->Cn[i][j]);
	    w     +=     y*u*log(R->WW[p][2]/L->Cn[i][j]);
	    v     +=     y*u;
	  }
	  e     +=     y*u;
	  E0    += 0.5*x*u;
	}
      }
    if(e>1e-14)
    {
      e = sqrt(e);
      E0    -= e;
      Wp[2] -= 0.5*v/e*R->WW[R->NSPC-1][3]/R->WW[R->NSPC-1][2];
      Wp[3] -= 0.5*w/e;
    }
  }

  c0 = 2.0*(E0-L->E)/(double)(L->N*L->N)*R->WE;

  for(n=0;n<R->NW;n++)
    R->WWp[R->NSPC-1][n] += Wp[n]*c0;

  return E0;

}
//==================================================
double ENE_LJ(Cell *C)
{
  int i,j;
  double r,x,a,A,B,E,e;
  double fc(double Rc, double r);

  LIST(C);
  a = C->LJa;
  E = C->LJe;
  A = -2.0*E*pow(a,6.0);
  B =  0.5*A*pow(a,6.0);

  for(i=0,e=0.0;i<C->N;i++)
    for(j=0;j<C->Nn[i];j++)
    {
      r = NDX(C,i,j);
      x = pow(r,-6.0);
      e += x*(-A + B*x)*fc(C->Rc,r);
    }
  return e;
}
//==================================================
void frc_lj(Cell *C)
{
  int i,q;
  double t,dx,Em,Ep;

  dx = 0.001;
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      t = C->X[i][q];
      C->X[i][q] += dx;
      Ep = ENE_LJ(C);
      C->X[i][q] = t - dx;
      Em = ENE_LJ(C);
      C->F[i][q] = -(Ep-Em)/(dx*2.0);
      C->X[i][q] = t;
    }
  
}
//==================================================
double FRC_LJ(Cell *C)
{
  int i,j,q;
  double r,x,a,A,B,E,t,e;
  double  fc(double Rc, double r);
  double dfc(double Rc, double r);

  LIST(C);
  a = C->LJa;
  E = C->LJe;
  A = -2.0*E*pow(a,6.0);
  B =  0.5*A*pow(a,6.0);

  for(i=0,e=0.0;i<C->N;i++)
  {
    for(q=0;q<3;q++)
      C->F[i][q] = 0.0;
    for(j=0;j<C->Nn[i];j++)
    {
      r =  NDX(C,i,j);
      x =  pow(r,-6.0);
      t =  x*(6.0*A - 12.0*B*x)/r*fc(C->Rc,r) + x*(-A + B*x)*dfc(C->Rc,r);
      e += x*(-A + B*x)*fc(C->Rc,r);
      for(q=0;q<3;q++)
	C->F[i][q] += 2.0*t*DX(C,i,j,q)/r;
    }
  }
  return e;
}
//===================================================
void SAVE_SC(ANN *R, double time)
{
  int i,iii,n,m;
  char s[200];

  FILE* out;

  sprintf(s,"%s/%s",R->otpt,"sc.dat");
  out=fopen(s,"w");
  
  for(n=0;n<R->NSPC;n++)
  {
    fprintf(out,"MODEL %3d %3d\n",R->MODT,R->SPCZ[n]);
    fprintf(out,"%d\n",R->NW);
    for(i=0;i<R->NW+1;i++)       // the last two are not adjusted                                                                                                                                                                         
      fprintf(out,"% 24.16lf\n",R->WW[2*n*n][i]);
    for(m=n+1;m<R->NSPC;m++)
    {
      fprintf(out,"MODEL %3d %3d %3d\n",R->MODT,R->SPCZ[n],R->SPCZ[m]);
      fprintf(out,"%d\n",R->NW);
      for(i=0;i<R->NW+1;i++)       // the last two are not adjusted                                                                                                                                                               
	fprintf(out,"% 24.16lf\n",R->WW[n*n+m*m][i]);      
    }
  }
  fprintf(out,"STR     %d\n",R->STR);
  fprintf(out,"N       %d\n",R->N);
  fprintf(out,"T       %d\n",R->TN);
  fprintf(out,"MAX_NA  %d\n",R->A);
  fprintf(out,"NW      %d\n",R->NW);
  fprintf(out,"NSPC    %d\n",R->NSPC);
  fprintf(out,"TSPC   ");for(iii=0;iii<R->NSPC;iii++) fprintf(out," %d",R->SPCZ[iii]); fprintf(out,"\n");
  fprintf(out,"--------------------------------------------------------\n");
  fprintf(out,"Total training time=       %24.16lf (s)   N of Iterations=%d      Optimizer=%d\n",time,R->ITER,R->MINT);
  fprintf(out,"Avg. Training Error=       %24.16lf      Avg. Test Error=%24.16lf\n",R->RE,R->EE);
  fprintf(out,"Average Energy(all strct.)=%24.16lf      Stndrd. Dev.=   %24.16lf\n",R->Eavg,R->Edev);
  fclose(out);
}
//===================================================
// only one species for now
//===================================================
void READ_MOD(Cell *C, char *file)
{
  int  k,n,m,i,Z[3],nw,z[3];
  FILE *in;
  char s[200],buf[200],f[200];

  if(C->MODT==0)
    return;

  n = 0;
  if(C->MODT==3)
  {
    sprintf(f,"%s/sc.dat",file);
    nw = 2;
  }
  if(C->MODT==4)
  {
    sprintf(f,"%s/gp.dat",file);
    nw = 2;
  }
  in = fopen(f,"r");
  if( in!= NULL )
  {
    n = 0;
    while(fgets(buf,200,in))
      if(strncmp(buf,"MODEL",5)==0) 
      {
	Z[0] = Z[1] = 0;
	sscanf(buf+5,"%d %d %d",&m,&Z[0],&Z[1]);

	if(Z[1]==0)  Z[1] = Z[0];
	C->SPCZ[0]=Z[0];
	C->SPCZ[1]=Z[1];

	for(z[0]=0;z[0]<C->NSPC&&C->SPCZ[z[0]]!=Z[0];z[0]++);
	for(z[1]=0;z[1]<C->NSPC&&C->SPCZ[z[1]]!=Z[1];z[1]++);

	if( z[0]==C->NSPC || z[1]==C->NSPC )
	{
	  sprintf(buf,"Species in setup and %s are inconsistent: %d %d %d %d",f,Z[0],Z[1],z[0],C->SPCZ[z[0]]);
	  fprintf(stderr,"Species in setup and %s are inconsistent: %d %d %d %d",f,Z[0],Z[1],z[0],C->SPCZ[z[0]]);
	  EXIT(buf);
	}
	fgets(s,200,in);
	sscanf(s,"%d",&C->NW);
	for(i=0;i<C->NW+nw;i++)       // the last one(s) are not adjusted
        {
	  fgets(s,200,in);
	  sscanf(s,"%lf",&C->WW[z[0]*z[0]+z[1]*z[1]][i]);
	}    
	n++;
      }	
    fclose(in);
    if( !( (n==C->NSPC) || (C->NSPC==1&&n==1) || (C->NSPC==2&&n==3)|| (C->NSPC==3&&n==6) ) )
    {
      fprintf(stderr,"The number of species in setup and are incosistent\n");
      EXIT("The number of species in setup and are incosistent");
    }

    if(n==C->NSPC)       // if A-B are not specified set them as average
       for(k=0;k<7;k++)
	 for(n=0;n<C->NSPC;n++)
	   for(m=n+1;m<C->NSPC;m++)
	     C->WW[n*n+m*m][k] = 0.5*(C->WW[2*n*n][k]+C->WW[2*m*m][k]);

    for(n=0, C->Rc = C->rc = 0.0; n<C->NSPC; n++)
      if( C->Rc < C->WW[2*n*n][5] )
	C->Rc = C->rc = C->WW[2*n*n][C->NW+nw-1];
  }
  else
  {
    if(C->MODT==3)
      fprintf(stderr,"Please provide starting parameters for Sutton-Chen potential in %s\n",f);
    if(C->MODT==4)
      fprintf(stderr,"Please provide starting parameters for Gupta potential in %s\n",f);
  }

  if(0)
    for(n=0;n<C->NSPC;n++)
      for(m=n;m<C->NSPC;m++)
      {
	printf("%3d %3d   ",C->SPCZ[n],C->SPCZ[m]);
	printf("NW = %d nw = %d\n",C->NW,nw);
	for(k=0;k<C->NW+nw;k++)
	  printf("% lf",C->WW[n*n+m*m][k]);
	printf("\n");
      }

}
//===================================================
void INIT_MOD(ANN *R, Cell *C)
{
  int  i,p;
  FILE *in;
  char s[200];

  if(R->MODT==3)
    sprintf(s,"%s/sc.dat",R->otpt);
  if(R->MODT==4)
    sprintf(s,"%s/gp.dat",R->otpt);
  in = fopen(s,"r");
  if( in!= NULL )
  {
    fclose(in);
    READ_MOD(C,R->otpt);
  }
  R->nw = R->NW = C->NW;

  for(p=0;p<2*C->NSPC-1;p++)
    for(i=0;i<C->NW+1;i++)         
      R->WW[p][i] = C->WW[p][i];
}
//===================================================
double scfunc(ANN *R, LNK *L)
{
  int n;

  for(n=0,R->RE=R->RT=0.0;n<R->N;n++)
    R->RE += pow( (L[n].E-ENE_ANN(R,&L[n]))/(double)L[n].N,2.0 );
  R->RT = (R->RE*R->WE + R->RT*R->WF)/(double)R->N;

  return R->RT;
}
//===================================================
void scfunc_num(ANN *R, LNK *L)
{
  int i;
  double o,dx;

  dx = 0.000001;

  for(i=0;i<R->NW;i++)
  {
    o = R->WW[R->NSPC-1][i];
    R->WW[ R->NSPC-1][i] = o + dx;
    R->WWp[R->NSPC-1][i] = scfunc(R,L);
    R->WW[ R->NSPC-1][i] = o - dx;
    R->WWp[R->NSPC-1][i] = (R->WWp[R->NSPC-1][i]-scfunc(R,L))/(2.0*dx);
    R->WW[ R->NSPC-1][i] = o;
    printf("NUM %3d % 24.16lf % 24.16lf\n",i,R->WW[R->NSPC-1][i],R->WWp[R->NSPC-1][i]);
  }  
  return;
}
//===================================================
void scdfunc(ANN *R, LNK *L)
{
  int p,n;

  R->N = 1;
  scfunc_num(R,&L[0]);
  for(p=0;p<2*R->NSPC-1;p++)
    for(n=0;n<R->NW;n++)
      R->WWp[p][n] = 0.0;
  for(n=0;n<1;n++)
    DE_SC_(R,&L[n]);
  for(p=0;p<2*R->NSPC-1;p++)
    for(n=0;n<R->NW;n++)
      R->WWp[p][n] /= 1.0;
}
//===================================================
void TRAN_SC(ANN *R, Cell *C, LNK *L)
{
  int n,m,i;
  double E,ERR[100];

  FILE *out;
  char s[200];
  double WW[20];

  sprintf(s,"%s/out.dat",R->otpt);

  printf("D  %3d\n",R->D);
  for(i=0;i<1;i++)
  {
    Random();
    printf("Total number of parameters: %3d\n",R->NW);
    printf("ERR = % lf\n",sqrt(TOT_ERR(R,L)));
    
    out=fopen(s,"w");
    fprintf(out,"Total number of parameters: %3d\n",R->NW);
    fprintf(out,"ERR = % lf\n",sqrt(TOT_ERR(R,L)));
    fclose(out);

    printf("% lf\n",sqrt(scfunc(R,L)));
    if(0)
    {
      scdfunc(R,L);
      for(i=0;i<R->NW;i++)
	printf("ANA %3d % 24.16lf % 24.16lf\n",i,R->WW[R->NSPC-1][i],R->WWp[R->NSPC-1][i]);
      exit(0);
    }

    for(n=0;n<R->NW;n++)
      WW[n] = R->WW[0][n];
    E = TOT_ERR(R,L);
    if(R->MITR>0)
    {
      printf("%3d\n",R->N);
      if(R->MINT>=1)
        bfgs_mlp(R,L);
      else
	for(m=0;m<3;m++)
        {
	  if(m>0)
	    for(n=0;n<R->NW;n++)
	      R->WW[R->NSPC-1][n] *= (0.7 + 0.6*Random());
	  for( n=0;n<10;n++)
	    bfgs_mlp(R,L);
	  if(TOT_ERR(R,L)>E)
	  {
	    E = TOT_ERR(R,L);
	    for(n=0;n<R->NW;n++)
	      WW[n] = R->WW[R->NSPC-1][n];
	  }
	}
    }
    ERR[i] = sqrt(TOT_ERR(R,L));
    for(i=0,printf("\n");i<R->NW;i++)
      printf("ANA %3d % 24.16lf % 24.16lf\n",i,R->WW[R->NSPC-1][i],R->WWp[R->NSPC-1][i]);
  }
  return;
}
//===================================================
