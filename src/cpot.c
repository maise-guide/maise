#include "cpot.h"

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
//
//==================================================    
double ENE_POT(Cell *C)
{
  int    i,j,p;
  double r,e,t,y,E,u;

  LIST(C,0);
  E = 0.0;

#pragma omp parallel private(p,r,u,e,y,t,j) num_threads(C->NP)
  {
#pragma omp for reduction (+:E) schedule(dynamic,C->NB)
  for(i=0;i<C->N;i++)
  {
    for(j=0,e=t=0.0;j<C->Nn[i];j++)
    {
      p  = C->ATMN[i] + C->ATMN[C->Ni[i][j]];
      r  = NDX(C,i,j);
      u  = FC(C->WW[p][5],C->WW[p][6],r);

      //===== Gupta potential =====
      if( C->MODT==2 )
      {
        y  = r/C->WW[p][4] - 1.0;
        e += exp(-y*C->WW[p][1])*u*C->WW[p][0];      
        t += exp(-y*C->WW[p][3]*2.0)*u*C->WW[p][2]*C->WW[p][2];
      }
      //===== Sutton-Chen potential =====
      if( C->MODT==3 )
      {
        e += pow(C->WW[p][0]/r,C->WW[p][1])*u*0.5*C->WW[2*C->ATMN[i]][4];
        t += pow(C->WW[p][0]/r,C->WW[p][3])*u*pow(C->WW[2*C->ATMN[i]][2]*C->WW[2*C->ATMN[i]][4],2.0);
      }
      //===== Lennard-Jones potential =====
      if( C->MODT==4 )
      {
	r  = C->WW[p][1]/r;
        e += C->WW[p][0]*u*(-pow(r,C->WW[p][2])+pow(r,C->WW[p][3]));
      }
      //===========================
    }
    C->EA[i] = e - sqrt(t);
    E += C->EA[i];
  }
  }
  return E;
}
//==================================================        
double FRC_POT(Cell *C)
{
  int i,j,q,p,nth;
  double r,y,t,e,*h,u,z,E,**us,**hs;

  h  = make_d1D(C->N);
  hs = make_d2D(C->NP,C->N);
  us = make_d2D(C->NP,6);

  for(j=0;j<C->N;j++)
    h[j] = 0.0;
  
  for(i=0;i<C->NP;i++)
    for(j=0;j<C->N;j++)
      hs[i][j] = 0.0;

  LIST(C,0);
  E = 0.0;

#pragma omp parallel private(p,r,u,e,y,t,j) num_threads(C->NP)
  {
    nth = omp_get_thread_num();
#pragma omp for reduction (+:E) schedule(dynamic,C->NB)
  for(i=0;i<C->N;i++)
  {
    for(j=0,e=t=0.0;j<C->Nn[i];j++)
    {
      p = C->ATMN[i] + C->ATMN[C->Ni[i][j]];
      r  = NDX(C,i,j);
      u  = FC(C->WW[p][5],C->WW[p][6],r);

      //===== Gupta potential =====
      if( C->MODT==2 )
      {
        y  = r/C->WW[p][4] - 1.0;
        e += exp(-y*C->WW[p][1])*u*C->WW[p][0];      
        t += exp(-y*C->WW[p][3]*2.0)*u*C->WW[p][2]*C->WW[p][2];
      }
      //===== Sutton-Chen potential ===== 
      if( C->MODT==3 )
      {
        e += pow(C->WW[p][0]/r,C->WW[p][1])*u*0.5*C->WW[2*C->ATMN[i]][4];
        t += pow(C->WW[p][0]/r,C->WW[p][3])*u*pow(C->WW[2*C->ATMN[i]][2]*C->WW[2*C->ATMN[i]][4],2.0);
      }
      //===== Lennard-Jones potential =====
      if( C->MODT==4 )
      {
        r  = C->WW[p][1]/r;
        e += C->WW[p][0]*u*(-pow(r,C->WW[p][2])+pow(r,C->WW[p][3]));
      }
      //===========================
    }
    hs[nth][i] = sqrt(t);
    C->EA[i] = e - sqrt(t);
    E += C->EA[i];
  }
  }

  for(i=0;i<C->NP;i++)
    for(j=0;j<C->N;j++)
      h[j] += hs[i][j];
  
  for(i=0;i<C->NP;i++)
    for(q=0;q<6;q++)
      us[i][q] = 0.0;

  for(q=0;q<6;q++)
    C->U[q] = 0.0;

  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      C->F[i][q] = 0.0;

#pragma omp parallel private(nth,p,r,u,z,y,t,j,q) num_threads(C->NP)
  {
    nth = omp_get_thread_num();
#pragma omp for schedule(dynamic,C->NB)
  for(i=0;i<C->N;i++)
  {
    for(j=0,t=0.0;j<C->Nn[i];j++)
    {
      p = C->ATMN[i] + C->ATMN[C->Ni[i][j]];
      r = NDX(C,i,j);
      u =  FC(C->WW[p][5],C->WW[p][6],r);
      z = dFC(C->WW[p][5],C->WW[p][6],r);

      //===== Gupta potential =====
      if( C->MODT==2 )
      {
        y = r/C->WW[p][4] - 1.0;      
        t = (  2.0*  C->WW[p][0]*            (     u*C->WW[p][1]/C->WW[p][4]-z)*exp(-y*C->WW[p][1])
	       +0.5*( C->WW[p][2]*C->WW[p][2]*(-2.0*u*C->WW[p][3]/C->WW[p][4]+z)*exp(-y*C->WW[p][3]*2.0)*(1.0/h[i]+1.0/h[C->Ni[i][j]]) ))/r;
      }
      //===== Sutton-Chen potential =====
      if( C->MODT==3 )
      {
        t =  0.5*((C->WW[2*C->ATMN[i]][4]+C->WW[2*C->ATMN[C->Ni[i][j]]][4]) * (C->WW[p][1]/r*u-z) * pow(C->WW[p][0]/r,C->WW[p][1]) -
	          (C->WW[p][3]/r*u-z)*pow(C->WW[p][0]/r,C->WW[p][3])*
		  (pow(C->WW[2*C->ATMN[i]]          [2]*C->WW[2*C->ATMN[i]]          [4],2.0)/h[i] +
		   pow(C->WW[2*C->ATMN[C->Ni[i][j]]][2]*C->WW[2*C->ATMN[C->Ni[i][j]]][4],2.0)/h[C->Ni[i][j]])) /r;
      }
      //===== Lennard-Jones potential =====
      if( C->MODT==4 )
      {
	t = C->WW[p][1]/r;
        t = 2.0*C->WW[p][0]*( -pow(t,C->WW[p][2])*(z+u*C->WW[p][2]/r) + pow(t,C->WW[p][3])*(z+u*C->WW[p][3]/r) )/r;
      }
      //===========================

      for(q=0;q<3;q++)
	C->F[i][q] -= t*DX(C,i,j,q);
      for(q=0;q<3;q++)
	us[nth][q]   += 0.5*t*DX(C,i,j,q)*DX(C,i,j,q);
      for(q=0;q<3;q++)
	us[nth][q+3] += 0.5*t*DX(C,i,j,q)*DX(C,i,j,(q+1)%3);
    }
  }
  }

  for(i=0;i<C->NP;i++)
    for(q=0;q<6;q++)
      C->U[q] += us[i][q];

  free_d2D(hs,C->NP);
  free_d2D(us,C->NP);
  free_d1D(h);
  return E;
}
//===================================================
// read a classical model from 'model'
//===================================================
void READ_POT(Cell *C, char *dir)
{
  int  k,n,m,i;
  FILE *in;
  char s[200],f[200];

  if( C->MODT==0 )
    return;

  sprintf(f,"%s/model",dir);
  if( (in=fopen(f,"r")) == 0 )
  {
    fprintf(stderr,"ERROR opening %s\n",f);
    exit(1);
  }

  fgets(s,200,in);
  fgets(s,200,in);
  sscanf(s+2,"%s",s);
  C->MODT = 0;
  if( strncmp(s,"neural",6) == 0 )
    C->MODT = 1;
  if( strncmp(s,"Gupta",5) == 0 )
    C->MODT = 2;
  if( strncmp(s,"Sutton-Chen",11) == 0 )
    C->MODT = 3;
  if( strncmp(s,"Lennard-Jones",13) == 0 )
    C->MODT = 4;
  if( C->MODT<2 )
  {
    fprintf(stderr,"ERROR in %s\n",s);
    exit(1);
  }
  C->NW = 7;
  
  while( fgets(s,200,in) )
    if( strncmp(s,"|  number of species    |",25) == 0 )
    {
      sscanf(s+26,"%d" ,&C->nspc);
      break;
    }
  while( fgets(s,200,in) )
    if( strncmp(s,"|  species types        |",25) == 0 )
    {
      for(i=0,n=0,k=26; i < C->nspc;i++,k+=n)
	sscanf(s+k,"%d%n",&C->spcz[i],&n);
      break;
    }

  while( fgets(s,200,in) )
    if( strncmp(s,"|                               parameters",42) == 0 )
      break;
  fgets(s,200,in);

  //===== insert check of species =====

  for(n=0;n<C->nspc;n++)
    for(m=0;m<C->nspc;m++)
      if( !(n==0&&m==1) )
      for(k=0;k<C->NW;k++)
      {
	fgets(s,200,in);
	sscanf(s+26,"%lf",&C->WW[n*n+m*m][k]);
      }

  fclose(in);
  
  for(n=0, C->Rc = C->rc = 0.0; n<C->nspc; n++)
    if( C->Rc < C->WW[2*n*n][C->NW-1] )
      C->Rc = C->rc = C->WW[2*n*n][C->NW-1];
  
  if(0)
    for(n=0;n<C->nspc;n++)
      for(m=n;m<C->nspc;m++)
      {
	printf("%3d %3d   ",C->spcz[n],C->spcz[m]);
	for(k=0;k<C->NW;k++)
	  printf("% 12.6lf",C->WW[n*n+m*m][k]);
	printf("\n");
      }
}
//===================================================
