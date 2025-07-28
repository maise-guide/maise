#include "eutl.h"

//==================================================================
//     Check if atoms are too close; MUCH faster than using List
//==================================================================
int CHCK_Rm(Cell *C, double *Rm, double s)
{
  int i,j,q,q1,ic[3],N[3];
  double x,r;

  for(q=0;q<3;q++)
    N[q] = 1;
  if(C->ND==2)
    N[2] = 0;

  for(i=0;i<C->N;i++)
  {
    for(j=0;j<C->N;j++)
      if(j!=i)
      {
	for(q=0,r=0.0;q<D3;q++)
	  r += (C->X[j][q]-C->X[i][q])*(C->X[j][q]-C->X[i][q]);
	if( sqrt(r) < (Rm[C->ATMN[i]]+Rm[C->ATMN[j]])*s )
	  return 0;
      }
    if(C->ND>0)
    for(j=0;j<C->N;j++)
      for(ic[0]=-N[0];ic[0]<=N[0];ic[0]++)
      for(ic[1]=-N[1];ic[1]<=N[1];ic[1]++)
      for(ic[2]=-N[2];ic[2]<=N[2];ic[2]++)
      if( !(ic[0]==0&&ic[1]==0&&ic[2]==0) )
      {
	for(q=0,r=0.0;q<D3;q++)
        {
	  for(q1=0,x=0.0;q1<D3;q1++)
	    x += (double)(ic[q1]) *C->L[q1][q];
	  x += C->X[j][q]-C->X[i][q];
	  r += x*x;
	}
	if( sqrt(r) < (Rm[C->ATMN[i]]+Rm[C->ATMN[j]])*s )
	{
	  return 0;
	}
      }	
  }
  return 1;
}
//==================================================================
//     Pull atoms apart to help satisty Rm condition
//==================================================================
int ADJT_NP(Cell *C, double *Rm, double RM, double *b, int N)
{
  int i,j,q,n;
  double a,r,R,c[3];

  a = 0.2;
  for(n=0;n<N;n++)
  {
    if(0&&CHCK_Rm(C,Rm,0.5)==0)    // 0.5 here is NOT hard-core factor
      return 0;
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->F[i][q] = 0.0;
    for(i=0;i<C->N;i++)
      for(j=0;j<C->N;j++)
      if(i!=j)
      {
	R = Rm[C->ATMN[i]]+Rm[C->ATMN[C->Ni[i][j]]];
	r = NDR(C,i,j)/R;
	r = pow(r,-13.0);
	for(q=0;q<3;q++)
	  C->F[i][q] -= (C->X[j][q]-C->X[i][q])/R*r;
	for(q=0;q<3;q++)
	  c[q] = C->X[i][q]*b[q];
	r = VectorLen(c,3);
	//===== use attractive spring beyond targeted RM sphere =====
	if( r > RM )
	  for(q=0;q<3;q++)
	    C->F[i][q] -= C->X[i][q]/r*(r-RM)/R;
      }
    for(i=0;i<C->N;i++)
      if( (R=VectorLen(C->F[i],3)) > 2.0 )
	for(q=0;q<3;q++)
	  C->F[i][q] *= 2.0/R;
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->X[i][q] += a*C->F[i][q];
  }
  return CHCK_Rm(C,Rm,1.0);
}
//==================================================================
//     Pull atoms apart to help satisty Rm condition
//==================================================================
int ADJT_CL(Cell *C, double *Rm, int N)
{
  int i,j,q,n;
  double a,r,R;

  a = 0.2;
  for(n=0;n<N;n++)
  {
    if(CHCK_Rm(C,Rm,0.7)==0)    // 0.7 here is NOT hard-core factor
      return 0;
    LIST(C,0);
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->F[i][q] = 0.0;
    for(i=0;i<C->N;i++)
      for(j=0;j<C->Nn[i];j++)
      {
	R = Rm[C->ATMN[i]]+Rm[C->ATMN[C->Ni[i][j]]];
	r = NDX(C,i,j)/R;
	r = pow(r,-13.0);
	for(q=0;q<3;q++)
	  C->F[i][q] -= DX(C,i,j,q)/R*r;
      }
    for(i=0;i<C->N;i++)
      if( (R=VectorLen(C->F[i],3)) > 2.0 )
	for(q=0;q<3;q++)
	  C->F[i][q] *= 2.0/R;
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->X[i][q] += a*C->F[i][q];
    LIST(C,0);
  }
  JAR(C);
  return CHCK_Rm(C,Rm,1.0);
}
//==================================================================
//  Mate two lattices
//==================================================================
void MATE_LT(Cell *C1, Cell *C2, Cell *C, double s)
{
  int i,k,q;
  double w,W,t1,t2;

  for(i=k=0,W=100000.0;i<3;i++)
  {
    for(q=0,w=0.0;q<3;q++)
    {
      t1  = ( C1->LAT[q]-C2->LAT[(q+i)%3] )/( C1->LAT[q]+C2->LAT[(q+i)%3] );
      t2  = ( C1->ANG[q]-C2->ANG[(q+i)%3] )/( C1->ANG[q]+C2->ANG[(q+i)%3] );
      w  += t1*t1+t2*t2;
    }
    if(w<W)
    {
      k = i;
      W = w;
    }
  }
  if(Random()<s)
    s = 0.0;
  else
    s = 1.0;

  if(C->ND==2)
    k = 0;

  for(q=0;q<3;q++)
  {
    C->LAT[q] = s*C1->LAT[q] + (1.0-s)*C2->LAT[(q+k)%3];
    C->ANG[q] = s*C1->ANG[q] + (1.0-s)*C2->ANG[(q+k)%3];
  }
  ABC_LT(C);
}
//==================================================================
//
//==================================================================
int SHKE_CL(Cell *C, double dL, double dX)
{
  int i,q,k;
  double e[3][3],L[3][3];

  Relative(C);
  e[0][0] = 1.0+dL*RANG();   //create lattice mutation matrix
  e[1][1] = 1.0+dL*RANG();
  e[2][2] = 1.0+dL*RANG();
  e[0][1] = e[1][0] = 0.5*dL*RANG();
  e[0][2] = e[2][0] = 0.5*dL*RANG();
  e[1][2] = e[2][1] = 0.5*dL*RANG();

  if(C->ND==2)
  {
    e[0][2] = e[1][2] = e[2][1] = e[2][0] = 0.0;
    e[2][2] = 1.0;
  }

  for(i=0;i<3;i++)           //creating a local version of the lattice matrix
    for(q=0;q<3;q++)
      L[i][q] = C->L[i][q];

  for(i=0;i<3;i++)           //applying mutation matrix to the lattice matrix
    for(q=0;q<3;q++)
      for(k=0,C->L[i][q]=0.0;k<3;k++) 
        C->L[i][q] += e[i][k]*L[k][q];
  abc(C);
  Real(C);
  if( C->LAT[0]<C->R0 || C->LAT[1]<C->R0 || C->LAT[2]<C->R0 )
    return 0;

  if(C->ND!=2)                      //perturb atomic positions
    for(i=0;i<C->N;i++)             //checking 
      for(q=0;q<3;q++)
        C->X[i][q] += dX*RANG();
  if(C->ND==2)
    for(i=0;i<C->N;i++)
      for(q=0;q<C->ND;q++)
        C->X[i][q] += dX*RANG();

  return 1;
}
//==================================================================
//     
//==================================================================
void RAND_VC(double *a)
{
  int q;
  double  r;

  r = 2.0;
  while( r>1.0 )
  {
    for(q=0;q<3;q++)
      a[q] = 2.0*(0.5-Random());
    r = VectorLen(a,3);
  }
}
//==================================================================
//     
//==================================================================
void RAND_LV(Cell *C)
{
  int     i,q;
  double  V,L[3];

  V = 0.0;
  while( V<0.4 )
  {
    for(i=0;i<3;i++)
    {
      L[i] = 2.0;
      while( L[i]>1.0 )
      {
	for(q=0;q<3;q++)
	  C->L[i][q] = 2.0*(0.5-Random());
	L[i] = VectorLen(C->L[i],3);
      }
      VectorNorm(C->L[i]);
    }
    V = fabs(CELL_VOL(C));
  }
  
  for(i=0;i<3;i++)
  {
    L[i] = pow( 1.2 + 0.55*(0.5-Random())*2.0,1.0/3.0 );
    for(q=0;q<3;q++)
      C->L[i][q] *= L[i];
  }
  SCALE_LATTICE(C,pow(CELL_VOL(C),-1.0/3.0));
  Lat_Order(C);
  Lat_Align(C);
  if(C->L[2][2]<0.0)
    C->L[2][2] = -C->L[2][2];
}
//==================================================================
//     Generate random Cell
//==================================================================
void RAND_CL(Tribe *T, Cell *C, Cell *D, int J)
{
  int i,q,s;
  double a;
  char buf[200];

  Copy_C(C,D);     //COPY STRUCTURE READ IN FROM FILE TO A SECOND C - to work from

  a = 1.0;
  for(s=0;s<T->Nm;s++)
  {
    if(J==0)       // J=2 keep LV constant start from given X;  J=1 keep LV constant; J=0 generate random LV
    {
      RAND_LV(C);
      if(T->ND==3)
	a = pow( 2.0*T->VOL*(1.0+0.25*(0.5-Random())*2.0),1.0/3.0);
      if(T->ND==2)
        a = pow( 2.0*T->VOL*(1.0+0.25*(0.5-Random())*2.0),1.0/2.0);

      SCALE_LATTICE(C,a);

      if(T->ND==2)
      { 
	C->L[2][0] = C->L[2][1] = C->L[1][2] = C->L[0][2] = 0.0;
	C->L[2][2] = T->B;
	for(i=0;i<C->N;i++)
          if(C->FF[i][2]!=0)
	    C->X[i][2] = 0.5;
      }
      else
	SHRT_LV(C);

      for(i=0;i<C->N;i++)
        for(q=0;q<T->ND;q++)
          if (C->FF[i][q]!=0)
            C->X[i][q] = Random();  //SET RANDOM COORDS
      if(T->ND==2)
        for(i=0;i<C->N;i++)
          if (C->FF[i][2]!=0)
            C->X[i][2] = 0.5;

      C->XT = 0;
      Real(C);

      if(CHCK_Rm(C,T->Rm,1.0)==1||ADJT_CL(C,T->Rm,10)==1)
	break;
    }
    if(J==1)                        //KEEP LV and random X (except for the fixed positions)
    {
      Relative(C);
      for(i=0;i<C->N;i++)
        for(q=0;q<T->ND;q++)
          if (C->FF[i][q]!=0)
	    C->X[i][q] = Random();
      if(T->ND==2)
        for(i=0;i<C->N;i++)
          if (C->FF[i][2]!=0)
            C->X[i][2] = 0.5;
      C->XT = 0;
      Real(C);
      if(CHCK_Rm(C,T->Rm,1.0)==1 || ADJT_CL(C,T->Rm,10)==1)
	break;
    }
    if(J==2 )                       //from given X
    {
      for(i=0;i<C->N;i++)
	for(q=0;q<D3;q++)
          if (C->FF[i][q]!=0)
	    C->X[i][q] = D->X[i][q]+T->ca*C->R0*RANG();;
      if(T->ND==2)
	for(i=0;i<C->N;i++)
          if (C->FF[i][2]!=0)
	    C->X[i][2] = D->X[i][2];

      if(CHCK_Rm(C,T->Rm,1.0)==1)
	break;
    }
  }
  if(s==T->Nm)
  {
    sprintf(buf,"exit in INIT_CL: exceeded %3d\n",T->Nm);
    fprintf(stderr,"exit in INIT_CL: exceeded %3d\n",T->Nm);
    Print_LOG(buf);
    exit(1);
  }
}
//==================================================================
//     Rank Tribe: Sorting is N^2 but N is small
//==================================================================
double COMP_CL(Cell *C, Cell *D)
{
 
  double V,V1,V2,DV,DE,AV,AE,AR,W;

  DE = 0.050;
  DV = 0.010;

  AE = 0.1;
  AV = 0.1;
  AR = 0.8;

  V1 = CELL_VOL(C);
  V2 = CELL_VOL(D);
  V  = 0.5*(V1-V2)/(V1+V2);

  W = AR*CxC(C,D) + AE*exp( -0.5*(C->P-D->P)*(C->P-D->P)/(DE*DE) ) + AV*exp( -0.5*V*V/(DV*DV) );

  return W;
}
//==================================================================
//     Randomize initial structures
//==================================================================
void TEMP_CL(Tribe *T, Cell *C, int p)
{
  int i,j,m,q,s;
  double V;
  char buf[200];

  for(m=0,s=0;m<T->Nm;m++)
  {
    Copy_C(C,&T->C[p]);
    V = CELL_VOL(&T->C[p]);
    if( ! (SHKE_CL(&T->C[p],T->cl,T->C[p].R0*T->ca)==0 || SHRT_LV(&T->C[p])==0 || fabs( CELL_VOL(&T->C[p])/V-1.0 ) > 0.5) )
    {
      if( T->JS!=2  && !(T->ml>(10^-10) || T->cl>(10^-10)) )              // added in 4.1
	SCALE_Cell(&T->C[p], pow(CELL_VOL(&T->C[p])/V,-1.0/3.0) + 0.1*(0.5-Random()));
      if(T->NSPC>1)
	for(i=s=0;i<T->C[p].N;i++)
	  if(Random()<T->cs)
	  {
	    for(j=i;T->C[p].ATMN[i]==T->C[p].ATMN[j%T->C[p].N];j++);
	    j = ( j + (int)(Random()*(double)(T->C[p].N-T->SPCN[T->C[p].ATMN[i]])) )%T->C[p].N;
	    for(q=0;q<3;q++)
            {
	      dSwap(&T->C[p].X[i][q],&T->C[p].X[j][q]);
	      iSwap(&T->C[p].FF[i][q],&T->C[p].FF[j][q]);
            }
	    s++;
	  }
      if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1)
      {
	T->P1[p] = T->P2[p] = T->C[p].P;
	sprintf(buf,"gen %4d: %2d %2d are proud parents of cell %2d with                                                 %3d swaps\n",T->n,p ,p ,p,               s);
	Print_LOG(buf);
	break;
      }
    }
  }
  
  if(m==T->Nm)
  {
    sprintf(buf,"failed to clone in %d tries\n",T->Nm);
    fprintf(stderr,"failed to clone in %d tries\n",T->Nm);
    Print_LOG(buf);
    exit(1);
  }
}
//==================================================================
