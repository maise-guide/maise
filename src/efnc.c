#include "efnc.h"

#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"

//==================================================================
// Read in block information for tetris
//==================================================================
int READ_BLOCK(Cell *C, int II[30], int JJ[30], int S[30][30], double RR[30][30], double AA[30][30], double ZZ[30][30], double *Rc, int *MIN, int *MAX)
{
  int  n,k,Q;
  FILE *in;
  char buf[200];

  if( in = fopen("INI/BLOCK","r") )
  {
    //===== read in blocks for TETRIS =====
    fgets(buf,200,in);
    fgets(buf,200,in);
    sscanf(buf,"%d",&Q);
    for(k=0;k<Q;k++)
    {
      fgets(buf,200,in);
      JJ[k] = 0;
      sscanf(buf,"%d %d",&II[k],&JJ[k]);
      for(n=0;n<II[k];n++)
      {
	fgets(buf,200,in);
	sscanf(buf,"%d %lf %lf %lf",&S[k][n],&RR[k][n],&AA[k][n],&ZZ[k][n]);
      }
    }
    //===== read in bonds for TETRIS and RANDOM =====
    fgets(buf,200,in);
    for(k=0;k<C->NSPC;k++)
    {
      fgets(buf,200,in);
      sscanf(buf,"%lf %d %d",&Rc[k],&MIN[k],&MAX[k]);
      Rc[k] *= Rc[k];
    }
    fclose(in);
  }
  //===== if INI/BLOCK is not specified just use single atoms of each species =====
  else
  {
    Q = C->NSPC;
    for(k=0;k<Q;k++)
    {
      II[k]    = 1;
      JJ[k]    = 0;
      S[k][0]  = k;
      RR[k][0] = 0.0;
      AA[k][0] = 0.0;
      ZZ[k][0] = 0.0;
      Rc[k]    = 0.0;
      MIN[k]   = 0;
      MAX[k]   = 0;
    }
  }

  return Q;
}
//==================================================================
void ROT_BLOCK(int N, double *RR, double *AA, double *ZZ, int *S)
{
  int    i,q,k;
  double R[3][3],X[30][3],r[3];

  //===== convert polar coordinates to cartesian =====
  for(i=0;i<N;i++)
  {
    X[i][0] = RR[i]*cos(AA[i]*Pi/180.0);
    X[i][1] = RR[i]*sin(AA[i]*Pi/180.0);
    X[i][2] = ZZ[i];
  }
  //===== shift the center of mass to the origin =====
  for(q=0;q<3;q++)
    for(i=0,r[q]=0.0;i<N;i++)
      r[q] += X[i][q]/(double)N;
  for(i=0;i<N;i++)
    for(q=0;q<3;q++)
      X[i][q] -= r[q];
  //===== generate a random rotation matrix with quaternions =====
  RAND_ROT(R);

  //===== rotate all atoms in the block by random angles =====
  for(i=0;i<N;i++)
  {
    for(q=0;q<3;q++)
      for(k=0,r[q]=0.0;k<3;k++)
	r[q] += R[q][k]*X[i][k];
    for(q=0;q<3;q++)
      X[i][q] = r[q];
  }

  //===== check the center of mass after rotation =====
  for(q=0;q<3;q++)
    for(i=0,r[q]=0.0;i<N;i++)
      r[q] += X[i][q]/(double)N;
  if( fabs(r[0])>1E-12 || fabs(r[1])>1E-12 || fabs(r[2])>1E-12 )
    printf("WARNING in ROTATION:% lf % lf % lf\n",r[0],r[1],r[2]);

  //===== convert cartesian coordinates back to polar =====
  for(i=0;i<N;i++)
  {
    RR[i] = sqrt(X[i][0]*X[i][0]+X[i][1]*X[i][1]);
    if( RR[i]<1e-12 )
      AA[i] = 0.0;
    else
      AA[i] = acos( X[i][0]/RR[i] )*180.0/Pi;
    if( X[i][1] < 0.0 )
      AA[i] *= -1.0;
    ZZ[i] = X[i][2];
  }
  //===== place the first atom to have the lowest z =====
  for(k=0,i=1,r[2]=ZZ[0];i<N;i++)
    if( ZZ[i]<r[2] )
    {
      r[2] = ZZ[i];
      k    = i;
    }
  dSwap(&RR[0],&RR[k]);
  dSwap(&AA[0],&AA[k]);
  dSwap(&ZZ[0],&ZZ[k]);
  iSwap( &S[0], &S[k]);

  for(i=0;i<N;i++)
    ZZ[i] -= r[2];

}
//==================================================================
int CHCK_BOND(Cell *C, double *Rc, int *MIN, int *MAX)
{
  int    i,j,q,q1,ic[3];
  double x,r;

  for(i=0;i<C->N;i++)
  if( MAX[C->ATMN[i]]>0 )  // if MAX for a species is zero, the condition is fullfilled automatically
  {
    C->Nn[i] = 0;
    for(j=0;j<C->N;j++)
      for(ic[0]=-1;ic[0]<=1;ic[0]++)
      for(ic[1]=-1;ic[1]<=1;ic[1]++)
      for(ic[2]=-1;ic[2]<=1;ic[2]++)
      if(! (i==j&&ic[0]==0&&ic[1]==0&&ic[2]==0) )
      {
	for(q=0,r=0.0;q<D3;q++)
	{
	  for(q1=0,x=0.0;q1<D3;q1++) 
	    x += (double)(ic[q1]) *C->L[q1][q]; 
	  x += C->X[j][q]-C->X[i][q]; 
	  r += x*x; 
	}
	if(r<Rc[C->ATMN[i]])
	  C->Nn[i]++;
      }
  }
    
  for(i=0;i<C->N;i++)
    if( MAX[C->ATMN[i]]>0 && ( C->Nn[i]<MIN[C->ATMN[i]] || C->Nn[i]>MAX[C->ATMN[i]] ) )
      return 0;
  return 1;
}
//==================================================================
// Generate starting order of atoms for tetris
//==================================================================
int GEN_TETR(Tribe *T, Cell *C, double *D, int *I, int *B, int *U, int Y[30], int S[30][30], int Q)
{
  int    i,j,m,k,n,N[10],O[10],M,K,q,p,F[10];
  double fi,fj,e,t;

  Random();

  for(p=0;p<10000;p++)
  {
    for(k=0;k<C->NSPC;k++)
      N[k] = 0;
    for(k=1,O[0]=0;k<T->NSPC;k++)
      O[k] = O[k-1]+T->SPCN[k-1];
    for(n=0;n<Q;n++)                 // any block is possible
      F[n] = 1;
    
    for(i=K=0,q=Q;i<C->N&&q>0;)
    {
      m = (int)(Random()*(double)q); // pick random block

      for(n=j=0;n<Q&&j<=m;n++,j++)
	if(F[n]==0)
	  j--;
      n--;
      for(k=0;k<Y[n];k++)
	N[S[n][k]]++;
      for(k=0;k<C->NSPC;k++)         // check if block addition makes any species exceed limit
	if( N[k]>C->SPCN[k] )
	  break;
      if( k<C->NSPC )                // if some species exceeds limit, exclude nth block
      {
	for(k=0;k<Y[n];k++)
	  N[S[n][k]]--;
	F[n] = 0;
	q--;
      }
      else
      {
	U[K] = i;                    // atom number after K blocks
	for(k=0;k<Y[n];k++,i++)      // species type for ith atom
	  I[i] = O[S[n][k]]++;       // type of Kth block
	B[K++] = n;
      }
      // printf("%4d   Qq %3d %3d   mn %3d %3d   N %3d   BL %3d %3d  SP %3d %3d\n",p,Q,q,m,n,i,F[0],F[1],N[0],N[1]);
    }
    if(q>0)
      break;
  }

  if(p==1000)
  {
    printf("Could not select any block combinations to match the number of each species\n");
    exit(0);
  }

  return K;
}
//==================================================================
// Generate starting order of atoms for random or core-shell NPs
//==================================================================
void GEN_TYPE(Tribe *T, Cell *C, double *D, int *I)
{
  int i,k,n,N[10],K[10];

  printf("%d\n",T->NORD);

  if(T->NSPC==1||T->NORD==1)
  {
    for(i=0;i<C->N;i++)
      I[i] = i;
    return;
  }
  if(T->NORD==0)
  {
    for(i=0;i<C->N;i++)
      D[i] = Random();
    Sort(D,I,C->N);      // pick random i not to bias atom type           
    return;
  }
  if(T->NORD==2)
  {
    //===== define random species number =====    
    for(k=0;k<C->NSPC;k++)
      D[k] = Random();
    Sort(D,K,C->NSPC);
    for(k=1,N[0]=0;k<T->NSPC;k++)
      N[k] = N[k-1]+T->SPCN[k-1];

    for(k=i=0;k<C->NSPC;k++)
      for(n=0;n<C->SPCN[K[k]];n++,i++)
	I[i] = N[K[k]] + n;
    for(i=0;i<C->N;i++)
      printf("%3d %3d\n",i,I[i]);
  }
  return;
  //===== meant for PACK but not used anymore =====
  if(T->NORD>0)
  {
    n = (int)(Random()*(double)T->NORD);  // pick random set from INI/ORDINI
    for(k=1,N[0]=0;k<T->NSPC;k++)
      N[k] = N[k-1]+T->SPCN[k-1];         // set counter for each species
    for(i=0;i<C->N;i++)
    {
      D[i] = T->CORD[n][i];               // set i to min distance along c
      I[i] = N[T->TORD[n][i]]++;          // set i to species and increase counter
    }
    return;
  }
}
//=========================================================================
// Gives distance between atoms i and j; with x-y boundary conditions
//=========================================================================
double NDB(Cell *C, int i, int j)
{
  int    m,n,q;
  double x,r,t,d;

  t = 0.0;
  x = C->X[j][2] - C->X[i][2];
  r = x*x;
  d = 1000.0;
  for(m=-1;m<=1;m++)
    for(n=-1;n<=1;n++)
    {
      for(q=0,t=0.0;q<2;q++)
      {
	x = C->X[j][q] - C->X[i][q] + (double)m*C->L[0][q] + (double)n*C->L[1][q];
	t += x*x;
      }
      if( t<d )
	d = t;
    }
  r += d;
  return sqrt(r);
}
//==================================================================
// Order atoms by type
//==================================================================
void ORDER_TYPE(Tribe *T, int J)
{
  int  p,i,j,k,q,qq;
  char buf[200];
  
  for(p=T->N;p<2*T->N;p++)
    if( p>=T->SES[J] && p<T->FES[J] )
    for(qq=i=0;qq<T->NSPC-1;qq++)
      for(j=0;j<T->C[p].SPCN[qq];j++)
      {
	for(k=i;k<T->C[p].N;k++)
	  if(T->C[p].ATMN[k]==qq)
	    break;
	if(k==T->C[p].N)
	{
	  sprintf(buf,"ERROR: atom types are wrong\n");
	  fprintf(stderr,"ERROR: atom types are wrong\n");
	  Print_LOG(buf);
	  exit(1);
	}
	for(q=0;q<D3;q++)
	  dSwap(&T->C[p].X[i][q],&T->C[p].X[k][q]);
	iSwap(&T->C[p].ATMN[i],&T->C[p].ATMN[k]);
	iSwap(&T->C[p].ATMZ[i],&T->C[p].ATMZ[k]);
	i++;
      }
}
//==================================================================
// Reshaping NP shells
//==================================================================
void NANO_CHOP(Tribe *T, int J, Cell *C, int P)
{
  int i,j,q,s,w,W,n,m,M,*I,N0,ok,k,p;
  double a[3],R,*D,L,dr,r[3],b[3],S,Q,t,**X;
  char buf[200];

  I = make_i1D(C->N);
  D = make_d1D(C->N);
  X = make_d2D(C->N,3);

  k = N0 = 0;
  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    if(P<0)
    {
      k = (int)(Random()*(double)T->N);
      Copy_C(&T->C[k],C);
    }

    for(i=0,R=0.0;i<C->N;i++)
      if( VectorLen(C->X[i],3) > R )
        R = VectorLen(C->X[i],3);

    dr = 0.02;
    W  = (int)pow(C->N,1.5);
    M  = 100;
    
    for(i=0;i<C->N;i++)
      I[i] = i;

    //===== Chop off a few atoms =====
    S = R;
    Q = sqrt((double)C->N);
    for(s=0,t=1.0;s<T->Nm;s++)
    {
      RAND_VC(b);
      VectorNorm(b);
      S = R*(0.65 + 0.35*Random())*t;
      for(i=0,N0=0;i<C->N;i++)
      if( DotProd(C->X[i],b,3) > S )
      {
        I[N0++] = i;
        for(q=0;q<3;q++)
          X[i][q] = C->X[i][q];
      }
      if( (int)(0.5*Q)<N0 && N0<(int)(1.5*Q) )
        break;
      if(s%(T->Nm/10)==0)
        t *= 0.9;
    }
    if(s==T->Nm)
    {
      sprintf(buf,"exit in CHOP_NP 1: exceeded %3d\n",T->Nm);
      fprintf(stderr,"exit in CHOP_NP 1: exceeded %3d\n",T->Nm);
      Print_LOG(buf);
      exit(1);
    }

    for(s=0;s<T->Nm;s++)
    {
      for(n=0;n<N0;n++)
        for(q=0,i=I[n];q<3;q++)
          C->X[i][q] = X[i][q];

      for(n=0,ok=1;n<N0;n++)
      {
        i = I[n];
        L = 10.0*R;

	for(q=0;q<3;q++)
	  r[q] = C->X[i][q];
        for(w=0;w<W;w++)
        {
	  for(q=0;q<3;q++)
	    a[q] = b[q];
	  //=====  pick a direction in the opposite hemisphere  =====
	  while( DotProd(a,b,3)>0.0 )	    
	    RAND_VC(a);
          VectorNorm(a);
          for(q=0;q<3;q++)
            C->X[i][q] = a[q]*(R+4.0*T->Rm[C->ATMN[i]])*(1.0+dr*Random());

	  for(m=0;m<M;m++)
	  {
	    for(q=0;q<3;q++)
	      C->X[i][q] *= (1.0-dr);
	    //=====  OK to go over all C->N  =====
	    for(j=0;j<C->N;j++)
              if(i!=j)
              {
                t = NDR(C,i,j);
                if( t < 1.5*(T->Rm[C->ATMN[i]]+T->Rm[C->ATMN[j]]) )
                {
                  j = C->N;
                  m = M;
                  if( t < 1.0*(T->Rm[C->ATMN[i]]+T->Rm[C->ATMN[j]]) )
                    ok = 0;
                }
              }
	  }
	  if( VectorLen(C->X[i],3) < L )
	  {
	    L = VectorLen(C->X[i],3);
	    for(q=0;q<3;q++)
	      r[q] = C->X[i][q];
	  }             

        }
        for(q=0;q<3;q++)
          C->X[i][q] = r[q];
      }
      if( ok==1 )
        break;
    }

    if(s==T->Nm)
    {
      sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      Print_LOG(buf);
      exit(1);
    }

    LIST(C,0);
    RDF(C,1);
    LIST(&T->C[k],0);
    RDF(&T->C[k],1);
    if( 1.0-CxC(C,&T->C[k])<1e-6 )
    {
      printf("%3d %3d\n",p,k);
      exit(0);      
    }
    
    Copy_C(C,&T->C[p]);
    if(P<0)
    {
      T->P1[p] = T->P2[p] = T->C[k].P;
      sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k,k);
      Print_LOG(buf);
    }
  }
  free_d2D(X,C->N);
  free_i1D(I);
  free_d1D(D);
}
//==================================================================
//     Generate NPs using randomized initial positions
//==================================================================
void NANO_RAND(Tribe *T, int J, Cell *C, int P)
{
  int i,q,s,n,*I,p,m;
  double V,R,*D,b[3],c[3],t;
  char buf[200];

  I = make_i1D(C->N);
  D = make_d1D(C->N);

  for(i=0,V=0.0;i<T->NSPC;i++)    
    V += pow(T->Rm[i],3.0)*(double)T->SPCN[i];
  R = pow(V,1.0/3.0)*1.5/T->COMP;

  for(i=0;i<3;i++)
    for(q=0;q<3;q++)
      C->L[i][q] = 0.0;
  for(i=0;i<3;i++)
    C->L[i][i] = T->B;

  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    printf("%s %3d\n",T->NES[J],p);
    //===== squash or stretch the sphere in x, y, or z =====
    for(q=0;q<3;q++)
      b[q] = 1.0;
    q = (int)(Random()*3.0);
    t = 1.0 + 2.0*T->te*(0.5-Random());
    b[(q+0)%3] = t;
    b[(q+1)%3] = b[(q+2)%3] = 1.0/sqrt(t);
    
    for(s=0;s<T->Nm;s++)
    {
      for(n=0;n<C->N;n++)
      {
	for(q=0;q<3;q++)
	  c[q] = 1.0;
	//===== generate an atom within radius 1.0 =====
	while( VectorLen(c,3)>1.0 )
	  for(q=0;q<3;q++)
	  {
	    T->C[p].X[n][q] = 1.0-2.0*Random();
	    c[q] = T->C[p].X[n][q]*b[q];
	  }
      }
      //===== for generation of core-shell structures =====
      if(0) 
      {
	for(i=0;i<T->C[p].N;i++)
	  D[i] = VectorLen(T->C[p].X[i],3);
	Sort(D,I,T->C[p].N);      
      }
      GEN_TYPE(T,C,D,I);
      for(i=0;i<T->C[p].N;i++)
	for(q=0;q<3;q++)
	  C->X[i][q] = R*T->C[p].X[I[i]][q];
      //===== adjust interatomic distances =====
      for(m=0;m<C->N/2;m++)
	if(ADJT_NP(C,T->Rm,R,b,1))
	  break;
      if( m<C->N/2 )
	break;
    }

    if(s==T->Nm)
    {
      sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      Print_LOG(buf);
      exit(1);
    }
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->X[i][q] += 0.5*C->L[q][q];
    Copy_C(C,&T->C[p]);
    sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],-1,-1);
    if(P>=0)
      Print_LOG(buf);
  }
  ORDER_TYPE(T,J);

  free_i1D(I);
  free_d1D(D);

}
//==================================================================
//     Generate random NPs using TETRIS algorithm
//==================================================================
void NANO_TETR(Tribe *T, int J, Cell *C, int P)
{
  int i,j,q,s,w,W,n,m,M,*I,p;
  double a[3],A,V,R,*D,L,dr,r[3],b[3],c[3];
  char buf[200];

  I = make_i1D(C->N);
  D = make_d1D(C->N);

  for(i=0,V=0.0;i<T->NSPC;i++)    
    V += pow(T->Rm[i],3.0)*(double)T->SPCN[i];
  R = pow(V,1.0/3.0)*3.0;

  dr = 0.02;
  W  = (int)pow(C->N,1.5);
  M  = (int)(R/dr)-1;

  for(i=0;i<3;i++)
    for(q=0;q<3;q++)
      C->L[i][q] = 0.0;
  for(i=0;i<3;i++)
    C->L[i][i] = T->B;
  for(q=0;q<3;q++)
    r[q] = 0.0;

  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    printf("%s %3d\n",T->NES[J],p);
    //===== squash or stretch the sphere in x, y, or z =====
    for(q=0;q<3;q++)
      b[q] = 1.0;
    b[(int)(Random()*3.0)] = 1.0 + 2.0*T->te*(0.5-Random());

    for(s=0;s<T->Nm;s++)
    {
      GEN_TYPE(T,C,D,I);
      
      for(n=0;n<C->N;n++)
      {
        i = I[n];
        L = R*(1.0+dr);
        for(w=0,j=0;w<W;w++)
        {
          RAND_VC(a);
          A = VectorLen(a,3);
          for(q=0;q<3;q++)
            C->X[i][q] = a[q]/A*L*(1.0+dr*Random())*b[q];
          for(m=0;m<M;m++)
          {
            for(q=0;q<3;q++)
              C->X[i][q] *= (1.0-dr);
            for(j=n-1;j>=0;j--)
              if( NDR(C,i,I[j]) < 1.5*(T->Rm[C->ATMN[i]]+T->Rm[C->ATMN[I[j]]]) )
              {
                m = M;
		break;
              }
	  }

          for(q=0;q<3;q++)
            c[q] = C->X[i][q]*b[q];

          if( VectorLen(c,3) < L )
          {
            for(;j>=0;j--)
              if( NDR(C,i,I[j]) < 1.0*(T->Rm[C->ATMN[i]]+T->Rm[C->ATMN[I[j]]]) )
		break;
	    if(j==-1)
	    {
	      L = VectorLen(c,3);
	      for(q=0;q<3;q++)
		r[q] = C->X[i][q];
	    }
          }
        }
        for(q=0;q<3;q++)
          C->X[i][q] = r[q];

	//===== shake the core N^1/3 atoms to avoid bias at the origin =====
	if( T->JITT = 1 && n == (int)pow((double)C->N,1.0/3.0) )
	{	  
	  //===== vector a is a random point inside a sphere of radius 1.0 =====
	  //===== to be a uniform shift, rely only on the cutoff of atom i =====
          RAND_VC(a);
	  for( j=0;j<=n;j++ )
	    for(q=0;q<3;q++)
	      C->X[I[j]][q] += a[q]*b[q]*T->Rm[C->ATMN[i]];
	}
      }
      for(i=0;i<C->N;i++)
        for(q=0;q<3;q++)
          C->X[i][q] += 0.5*C->L[q][q];
      if( CHCK_Rm(C,T->Rm,1.0)==1 )
        break;
    }

    if(s==T->Nm)
    {
      sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      Print_LOG(buf);
      exit(1);
    }
    Copy_C(C,&T->C[p]);
    sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],-1,-1);
    if(P>=0)
      Print_LOG(buf);
  }
  ORDER_TYPE(T,J);

  free_i1D(I);
  free_d1D(D);

}
//==================================================================
// Plant a few NP seeds from INI/POSCARXXX or TETR
//==================================================================
void NANO_PLNT(Tribe *T, int J, int P)
{
  int i,j,q,s,w,W,m,M,*I,ok,k,A,N,N0,p;
  double a[3],V,R,R0,*D,L,dr,r[3],t,c[3],b[3],d[3];
  char buf[200];

  N  = 2*T->N;

  for(q=0;q<3;q++)
    r[q] = 0.0;

  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    k = (int)((double)T->pos*Random());
    //=====  always plant at least one best seed in the first generation  =====
    if( p==T->SES[J] && T->n==0 )
      k = 0;

    sprintf(buf,"INI/POSCAR%03d",k);
    if( READ_CELL(&T->C[N],buf)==0 )
    {
      NANO_TETR(T,0,&T->C[N],p);      
      k = -1;
    }
    else
    {
    printf("%s %3d\n",T->NES[J],p);
    for(i=0;i<3;i++)
      for(q=0;q<3;q++)
        T->C[p].L[i][q] = 0.0;
    for(i=0;i<3;i++)
      T->C[p].L[i][i] = T->B;
    
    A  = T->C[p].N;
    I  = make_i1D(T->C[N].N);
    D  = make_d1D(T->C[N].N);
    
    for(i=0,R0=0.0;i<T->NSPC;i++)
      R0 += T->Rm[i]*(double)T->SPCN[i]/(double)T->C[p].N;
    R0 *= 1.5*2.0;

    for(i=0,V=0.0;i<T->NSPC;i++)
      V += pow(T->Rm[i],3.0)*(double)T->SPCN[i];
    R = pow(V,1.0/3.0)*3.0;
    
    RAND_VC(d);
    VectorNorm(d);

    if(T->n>0)
    for(i=0;i<T->C[N].N;i++)
      for(q=0;q<3;q++)
        T->C[N].X[i][q] += d[q]*R0 + T->ca*RANG();

    for(i=0;i<T->C[N].N;i++)
      D[i] = VectorLen(T->C[N].X[i],3);
    Sort(D,I,T->C[N].N);

    //===== Currently only for 1 species =====

    for(i=0;i<T->C[N].N;i++)
      for(q=0;q<3;q++)
        T->C[p].X[i][q] = T->C[N].X[I[i]][q];

    dr = 0.01;
    W  = (int)pow(A,1.5);
    M  = (int)(R/dr)-1;
    
    N0 = (int)((0.00+Random())*sqrt(A));
    if(T->n==0)
      N0 = 0;
    if( A - T->C[N].N > N0 )
      N0 = A - T->C[N].N;
    Copy_C(&T->C[p],&T->C[N]);

    for(s=0;s<T->Nm;s++)
    {
      for(q=0;q<3;q++)
        b[q] = 1.0;
      b[(int)(Random()*3.0)] = 1.0 + 2.0*T->te*(0.5-Random());
      
      //===== Reshape the outer ~ sqrt(N) shell of the NP =====
      for(i=A-N0,ok=1;i<A;i++)
      {
        L = 2.0*R;
        for(w=0;w<W;w++)
        {
          RAND_VC(a);
          VectorNorm(a);
          for(q=0;q<3;q++)
            T->C[p].X[i][q] = a[q]*R*(1.0+dr*Random());
          
          for(m=0;m<M;m++)
          {
            for(q=0;q<3;q++)
              T->C[p].X[i][q] *= (1.0-dr);
            for(j=0;j<i;j++)
            {
              t = NDR(&T->C[p],i,j);
              if( t < 1.5*(T->Rm[T->C[p].ATMN[i]]+T->Rm[T->C[p].ATMN[j]]) )
              {
                j = T->C[p].N;
                m = M;
                if( t < 1.0*(T->Rm[T->C[p].ATMN[i]]+T->Rm[T->C[p].ATMN[j]]) )
                  ok = 0;
              }
            }
          }
          for(q=0;q<3;q++)
            c[q] = T->C[p].X[i][q]*b[q];
          
          if( VectorLen(c,3) < L )
          {
            L = VectorLen(c,3);
            for(q=0;q<3;q++)
              r[q] = T->C[p].X[i][q];
          }
        }

        for(q=0;q<3;q++)
          T->C[p].X[i][q] = r[q];
      }
      if( ok==1 )
        break;
    }
    if(s==T->Nm)
    {
      sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      Print_LOG(buf);
      exit(1);
    }
    free_i1D(I);
    free_d1D(D);
    sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k,k);
    Print_LOG(buf);
    }
  }
} 
//==================================================================
//     Generate biased Particle 
//==================================================================
void NANO_PACK(Tribe *T, int J, Cell *C, int P)
{
  int i,j,q,m,M,*I,*II,N0,N,ic[3],S,p;
  double a[3],R,*D,*DD,L[3][3],b[3],c[3],B[4][3],**X,t;
  char buf[200];

  for(i=0,R=0.0;i<T->NSPC;i++)
    R += T->Rm[i]*(double)T->SPCN[i]/(double)C->N;
  R *= 1.5*2.0;

  II = make_i1D(C->N);
  DD = make_d1D(C->N);

  N = N0 = 0;
  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    for(i=0;i<3;i++)
      for(q=0;q<3;q++)
	L[i][q] = 0.0;
    S = (int)(3.0*Random());
    if(S==0) //===== bcc =====
    {
      N = 2;
      L[0][0] = L[1][1] = L[2][2] = R*2.0/sqrt(3.0);
      B[0][0] = -0.5; B[0][1] = -0.5; B[0][2] = -0.5;
      B[1][0] =  0.0; B[1][1] =  0.0; B[1][2] =  0.0;
      N0 = 3;
    }
    if(S==1) //===== fcc =====
    {
      N = 4;
      L[0][0] = L[1][1] = L[2][2] = R*sqrt(2.0);
      B[0][0] = -0.5; B[0][1] = -0.5; B[0][2] = -0.5;
      B[1][0] = -0.5; B[1][1] =  0.0; B[1][2] =  0.0;
      B[2][0] =  0.0; B[2][1] = -0.5; B[2][2] =  0.0;
      B[3][0] =  0.0; B[3][1] =  0.0; B[3][2] = -0.5;
      N0 = 2;
    }
    if(S==2) //===== hcp =====
    {
      N = 2;
      L[0][0] =  R;
      L[1][0] = -R/2.0;
      L[1][1] =  R*sqrt(3.0)/2.0;
      L[2][2] =  R*sqrt(8.0/3.0);
      B[0][0] = -1.0/3.0; B[0][1] = -2.0/3.0; B[0][2] = -0.5;
      B[1][0] =  0.0/3.0; B[1][1] =  0.0/3.0; B[1][2] =  0.0;
      N0 = 4;
    }
    M  = (2*N0+1)*(2*N0+1)*(2*N0+1)*N;
    
    X = make_d2D(M,3);
    I = make_i1D(M);
    D = make_d1D(M);
  
    for(i=0;i<3;i++)
      for(q=0;q<3;q++)
	C->L[i][q] = 0.0;
    for(i=0;i<3;i++)
    C->L[i][i] = T->B;
    abc(C);
    
    m = 0;
    for(ic[0]=-N0;ic[0]<=N0;ic[0]++)
    for(ic[1]=-N0;ic[1]<=N0;ic[1]++)
    for(ic[2]=-N0;ic[2]<=N0;ic[2]++)
    {
      for(j=0;j<N;j++,m++)
	for(q=0;q<3;q++)
        {
	  for(i=0,t=0.0;i<3;i++)
	    t += (B[j][i]+(double)ic[i])*L[i][q];
	  X[m][q] = t;
	}
    }
    for(q=0;q<3;q++)
      b[q] = 1.0;
    b[(int)(Random()*3.0)] = 1.0 + 2.0*T->te*(0.5-Random());

    for(m=0;m<M;m++)
    {
      for(q=0;q<3;q++)
	c[q] = X[m][q]*b[q];
      D[m] = VectorLen(c,3);    
    }
    Sort(D,I,M);
    
    RAND_VC(a);
    VectorNorm(a);
    GEN_TYPE(T,C,DD,II);
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->X[II[i]][q] = X[I[i]][q] + a[q]*R;
    
    if(0)
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->X[i][q] += T->ca*RANG();
    if(0)
      NANO_CHOP(T,8,C,p);    
    
    free_d2D(X,M);
    free_i1D(I);
    free_d1D(D);
    Copy_C(C,&T->C[p]);
    sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],S,S);
    Print_LOG(buf);
  }

  ORDER_TYPE(T,J);

  free_i1D(II);
  free_d1D(DD);
}
//==================================================================
//  Mate two parents by swapping core-shells
//==================================================================
void NANO_SWAP(Tribe *T, int J)
{
  int i,m,q,k1,k2,p,N,*TN1,*TN2,qq,po,s,*PM,*SM,**I,N1,qr,o,O;
  double t,r,**R,ar,x,y,dr;
  char buf[200];

  N  = 2*T->N;
  TN1 = make_i1D(T->NSPC);
  TN2 = make_i1D(T->NSPC);
  PM  = make_i1D(2*T->N);
  SM  = make_i1D(2*T->N);
  R   = make_d2D(T->N,T->C[0].N);
  I   = make_i2D(T->N,T->C[0].N);

  for(p=0;p<T->N;p++)
  {
    PM[p+T->N] = 0;
    if(Random()<T->pm)
      PM[p+T->N] = 1;
    SM[p] = 0;
  }

  k1 = k2 = 0;
  dr = 0.02;
  O  = (int)(1.0/dr);

  for(m=0,p=po=T->N;p<2*T->N&&m<T->Nm;p++,m++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {    
    // -----------------------------------------------------------------------
    //   select two parents and mate at least T->Nc times
    // -----------------------------------------------------------------------
    if( m%T->Nc==0 || p>po )
    {
      k1 = (int)(Random()*(double)T->N);
      k2 = (k1+1+(int)(Random()*(double)(T->N-1)))%T->N;
      po = p;
      //===== Sort parents by distance just once =====
      if(SM[k1]==0)
      {
	for(i=0;i<T->C[k1].N;i++)
	  R[k1][i] = VectorLen(T->C[k1].X[i],3);
	Sort(R[k1],I[k1],T->C[k1].N);
	SM[k1] = 1;
      }
      if(SM[k2]==0)
      {
	for(i=0;i<T->C[k2].N;i++)
	  R[k2][i] = VectorLen(T->C[k2].X[i],3);
	Sort(R[k2],I[k2],T->C[k2].N);
	SM[k2] = 1;
      }
    }
    // -----------------------------------------------------------------------
    //   take the core with 40-60% atoms from the first parent
    // -----------------------------------------------------------------------
    for(qq=0;qq<T->NSPC;qq++)
      TN1[qq] = 0;
    N1 = (int)((0.4+0.2*Random())*(double)T->C[k1].N);

    for(i=0;i<N1;i++)
    {
      for(q=0;q<3;q++)
	T->C[p].X[i][q] = T->C[k1].X[I[k1][i]][q];
      TN1[T->C[k1].ATMN[I[k1][i]]]++;
      T->C[p].ATMN[i] =  T->C[k1].ATMN[I[k1][i]];
      T->C[p].ATMZ[i] =  T->C[k1].ATMZ[I[k1][i]];
    }
    
    // -----------------------------------------------------------------------
    //   take the shell with 60-40% atoms from the second parent
    // -----------------------------------------------------------------------
    for(qq=0;qq<T->NSPC;qq++)
      TN2[qq] = 0;
    for(i=N1;i<T->C[k2].N;i++)
      TN2[T->C[k2].ATMN[I[k2][i]]]++;

    // -----------------------------------------------------------------------
    //   check if the composition are correct
    // -----------------------------------------------------------------------
    for(qq=0;qq<T->NSPC;qq++)
      if( TN1[qq]+TN2[qq] != T->C[k1].SPCN[qq] )
	break;

    // -----------------------------------------------------------------------
    //   if the composition is fine copy the second slice into C[p]
    // -----------------------------------------------------------------------
    if( qq==T->NSPC )
    {
      t = R[k2][I[k2][N1]]/R[k1][I[k1][N1-1]];
      for(i=N1;i<T->C[k2].N;i++)
      {
	// -----------------------------------------------------------------------
	//   Rescale the radii for the second parent to avoid core-shell overlap
	// -----------------------------------------------------------------------
	r = t*( R[k2][I[k2][i]]+T->C[k2].R0 ) / R[k2][I[k2][i]];
	for(q=0;q<3;q++)
	  T->C[p].X[i][q] = T->C[k2].X[I[k2][i]][q]*r;
	T->C[p].ATMN[i] =  T->C[k2].ATMN[I[k2][i]];
        T->C[p].ATMZ[i] =  T->C[k2].ATMZ[I[k2][i]];
      }
      T->C[p].N = T->C[k1].N;
      // -----------------------------------------------------------------------
      //   rotate the shell along x, y, or z by a random angle
      // -----------------------------------------------------------------------
      Copy_C(&T->C[p],&T->C[N]);      
      for(s=0;s<T->Ns;s++)
      {
	qr = (int)(Random()*3.0);
	ar = 2.0*Pi*Random();
	for(i=N1;i<T->C[p].N;i++)
	{
	  x = T->C[p].X[i][(qr+1)%3];
	  y = T->C[p].X[i][(qr+2)%3];
	  T->C[p].X[i][(qr+1)%3] =  x*cos(ar) + y*sin(ar);
	  T->C[p].X[i][(qr+2)%3] = -x*sin(ar) + y*cos(ar);
	}

	if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1)
	{
	  for(i=N1;i<T->C[p].N;i++)
	    for(o=0;o<O;o++)
	    {
	      for(q=0;q<3;q++)
		T->C[p].X[i][q] *= (1.0-dr);
	      if(CHCK_Rm(&T->C[p],T->Rm,1.0)==0)
	      {
		for(q=0;q<3;q++)
		  T->C[p].X[i][q] /= (1.0-dr);
		break;
	      }
	    }
	  break;
	}
	// -----------------------------------------------------------------------
	//   If distances are bad restore the offspring and re-try
	// -----------------------------------------------------------------------
	Copy_C(&T->C[N],&T->C[p]);
      }
      if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1) 
      {
	T->P1[p] = T->C[k1].P;
	T->P2[p] = T->C[k2].P;
	sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k1,k2);
	Print_LOG(buf);
      }
      else
	p--;    // unfit distances 
    }
    else
      p--;      // unfit composition
  }
  
  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
  // -----------------------------------------------------------------------
  //   order atoms by type
  // -----------------------------------------------------------------------
  ORDER_TYPE(T,J);

  free_d2D(R,T->N);
  free_i2D(I,T->N);
  free_i1D(PM);
  free_i1D(SM);
  free_i1D(TN1);
  free_i1D(TN2);
}
//==================================================================
//  Mate two parents by swapping two halves
//==================================================================
void NANO_MATE(Tribe *T, int J)
{
  int i,j,m,q,k1,k2,p,o,N,*TN1,*TN2,qq,po,u,s,*PM;
  double a[3],b[3],c[3];
  char buf[200];

  N  = 2*T->N;
  TN1 = make_i1D(T->NSPC);
  TN2 = make_i1D(T->NSPC);
  PM  = make_i1D(2*T->N);

  for(p=T->N;p<2*T->N;p++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    PM[p] = 0;
    if(Random()<T->pm)
      PM[p] = 1;
  }

  k1 = k2 = 0;
  for(m=0,p=po=T->N;p<2*T->N&&m<T->Nm;p++,m++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {    
    // -----------------------------------------------------------------------
    //   select two parents and mate at least T->Nc times
    // -----------------------------------------------------------------------
    if( m%T->Nc==0 || p>po )
    {
      k1 = (int)(Random()*(double)T->N);
      k2 = (k1+1+(int)(Random()*(double)(T->N-1)))%T->N;
      po = p;
      //      sprintf(buf,"%8d %3d %3d %d %2d %2d\n",m,p,po,PM[p],k1,k2);
      //      Print_LOG(buf);
    }
    // -----------------------------------------------------------------------
    //   slice the first parent
    // -----------------------------------------------------------------------
    for(qq=0;qq<T->NSPC;qq++)
      TN1[qq] = 0;

    RAND_VC(a);
    VectorNorm(a);

    for(q=0;q<3;q++)
      a[q] = a[q]*T->C[p].R0;
    
    RAND_VC(b);
    VectorNorm(b);

    for(i=j=0;i<T->C[k1].N;i++)
    {
      for(q=0;q<3;q++)
	c[q] = T->C[k1].X[i][ q] - a[q];
      
      if(DotProd(c,b,3)<0.0)
      {
	for(q=0;q<3;q++)
	  T->C[p].X[j][ q] = T->C[k1].X[i][ q];
	TN1[T->C[k1].ATMN[i]]++;
	T->C[p].ATMN[j++] =  T->C[k1].ATMN[i];
      }
    }
    T->C[p].N = j;
    
    if(1)
    {
      // -----------------------------------------------------------------------
      //   slice the second parent and try to have proper N
      // -----------------------------------------------------------------------
      for(o=-T->No/2;o<=T->No/2;o++)
      { 
	for(qq=0;qq<T->NSPC;qq++)
	  TN2[qq] = 0;
	for(i=j=0;i<T->C[k2].N;i++)
	{
          for(q=0;q<3;q++)
            c[q] = T->C[k2].X[i][q] - a[q] - b[q]*0.2*(double)o/(double)T->No;

	  if(DotProd(c,b,3)>0.0)
	  {
	    TN2[T->C[k2].ATMN[i]]++;
	    j++;
	  }
	}
	// -----------------------------------------------------------------------
	//   check if N and the composition are correct
	// -----------------------------------------------------------------------
	if( T->C[p].N+j==T->C[k2].N )
	{
	  for(qq=0;qq<T->NSPC;qq++)
	    if( TN1[qq]+TN2[qq] != T->C[k1].SPCN[qq] )
	      break;
	  if( qq==T->NSPC )
	    break;
	  else
	    j = 0;
	}
      }
      // -----------------------------------------------------------------------
      //   if the composition is fine copy the second slice into C[p]
      // -----------------------------------------------------------------------
      if(j>0&&T->C[p].N+j==T->C[k2].N)
      {
	for(i=j=0;i<T->C[k2].N;i++)
	{
          for(q=0;q<3;q++)
            c[q] = T->C[k2].X[i][q] - a[q] - b[q]*0.2*(double)o/(double)T->No;
	  if(DotProd(c,b,3)>0.0)
	  {
	    for(q=0;q<3;q++)
	      T->C[p].X[T->C[p].N+j][q] = T->C[k2].X[i][q];
	    T->C[p].ATMN[T->C[p].N+j] =  T->C[k2].ATMN[i];
	    j++;
	  }
	}
	j = T->C[p].N;
	T->C[p].N = T->C[k1].N;

	// -----------------------------------------------------------------------
	//   shift the part from parent 2 randomly to try to have good NN distances
	//   even if the distances are unfit proceed to mutation then check again
	// -----------------------------------------------------------------------
	Copy_C(&T->C[p],&T->C[N]);
	for(s=0;s<T->Ns;s++)
	{
	  RAND_VC(a);
	  VectorNorm(a);
	  for(i=j;i<T->C[p].N;i++)
	    for(q=0;q<3;q++)
	      T->C[p].X[i][q] += 0.7*a[q]*T->C[p].R0;

	  if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1)
	    break;
	  Copy_C(&T->C[N],&T->C[p]);
	}
	// -----------------------------------------------------------------------
	//   mutate with probability T->pm and swap with probability T->ps
	// -----------------------------------------------------------------------
	s = u = -1;
	if(PM[p])
	{
	  for(u=0;u<T->Nu;u++)
	  {
	    Copy_C(&T->C[p],&T->C[N]);
	    if( SHKE_CL(&T->C[N],T->ml,T->C[p].R0*T->ma)!=0 )  //SHKE_CL already made not to shift
	    {
	      if(T->NSPC>1)
		for(i=s=0;i<T->C[N].N;i++)
		  if(Random()<T->ps)
	          {
		    for(j=i;T->C[N].ATMN[i]==T->C[N].ATMN[j%T->C[N].N];j++);
		    j = ( j + (int)(Random()*(double)(T->C[N].N-T->SPCN[T->C[N].ATMN[i]])) )%T->C[N].N;  // swap only different species
                    for(q=0;q<3;q++)
                      dSwap(&T->C[N].X[i][q],&T->C[N].X[j][q]);
                    s++;
                  }
	      if(CHCK_Rm(&T->C[N],T->Rm,1.0)==1)
	      {
		Copy_C(&T->C[N],&T->C[p]);
		break;
	      }
	    }
	  }
	}
	if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1) 
        {
	  T->P1[p] = T->C[k1].P;
	  T->P2[p] = T->C[k2].P;
          sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k1,k2);
          Print_LOG(buf);
	}
	else
	  p--;    // unfit distances 
      }
      else
	p--;      // unfit composition
    }
    else
      p--;        // unfit volume: always fine for NPs
  
  }
  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
  // -----------------------------------------------------------------------
  //   order atoms by type
  // -----------------------------------------------------------------------
  ORDER_TYPE(T,J);

  free_i1D(PM);
  free_i1D(TN1);
  free_i1D(TN2);

}
//==================================================================
//     Rubic's cube mutation
//==================================================================
void NANO_RUBE(Tribe *T, int J)
{
  int i,j,m,q,k,k1,k2,p,o,N,*TN1,qr,po,u,s,*PM;
  double a[3],b[3],Rm,x,y,ar,R;
  char buf[200];

  N  = 2*T->N;
  TN1 = make_i1D(T->NSPC);
  PM  = make_i1D(2*T->N);

  for(p=T->N;p<2*T->N;p++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    PM[p] = 0;
    if(Random()<T->pm)
      PM[p] = 1;
  }
  for(m=0,Rm=0.0;m<T->NSPC;m++)
    if( Rm < T->Rm[m] )
      Rm = T->Rm[m];           //determing the minimum hard sphere radius

  k1 = k2 = 0;
  for(m=0,p=po=T->N;p<2*T->N&&m<T->Nm;p++,m++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    // -----------------------------------------------------------------------
    //   select one parent at least T->Nc times
    // -----------------------------------------------------------------------
    R = 0.0;
    if( m%T->Nc==0 || p>po )
    {
      k1 = k2 = (int)(Random()*(double)T->N);
      po = p;
      for(i=0,R=0.0;i<T->C[k1].N;i++)
	if( VectorLen(T->C[k1].X[i],3) > R )
	  R = VectorLen(T->C[k1].X[i],3);

      //      sprintf(buf,"%8d %3d %3d %d %2d %2d\n",m,p,po,p,k1,k2);
    }
    // -----------------------------------------------------------------------
    //   slice the parent
    // -----------------------------------------------------------------------
    for(o=0;o<=T->No;o++)
    {    
      Copy_C(&T->C[k1],&T->C[N]);
      Copy_C(&T->C[k1],&T->C[p]);   //defining two C's one as the offspring and one for messing around
      // ---------------------------------------------------------------------
      //  select an axis for Rubik's rotation and a random shift along it
      // ---------------------------------------------------------------------

      qr = (int)(Random()*3.0);   //values are 0,1,2 os selecting axis a, b, or c
      b[qr] = (0.5-Random())*R;

      // ---------------------------------------------------------------------
      //  Do Rubik's rotation, only by 180 degrees for 3D but any for 0D
      // ---------------------------------------------------------------------

      ar = (0.25+1.50*Random())*Pi;                 // select random angle between 45 and 315

      //===== Rotate NP to put a random atom k along z =====
      k = (int)(Random()*(double)T->C[p].N);
      ROT_CELL(&T->C[p],k,2,1);
      ROT_CELL(&T->C[p],k,0,2);

      for(i=0;i<T->C[p].N;i++)
      if(T->C[p].X[i][qr]>b[qr])
      {
	x = T->C[p].X[i][(qr+1)%3];
	y = T->C[p].X[i][(qr+2)%3];
	T->C[p].X[i][(qr+1)%3]  =  x*cos(ar) + y*sin(ar);
	T->C[p].X[i][(qr+2)%3]  = -x*sin(ar) + y*cos(ar);
	T->C[p].X[i][(qr+0)%3] += Rm;
      }
      Copy_C(&T->C[p],&T->C[N]);
      // -----------------------------------------------------------------------
      //   shift the rotated part randomly to try to have good NN distances
      //   even if the distances are unfit proceed to mutation then check again
      // -----------------------------------------------------------------------
      for(s=0;s<T->Ns;s++)
      {
	RAND_VC(a);
	VectorNorm(a);
	for(i=0;i<T->C[p].N;i++)
	  if(T->C[p].X[i][qr]>b[qr])
	    for(q=0;q<3;q++)
	      T->C[p].X[i][q] += 0.5*a[q]*Rm;
	if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1)
	  break;
	Copy_C(&T->C[N],&T->C[p]);           //If distances are bad overwrite offspring with playing and re-try
      }
      // -----------------------------------------------------------------------
      //   mutate with probability T->pm and swap with probability T->ps
      // -----------------------------------------------------------------------
      s = u = -1;
      if(PM[p])
      {
	for(u=0;u<T->Nu;u++)
	{
	  Copy_C(&T->C[p],&T->C[N]);
	  if( SHKE_CL(&T->C[N],T->ml,T->C[p].R0*T->ma)!=0 )       //This function changes the lattice vectors and positions of atoms randomly
	  {
	    if(T->NSPC>1)
	      for(i=s=0;i<T->C[N].N;i++)
		if(Random()<T->ps)                            //if probability is met (and there is more than 1 atom type) then swap atoms
		{
		  for(j=i;T->C[N].ATMN[i]==T->C[N].ATMN[j%T->C[N].N];j++);
		  j = ( j + (int)(Random()*(double)(T->C[N].N-T->SPCN[T->C[N].ATMN[i]])) )%T->C[N].N;  // swap only different species
		  for(q=0;q<3;q++)
		    dSwap(&T->C[N].X[i][q],&T->C[N].X[j][q]);
		  s++;
		}
	    if(CHCK_Rm(&T->C[N],T->Rm,1.0)==1)                 //if distances are good copy plating to offspring and move on
	    {
	      Copy_C(&T->C[N],&T->C[p]);
	      break;
	    }
	  }
	}
      }
      // -----------------------------------------------------------------------
      //   check distances
      // -----------------------------------------------------------------------
      if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1)
      {
	T->P1[p] = T->C[k1].P;
	T->P2[p] = T->C[k2].P;
	sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k1,k2);
	Print_LOG(buf);
	break;            // distances are good, leave To loop
      }
    }
    if(o==T->No+1)
      p--;                  // unfit distances
  }
  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
  // -----------------------------------------------------------------------
  //   order atoms by type
  // -----------------------------------------------------------------------
  ORDER_TYPE(T,J);
  
  free_i1D(PM);
  free_i1D(TN1);

}
//==================================================================
//     Inversion
//==================================================================
void NANO_SYMM(Tribe *T, int J)
{
  int i,j,jj,m,q,k1,k2,p,o,N,*TN1,qq,po,M,U;
  double b[3],c[3],d[3],t,Rm,R;
  char buf[200];

  N  = 2*T->N;                             //Set N to twice as many as C's in Tribe 
  U  = N + 1;
  TN1 = make_i1D(T->NSPC);                //Vector of length = # of species

  for(m=0,Rm=0.0;m<T->NSPC;m++)              //Sets minimum radius of the various atom types
    if( Rm < T->Rm[m] )
      Rm = T->Rm[m];

  k1 = k2 = 0;
  for(m=0,p=po=T->N;p<2*T->N&&m<T->Nm;p++,m++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    // -----------------------------------------------------------------------
    //   select one parent at least T->Nc times
    // -----------------------------------------------------------------------
    if( m%T->Nc==0 || p>po )
    {
      k1 = k2 = (int)(Random()*(double)T->N);
      po = p;
    }
    // -----------------------------------------------------------------------
    //   slice the parent
    // -----------------------------------------------------------------------
    for(o=0;o<=T->No;o++)  //loops over the number of tries to find a good slice
    {    
      Copy_C(&T->C[k1],&T->C[N]);
      Copy_C(&T->C[k1],&T->C[U]); // Generate two new C's from the parent one to be offspring and one to mess with
      // ---------------------------------------------------------------------
      //  always plane cut through middle for inversion (but shake atoms to keep slice location (random)
      // ---------------------------------------------------------------------

      t = Random()*Rm;//(double)o/(double)T->No*Rm; //determine length of vector to shake atoms by

      RAND_VC(b);
      VectorNorm(b);
      RAND_VC(d);
      VectorNorm(d);                           //generate two random vectors of length 1
      for(i=0;i<T->C[N].N;i++)                 //shake all atoms by |d|*t
        for(q=0;q<3;q++)
	  T->C[N].X[i][q] += d[q]*t;         // shake with increasing amplitude
     
      for(qq=0;qq<T->NSPC;qq++)                  //Keeping track of number of each type of atom in the slice (initialize)
	TN1[qq] = 0;
      for(i=j=0;i<T->C[N].N;i++)               // for each atom in C[messing around] selects if atom is in slice
      {
	for(q=0;q<3;q++)
	  c[q] = T->C[N].X[i][ q];
	if(DotProd(c,b,3)<0.0)               //check if atom dist to cell center is parallel 
	{
	  for(q=0;q<3;q++)
	    T->C[U].X[j][ q] = T->C[N].X[i][ q];
	  TN1[T->C[N].ATMN[i]]++;                //Keeping track of number of each type of atom in the slice
	  T->C[U].ATMN[j  ] =  T->C[N].ATMN[i];    //sets C[U].ATMN[j] to C[U].ATMN[i] and then j++
	  T->C[U].ATMZ[j++] =  T->C[N].ATMZ[i];
	}
      }
      for(qq=0;qq<T->NSPC;qq++)
	if( TN1[qq] < T->C[k1].SPCN[qq]/2 )      //checks if at least 1/2 of each type of atom is selected //SHOULD I WORRY ABOUT THE FIXED ATOMS THAT AREN"T DUPLICATING?
	  break;
      M = T->C[U].N = j;                       //Number of atoms is slice

      i = 0;
      if( T->C[U].N>=T->C[k1].N/2 && qq==T->NSPC )  // atoms may be merged, hence try >= N/2  //can only be merged if more than 1/2 of the atoms are here (can delete atoms not add)
      {
        for(i=0;i<M;i++)
	{
	  if(J==7)                                                      // reflection
	  {
	    for(q=0;q<3;q++)
	      c[q] = T->C[U].X[i][q];
	    t = DotProd(c,b,3);
	    for(q=0;q<3;q++)
	      T->C[U].X[T->C[U].N][q] = c[q] - 2.0*b[q]*t;	    
	  }	    
          if(J==8)                                                      // inversion
	    for(q=0;q<3;q++)
	      T->C[U].X[T->C[U].N][q] = -T->C[U].X[i][q];      
	  
	  T->C[U].ATMN[T->C[U].N  ] = T->C[U].ATMN[i];
	  T->C[U].ATMZ[T->C[U].N++] = T->C[U].ATMZ[i];
	}
	for(i=0;i<M;i++)                        //for all atoms in original slice
        {
	  for(j=jj=M,R=T->B;jj<T->C[U].N;jj++)
	  {
	    t = NDR(&T->C[U],i,jj);
	    if( t < R )
	    {
	      R = t;
	      j = jj;
	    }		
	  }
	  t = R/(T->Rm[T->C[U].ATMN[i]]+T->Rm[T->C[U].ATMN[j]]);  //dist between i and j / (minimum distant between the two atoms)
	  if( T->C[U].ATMN[i]==T->C[U].ATMN[j] && t<1.0 )                // merge the two atoms if they are too close together
          {
            for(q=0;q<3;q++)
              T->C[U].X[i][q] = (T->C[U].X[i][q]+T->C[U].X[j][q])*0.5;
            T->C[U].ATMN[j] = 99;
	  }
	}
	for(i=j=0;j<T->C[U].N;i++,j++)         //deletes atoms with type 99 and shifts array up
        { 
	  for(q=0;q<3;q++)
	    T->C[U].X[i][q] = T->C[U].X[j][q];
	  T->C[U].ATMN[i] = T->C[U].ATMN[j];
	  T->C[U].ATMZ[i] = T->C[U].ATMZ[j];
	  if( T->C[U].ATMN[j] == 99 )
	    i--;
	}	  
	T->C[U].N = i;                        //updates number of atoms
        // -----------------------------------------------------------------------
        //   check if N and the composition are correct
        // -----------------------------------------------------------------------
        if( T->C[U].N==T->C[k1].N )
        {
	  for(qq=0;qq<T->NSPC;qq++)
	    TN1[qq] = 0;
	  for(i=0;i<T->C[U].N;i++)
	    TN1[T->C[U].ATMN[i]]++;
          for(qq=0;qq<T->NSPC;qq++)
            if( TN1[qq] != T->C[k1].SPCN[qq] )
              break;
          if( qq<T->NSPC )
            i = 0;
        }
      }
      // -----------------------------------------------------------------------
      //   if the composition is fine check distances
      // -----------------------------------------------------------------------
      if(i>0&&T->C[U].N==T->C[k1].N)
      {
        if(CHCK_Rm(&T->C[U],T->Rm,1.0)==1)
        {
          T->P1[p] = T->C[k1].P;
          T->P2[p] = T->C[k2].P;
	  Copy_C(&T->C[U],&T->C[p]);
	  printf("%s %3d\n",T->NES[J],p);
	  sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k1,k2);
	  Print_LOG(buf);
	  break;            // composition and distances are good, leave To loop
        }
      }
    }
    if(o==T->No+1)
      p--;                  // unfit composition or distances
  }
  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
  // -----------------------------------------------------------------------
  //   order atoms by type
  // -----------------------------------------------------------------------
  ORDER_TYPE(T,J);

  free_i1D(TN1);

}
//==================================================================
//     Evolve Tribe
//==================================================================
void NANO_MUTE(Tribe *T, int J)
{
  int i,j,k,p,m,q,s,o;
  double r,*f,*D,R,V,b[3];
  char buf[200];

  s = 0;
  f = make_d1D(T->N+1);
  D = make_d1D(T->C[0].N);

  for(i=0,V=0.0;i<T->NSPC;i++)
    V += pow(T->Rm[i],3.0)*(double)T->SPCN[i];
  R = pow(V,1.0/3.0)*1.5;
  for(q=0;q<3;q++)
    b[q] = 1.0;

  for(m=0,p=T->N;p<2*T->N&&m<T->Nm;p++,m++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    r = Random();
    for(j=1,f[0]=0.0;j<T->N+1;j++)
      f[j] = f[j-1] + T->f[j-1];
    for(j=0;j<T->N+1;j++)
      f[j] /= f[T->N];
    for(k=0;k<T->N-1&&f[k+1]<r;k++); //selects random k

    Copy_C(&T->C[k],&T->C[p]);

    if( SHKE_CL(&T->C[p],T->cl,T->C[p].R0*T->ca)==0 )
      p--;
    else
    {
      if(T->NSPC>1)
      {
	for(i=0;i<T->C[p].N;i++)
	  D[i] = VectorLen(T->C[p].X[i],3);
	for(i=s=0;i<T->C[p].N;i++)
	  if(Random()<T->cs)
	  {
	    for(j=i;T->C[p].ATMN[i]==T->C[p].ATMN[j%T->C[p].N];j++);	    
	    for(o=0;o<T->C[p].N;o++)
	    {
	      j = ( j + (int)(Random()*(double)(T->C[p].N-T->SPCN[T->C[p].ATMN[i]])) )%T->C[p].N;
	      if( fabs(D[i]-D[j]) < 0.5*(T->Rm[T->C[p].ATMN[i]]+T->Rm[T->C[p].ATMN[j]]) )
	      {
		for(q=0;q<3;q++)
		  dSwap(&T->C[p].X[i][q],&T->C[p].X[j][q]);
		s++;
		break;
	      }
	    }
	  }
      }
      if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1||ADJT_NP(&T->C[p],T->Rm,R,b,3))
      {
	T->P1[p] = T->P2[p] = T->C[k].P;
	sprintf(buf,"%3d %3d %s %3d %3d %3d swaps\n",T->n,p,T->NES[J],k,k,s);
	Print_LOG(buf);
      }
      else
	p--;
    }
    sprintf(T->C[p].TAG,"%3d %3d",T->n,p);
  }
  free_d1D(f);
  free_d1D(D);

  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
}
//==================================================================
//     Generate crystals or films using randomized initial positions
//==================================================================
void BULK_RAND(Tribe *T, int J, Cell *C, int P)
{
  int    p,N,q,S[30][30],II[30],JJ[30],MIN[5],MAX[5];;
  double RR[30][30],AA[30][30],ZZ[30][30],Rc[5];
  char   buf[200],s[200],w[3][200],*t;
  FILE   *in;

  N = 2*T->N;

  READ_BLOCK(C,II,JJ,S,RR,AA,ZZ,Rc,MIN,MAX);

  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    //===== generate random structures =====
    if(T->JS==0)
    {
      RAND_CL(T,&T->C[p],&T->C[N],0,Rc,MIN,MAX);
      SHRT_LV(&T->C[p]);
    }
    //===== randomize specified structures =====
    else if(T->JS==1)
    {
      sprintf(buf,"INI/POSCAR%03d",0);
      READ_CELL(&T->C[N],buf);
      TEMP_CL(T,&T->C[N],p);
    }
    //===== use specified LVs and generate atoms randomly =====
    else if(T->JS==2 && T->ND==3)
    {
      sprintf(buf,"INI/POSCAR%03d",0);
      READ_CELL(&T->C[p],buf);
      RAND_CL(T,&T->C[p],&T->C[N],1,Rc,MIN,MAX);
    }
    else if(T->JS==3 && T->ND==3)
    {
      sprintf(buf,"INI/POSCAR%03d",0);
      READ_CELL(&T->C[N],buf);
      RAND_SYMM(T,&T->C[N],p);
    }
  }
}
//==================================================================
void CELL_PLOT(Cell *C, double s)
{
  int   i,q,m,n,**X,N,M;
  char  buf[200];

  return;

  N = 50;
  M = 100;

  Relative(C);
  sprintf(buf,"sleep %6.2lf; clear",s);
  system(buf);
  printf("% 10.4lf  % 10.4lf\n",C->L[0][0],C->L[2][2]);
  for(m=0;m<M;m++)
    printf("-");
  printf("\n\n");
  for(n=0;n<N;n++,printf("\n"))
  {
    for(m=0;m<M;m++)
    {
      sprintf(buf," ");
      for(i=0;i<C->N;i++)
	if( (int)((double)N*(1.0-C->X[i][2]))==n && (int)((double)M*C->X[i][0])==m )
	  sprintf(buf,"%d",C->ATMN[i]);
//      printf("%s%s%s",(C->ATMN[0] == 1) ? RED : (C->ATMN[0] == 2) ? GREEN : (C->ATMN[0] == 3) ? YELLOW : RESET, buf, RESET);
      printf("%s",buf);
    }
  }

  for(m=0;m<M;m++)
    printf("-");
  printf("\n\n");
  Real(C);
}
//==================================================================
//     Generate random crystals using TETRIS algorithm
//==================================================================
void BULK_TETR(Tribe *T, int J, Cell *C, int P)
{
  int    i,j,q,s,w,W,n,m,M,*I,p,Z,k,S[30][30],K,Q,h,ii,o,F,O,ww,*B,*U,*Y,II[30],JJ[30],MIN[5],MAX[5];
  double t,a[3],A,V,*D,L,R,dr,r[3],b[3],c[3],X[3],Rm,z,RR[30][30],AA[30][30],ZZ[30][30],rr[30][3],Rc[5],x[3],oo[30][3];
  char   buf[200];

  I = make_i1D(C->N);
  U = make_i1D(C->N);
  Y = make_i1D(C->N);
  B = make_i1D(C->N);
  D = make_d1D(C->N);

  Q = READ_BLOCK(C,II,JJ,S,RR,AA,ZZ,Rc,MIN,MAX);

  for(i=0,V=0.0;i<T->NSPC;i++)    
    V += (4.0/3.0*Pi)*pow(T->Rm[i]/0.7,3.0)*(double)T->SPCN[i];
  for(i=0,Rm=0.0;i<T->NSPC;i++)
    Rm += T->Rm[i]/0.7*(double)T->SPCN[i];
  Rm /= (double)C->N;

  dr = 0.02;

  for(p=T->N;p<2*T->N;p++)
  if( p==P || ( (P<0) && (p>=T->SES[J]) && (p<T->FES[J] ) ) )
  {
    printf("%s %3d\n",T->NES[J],p);
    Z = 0;
    z = 1000.0;
    for(s=0;s<T->Nm;s++)
    {
      K = GEN_TETR(T,C,D,I,B,U,II,S,Q);

      if(T->JS==0)
      {
	RAND_LV(C);
	t = pow( 2.0*T->VOL*(1.0+0.25*(0.5-Random())*2.0),1.0/3.0);
	SCALE_LATTICE(C,t);
	SHRT_LV(C);
	Lat_Align(C);
	abc(C);
      }
      if(T->JS==2)
	READ_CELL(C,"INI/POSCAR000");

      //===== atoms at drop point should not see neighbors along the c-axis =====      
      for(q=0;q<3;q++)
	C->L[2][q] *= 3.0;
      R = C->L[2][2];
      Reciprocal(C);
      
      CrossProd(C->L[0],C->L[1],a);                // find the unit cell base area
      W = 5*(1+(int)(VectorLen(a,2)/(Pi*Rm*Rm)));  // 5 times the number of atoms per base
      M = (int)(R*0.5/dr)-1; 

      for(i=0;i<C->N;i++)
	for(q=0;q<3;q++)
	  C->X[i][q] = 0.0;

      for(n=0;n<K;n++)                             // loop over blocks
      {
	L = C->L[2][2];	
	if( II[B[n]]>1 && JJ[B[n]]>0 )
	  ROT_BLOCK(II[B[n]],RR[B[n]],AA[B[n]],ZZ[B[n]],S[B[n]]);
	for(w=j=0;w<W;w++)
	{
	  X[0] = Random();
	  X[1] = Random();
	  X[2] = 0.5+0.1*Random();                 // start at no lower than 1.5 = 3.0*0.5 fractional height

	  for(m=0;m<M;m++)
	  {
	    X[2] -= dr/R;
            for(o=O=0;o<36;o++)
            {
              for(ii=0;ii<II[B[n]];ii++)
              {
                i = I[U[n]+ii];
                for(q=0;q<D3;q++)
                  for(k=0,C->X[i][q]=0.0;k<D3;k++)
                    C->X[i][q] += X[k]*C->L[k][q];
                C->X[i][0] += RR[B[n]][ii]*cos((AA[B[n]][ii]+5.0*(double)o)*Pi/180.0);
                C->X[i][1] += RR[B[n]][ii]*sin((AA[B[n]][ii]+5.0*(double)o)*Pi/180.0);
		C->X[i][2] += ZZ[B[n]][ii];

		//===== put the atom in the cell, only lateral shifts are needed =====
		for(q=0;q<3;q++)
		  for(k=0,x[q]=0.0;k<3;k++)
		    x[q] += C->X[i][k]*C->R[q][k];
		for(q=0;q<2;q++)
		{
		  while(x[q] <  0.0)
		    x[q] += 1.0;
		  while(x[q] >= 1.0)
		    x[q] -= 1.0;
                  for(k=0,C->X[i][q]=0.0;k<3;k++)
                    C->X[i][q] += x[k]*C->L[k][q];
		}
              }
              for(ii=F=0;ii<II[B[n]];ii++)
              {
                i = I[U[n]+ii];
                if( U[n]+ii > 0 )
                  for(j=U[n]+ii-1;j>=0;j--)
                  {
                    if( NDB(C,i,I[j]) < 1.2*(T->Rm[C->ATMN[i]]+T->Rm[C->ATMN[I[j]]]) )
                    {
                      F = 1;                         // if at least one distance is too short, go to next angle o 
                      j = -1;
                      ii = II[B[n]];
                    }
                  }
              }
              if( F==0 || (m==0&&o==0) )             // save positions of atoms in block n that do not crash into others 
              {
		for(ii=0;ii<II[B[n]];ii++)
		  for(q=0;q<D3;q++)
		    oo[ii][q] = C->X[I[U[n]+ii]][q];;
		if(F==0)
		  O = 1;
              }
              if( II[B[n]]==1 )                     // do not rotate if there is just one atom 
                break;
	    }

            if( O==0 )                             // if distances too short at any angle, freeze z
	      if(1)                                // if 0, allow blocks to explore full range of z
		break;
	  }
	  if( oo[0][2] < L )                       // keep the block with the lowest z coordinate
	  {
	    L = oo[0][2];
	    for(ii=0;ii<II[B[n]];ii++)
	      for(q=0;q<D3;q++)
		rr[ii][q] = oo[ii][q];
	  }
	}
	for(ii=0;ii<II[B[n]];ii++)
          for(q=0;q<D3;q++)
            C->X[I[U[n]+ii]][q] = rr[ii][q];
      }

      if( CHCK_Rm(C,T->Rm,0.95)==1 && CHCK_BOND(C,Rc,MIN,MAX)==1 )
      {
	L = C->L[2][2];
	for(q=0;q<3;q++)
	  C->L[2][q] *= 0.99;
	while( CHCK_Rm(C,T->Rm,0.95)==1 && CHCK_BOND(C,Rc,MIN,MAX)==1 )
	  for(q=0;q<3;q++)
	    C->L[2][q] *= 0.99;
	for(q=0;q<3;q++)
	  C->L[2][q] /= 0.99;
	L = C->L[2][2]/L;

	if( L < z )
	{
	  Copy_C(C,&T->C[2*T->N]);
	  z = L;
	}
	//===== accept cell with c axis no longer than the original or best out of 10 =====
	if( z < 1.0/3.0 || Z==10 )
	{
	  Z = 0;
	  L = z;
	  z = 1000.0;
	  break;
	}
	Z++;
      }
    }

    Copy_C(&T->C[2*T->N],C);

    if(s==T->Nm)
    {
      sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
      Print_LOG(buf);
      exit(1);
    }

    //===== shorten the c lattice constant to remove gap =====
    Copy_C(C,&T->C[p]);
    if( ADJT_CL(&T->C[p],T->Rm,10)==1 )         // It was not adjusted before to caluclate the c-axis mismatch easily
    {
      JAR(&T->C[p]);
      ORDER_Z(&T->C[p]);
      abc(C);
      t = ( C->LAT[0]<C->LAT[1] ) ? C->LAT[0] : C->LAT[1];
      sprintf(buf,"%3d %3d %s %3d %3d %3.2lf %3.3lf %3.3lf\n",T->n,p,T->NES[J],-1,-1,C->LAT[2]/t,L,CELL_VOL(&T->C[p])/(double)T->C[p].N);
      Print_LOG(buf);
    }
    else
    {
      printf("WARNING in TETRIS: cell discarded because the interatomic distances are too short\n");
      p--;
    }
  }
  ORDER_TYPE(T,J);

  free_i1D(I);
  free_i1D(U);
  free_i1D(Y);
  free_i1D(B);
  free_d1D(D);

}
//==================================================================
//     Evolve Tribe for 2D or 3D crystals
//==================================================================
void BULK_MATE(Tribe *T, int J)
{
  int i,j,m,q,k,k1,k2,p,o,N,*TN1,*TN2,qq,po,u,s,*PM,Q;
  double a[3],b[3],c[3],V,frac;
  char buf[200];

  N  = 2*T->N;
  TN1 = make_i1D(T->NSPC);
  TN2 = make_i1D(T->NSPC);
  PM  = make_i1D(2*T->N);

  for(p=T->N;p<2*T->N;p++)
    PM[p] = 0;
  for(p=T->N;p<2*T->N;p++)
    if(Random()<T->pm)
      PM[p] = 1;
  for(k=0;k<T->N;k++)
    abc(&T->C[k]);
  for(k=0;k<T->N;k++)
    Relative(&T->C[k]);

  k1 = k2 = 0;
  for(m=0,p=po=T->N;p<2*T->N&&m<T->Nm;p++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {    
    T->C[p].XT = 0;
    // -----------------------------------------------------------------------
    //   select two parents and mate at least T->Nc times
    // -----------------------------------------------------------------------
    if( m%T->Nc==0 || p>po )
    {
      k1 = (int)(Random()*(double)T->N);
      k2 = (k1+1+(int)(Random()*(double)(T->N-1)))%T->N;
      po = p;
      sprintf(buf,"%8d %3d %3d %d %2d %2d\n",m,p,po,PM[p],k1,k2);
      Print_LOG(buf);
    }
    // -----------------------------------------------------------------------
    //   apply a random shift to all atoms in each parent
    // -----------------------------------------------------------------------
    for(q=0;q<3;q++)
    {
      a[q] = Random();
      b[q] = Random();
    }
    for(i=0;i<T->C[k1].N;i++)
      for(q=0;q<3;q++)
      {
	T->C[k1].X[i][q] += a[q];
	T->C[k2].X[i][q] += b[q];
	if( T->C[k1].X[i][q]>1.0 )
	  T->C[k1].X[i][q] -= 1.0;
	if( T->C[k2].X[i][q]>1.0 )
          T->C[k2].X[i][q] -= 1.0;
      }
    // -----------------------------------------------------------------------
    //   slice the first parent
    // -----------------------------------------------------------------------
    for(qq=0;qq<T->NSPC;qq++)
      TN1[qq] = 0;

    RAND_VC(a);
    if(T->ND==3)     // for crystal select center near middle
      for(q=0;q<3;q++)
	a[q] = 0.5 + a[q]*0.25;
    if(T->ND==2)     // for film select center near middle
      a[2] = 0.5;

    RAND_VC(b);
    if(T->ND==2)
      b[2] = 0.0;
    VectorNorm(b);

    for(i=j=0;i<T->C[k1].N;i++)
    {
      Q = 0;
      for(q=0;q<3;q++)
	c[q] = T->C[k1].X[i][ q] - a[q];
      
      if(DotProd(c,b,3)<0.0)
	Q = 1;
      if(Q==1)
      {
	for(q=0;q<3;q++)
	  T->C[p].X[j][ q] = T->C[k1].X[i][ q];
	TN1[T->C[k1].ATMN[i]]++;
	T->C[p].ATMN[j++] =  T->C[k1].ATMN[i];
      }
    }
    T->C[p].N = j;
    // -----------------------------------------------------------------------
    //   mix lattice vectors, discard if the volume is too different
    // -----------------------------------------------------------------------
    frac = (double)T->C[p].N/(double)T->C[k1].N;
    MATE_LT(&T->C[k1],&T->C[k2],&T->C[p],frac);
    V = (CELL_VOL(&T->C[k1])+CELL_VOL(&T->C[k2]))*0.5;
    
    if( fabs(CELL_VOL(&T->C[p])-V)/V < 0.1 )
    {
      // -----------------------------------------------------------------------
      //   slice the second parent and try to have proper N
      // -----------------------------------------------------------------------
      for(o=-T->No/2;o<=T->No/2;o++)
      { 
	for(qq=0;qq<T->NSPC;qq++)
	  TN2[qq] = 0;
	for(i=j=0;i<T->C[k2].N;i++)
	{
	  Q = 0;
          for(q=0;q<3;q++)
            c[q] = T->C[k2].X[i][q] - a[q] - b[q]*0.2*(double)o/(double)T->No;

	  if(DotProd(c,b,3)>0.0)
	    Q = 1;
	  if(Q==1)
	  {
	    TN2[T->C[k2].ATMN[i]]++;
	    j++;
	  }
	}
	// -----------------------------------------------------------------------
	//   check if N and the composition are correct
	// -----------------------------------------------------------------------
	if( T->C[p].N+j==T->C[k2].N )
	{
	  for(qq=0;qq<T->NSPC;qq++)
	    if( TN1[qq]+TN2[qq] != T->C[k1].SPCN[qq] )
	      break;
	  if( qq==T->NSPC )
	    break;
	  else
	    j = 0;
	}
      }
      // -----------------------------------------------------------------------
      //   if the composition is fine copy the second slice into C[p]
      // -----------------------------------------------------------------------
      if(j>0&&T->C[p].N+j==T->C[k2].N)
      {
	for(i=j=0;i<T->C[k2].N;i++)
	{
          Q = 0;
          for(q=0;q<3;q++)
            c[q] = T->C[k2].X[i][q] - a[q] - b[q]*0.2*(double)o/(double)T->No;
	  if(DotProd(c,b,3)>0.0)
	    Q = 1;
          if(Q==1)
	  {
	    for(q=0;q<3;q++)
	      T->C[p].X[T->C[p].N+j][q] = T->C[k2].X[i][q];
	    T->C[p].ATMN[T->C[p].N+j] =  T->C[k2].ATMN[i];
	    j++;
	  }
	}
	j = T->C[p].N;
	T->C[p].N = T->C[k1].N;
	abc(&T->C[p]);
	// -----------------------------------------------------------------------
	//   shift the part from parent 2 randomly to try to have good NN distances
	//   even if the distances are unfit proceed to mutation then check again
	// -----------------------------------------------------------------------
	Copy_C(&T->C[p],&T->C[N]);
	for(s=0;s<T->Ns;s++)
	{
	  RAND_VC(a);
	  T->C[p].XT = 0;
	  if(T->ND==2)   // because ND=2 the atoms are *constrained* within the plane
	    a[2] = 0.0;
	  for(i=j;i<T->C[p].N;i++)
	    for(q=0;q<3;q++)
	      T->C[p].X[i][q] += 0.5*a[q]*T->C[p].R0/T->C[p].LAT[q];

	  Real(&T->C[p]);
	  if( CHCK_Rm(&T->C[p],T->Rm,1.0)==1 || ADJT_CL(&T->C[p],T->Rm,T->C[p].N/4)==1 )
	  {
	    Relative(&T->C[p]);
	    break;
	  }
	  Relative(&T->C[p]);
	  Copy_C(&T->C[N],&T->C[p]);
	}
	// -----------------------------------------------------------------------
	//   mutate with probability T->pm and swap with probability T->ps
	// -----------------------------------------------------------------------
	s = u = -1;
	if(PM[p])
	{
	  V = CELL_VOL(&T->C[p]);
	  for(u=0;u<T->Nu;u++)
	  {
	    Copy_C(&T->C[p],&T->C[N]);
	    Real(&T->C[N]);
	    if( SHKE_CL(&T->C[N],T->ml,T->C[p].R0*T->ma)!=0 )  //SHKE_CL already made not to shift
	    if( SHRT_LV(&T->C[N])!=0 )                         //SHRT_LV only deals with lattice - should be fine
	    if( CELL_VOL(&T->C[N]) > 0.5*V )
	    {
	      if( T->JS!=2 && !(T->ml>(10^-10) || T->cl>(10^-10)) )
	      {
		if(T->ND==3)
		  SCALE_Cell(&T->C[N], pow(CELL_VOL(&T->C[N])/V,-1.0/3.0) + 0.1*Random());
		if(T->ND==2)
		  SCALE_Cell(&T->C[N], pow(CELL_VOL(&T->C[N])/V,-1.0/2.0) + 0.1*Random());
	      }
	      if(T->NSPC>1)
		for(i=s=0;i<T->C[N].N;i++)
		  if(Random()<T->ps)
	          {
		    for(j=i;T->C[N].ATMN[i]==T->C[N].ATMN[j%T->C[N].N];j++);
		    j = ( j + (int)(Random()*(double)(T->C[N].N-T->SPCN[T->C[N].ATMN[i]])) )%T->C[N].N;  // swap only different species
                    for(q=0;q<3;q++)
                      dSwap(&T->C[N].X[i][q],&T->C[N].X[j][q]);
                    s++;
                  }
	      if( CHCK_Rm(&T->C[N],T->Rm,1.0)==1 || ADJT_CL(&T->C[N],T->Rm,T->C[p].N/4)==1 )
	      {
		Copy_C(&T->C[N],&T->C[p]);
		break;
	      }
	    }
	  }
	  if(T->C[p].XT==1)
	    Relative(&T->C[p]);
	}
	Real(&T->C[p]);	
	if(CHCK_Rm(&T->C[p],T->Rm,1.0)==1) 
        {
	  T->P1[p] = T->C[k1].P;
	  T->P2[p] = T->C[k2].P;
	  Relative(&T->C[p]);
          printf("%s %3d\n",T->NES[J],p);
          sprintf(buf,"%3d %3d %s %3d %3d % lf\n",T->n,p,T->NES[J],k1,k2,frac);
          Print_LOG(buf);
	}
	else
	  p--;    // unfit distances 
      }
      else
	p--;      // unfit composition
    }
    else
      p--;        // unfit volume
    m++;
  }
  sprintf(buf,"% 5d mating   tries\n",m);
  Print_LOG(buf);
  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
  // -----------------------------------------------------------------------
  //   order atoms by type
  // -----------------------------------------------------------------------
  for(p=T->N;p<2*T->N;p++)
  if( p>=T->SES[J] && p<T->FES[J] )
    for(qq=i=0;qq<T->NSPC-1;qq++)
      for(j=0;j<T->C[p].SPCN[qq];j++)
      {
	for(k=i;k<T->C[p].N;k++)
	  if(T->C[p].ATMN[k]==qq)
	    break;
	if(k==T->C[p].N)
	{
	  sprintf(buf,"ERROR: atom types are wrong\n");
	  fprintf(stderr,"ERROR: atom types are wrong\n");
	  Print_LOG(buf);
	  exit(1);
	}
	for(q=0;q<D3;q++)
	  dSwap(&T->C[p].X[i][q],&T->C[p].X[k][q]);
	iSwap(&T->C[p].ATMN[i],&T->C[p].ATMN[k]);
	i++;
      }

  for(k=0;k<T->N;k++)
    Real(&T->C[k]);

  for(p=T->N;p<2*T->N;p++)   //cleaning up atomic positions
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    Real(&T->C[p]);
    SHRT_LV(&T->C[p]);
    JAR(&T->C[p]);
    LIST(&T->C[p],0);
  }

  free_i1D(PM);
  free_i1D(TN1);
  free_i1D(TN2);
}
//==================================================================
//     Evolve Tribe
//==================================================================
void BULK_MUTE(Tribe *T, int J)
{
  int i,j,k,p,m,q,s;
  double V,r,*f;
  char buf[200];
  long SEED;

  s = 0;
  f = make_d1D(T->N+1);

  for(m=0,p=T->N;p<2*T->N&&m<T->Nm;p++)
  if( p>=T->SES[J] && p<T->FES[J] )
  {
    r = Random();
    for(j=1,f[0]=0.0;j<T->N+1;j++)
      f[j] = f[j-1] + T->f[j-1];
    for(j=0;j<T->N+1;j++)
      f[j] /= f[T->N];
    for(k=0;k<T->N-1&&f[k+1]<r;k++); //selects random k

    Copy_C(&T->C[k],&T->C[p]);
    V = CELL_VOL(&T->C[p]);

    if( SHKE_CL(&T->C[p],T->cl,T->C[p].R0*T->ca)==0 || SHRT_LV(&T->C[p])==0 || fabs( CELL_VOL(&T->C[p])/V-1.0 ) > 0.5 )
      p--;
    else
    {
      if( T->JS!=2 && !(T->ml>(10^-10) || T->cl>(10^-10)) )
      {
	if(T->ND==3)
	  SCALE_Cell(&T->C[p], pow(CELL_VOL(&T->C[p])/V,-1.0/3.0) + 0.1*(0.5-Random()));
	if(T->ND==2)
	  SCALE_Cell(&T->C[p], pow(CELL_VOL(&T->C[p])/V,-1.0/2.0) + 0.1*(0.5-Random()));
      }
      if(T->NSPC>1)
	for(i=s=0;i<T->C[p].N;i++)
	  if(Random()<T->cs)
	  {
	    for(j=i;T->C[p].ATMN[i]==T->C[p].ATMN[j%T->C[p].N];j++);
	    j = ( j + (int)(Random()*(double)(T->C[p].N-T->SPCN[T->C[p].ATMN[i]])) )%T->C[p].N;  // swap only different species
	    for(q=0;q<3;q++)
	      dSwap(&T->C[p].X[i][q],&T->C[p].X[j][q]);
	    s++;
	  }
      if( CHCK_Rm(&T->C[p],T->Rm,1.0)==1 || ADJT_CL(&T->C[p],T->Rm,T->C[p].N/4)==1 )
      {
	T->P1[p] = T->P2[p] = T->C[k].P;
	printf("%s %3d\n",T->NES[J],p);
	sprintf(buf,"%3d %3d %s %3d %3d\n",T->n,p,T->NES[J],k,k);
	Print_LOG(buf);
      }
      else
	p--;
    }
    GetSeed(&SEED);
    sprintf(T->C[p].TAG,"%3d %3d",T->n,p);
    sprintf(buf,"mutated %5d with the %ld seed\n",p,SEED);
    m++;
  }
  sprintf(buf,"% 5d mutation tries out of % 5d\n",m,T->Nm);
  Print_LOG(buf);
  free_d1D(f);

  if(m==T->Nm)
  {
    sprintf(buf,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    fprintf(stderr,"EXIT in %s: exceeded %3d tries\n",T->NES[J],T->Nm);
    Print_LOG(buf);
    exit(1);
  }
}
//==================================================================
