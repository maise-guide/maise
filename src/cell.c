#include "cell.h"

//=========================================================================
//
//=========================================================================
void Build_Cell(Cell *C, int J)
{
  int i,q;

  C->N += 2;

  C->NH   = 1000;
  C->VH   = make_i1D(C->NH);
  C->R1   = make_d2D(C->N,C->NM);
  C->R2   = make_d2D(C->N,C->NM);
  C->ATMN = make_i1D(C->N);
  C->ATMZ = make_i1D(C->N);
  C->L    = make_d2D(D3,D3);                 // lattice vectors
  C->R    = make_d2D(D3,D3);                 // lattice vectors 
  C->X    = make_d2D(C->N,D3);               // coordinates
  C->W    = make_d2D(C->N,D3);               // Wyckoff positions
  C->V    = make_d2D(C->N+3,D3);             // velocities
  C->F    = make_d2D(C->N+3,D3);             // forces: atoms + lattice vectors
  C->U    = make_d1D(6);                     // stresses
  C->Nn   = make_i1D(C->N);                  // number of neighbors for C->Rc
  C->nn   = make_i1D(C->N);                  // number of neighbors for C->rc 
  C->BC   = make_d1D(C->N);  
  C->EA   = make_d1D(C->N);
  C->MM   = make_i1D(C->N);                  // initial magnetic moment
  C->FRC  = make_i1D(C->N);                  // mask for training/testing NN
  C->FF   = make_i2D(C->N,D3);               // fixed: atoms + lattice vectors
  C->min  = make_d1D(D3);

  C->max  = make_d1D(D3);
  C->Ni   = make_i2D(C->N,C->NM);
  C->Nj   = make_i2D(C->N,C->NM);            // to avoid considering i-j-i triplets
  C->MNT  = 9;
  C->NRDF = 3000;
  C->RDF  = make_d3D(C->NRDF,C->MNT,C->MNT); // max number of species
  C->S    = make_d3D(C->N,C->NM,D3);
  C->SG   = make_d2D(200,12);
  C->SL   = make_i1D(2);
  C->NDX  = make_d2D(C->N,C->NM);
  C->EVOK = 0;

  //===== needed for NNs, not needed for potentials or the ES =====
  if(J==1&&C->MODT==1)
  {
    C->ndx = make_d3D(C->N,C->NM,C->NM);
    C->cos = make_d3D(C->N,C->NM,C->NM);
    C->fc  = make_d2D(C->N,C->NM);
  }
  C->R0   = 1.0;
  for(i=0;i<C->N;i++)
    C->MM[i] = 0;
  for(i=0;i<C->N;i++)
    C->ATMN[i] = 1;
  for(i=0;i<C->N;i++) 
    for(q=0;q<D3;q++) 
      C->FF[i][q] = 1; 
  for(i=0;i<C->N;i++)
    for(q=0;q<D3;q++)
      C->X[i][q] = C->V[i][q] = C->F[i][q] = 0.0;
  for(i=0;i<D3;i++) 
    for(q=0;q<D3;q++) 
      C->L[i][q] = 0.0; 

  C->N -= 2;

  if(C->EVOK==0)
  {
    C->ev   = make_d1D(C->N*D3);
    C->EV   = make_d2D(C->N*D3,C->N*D3);
    C->EVOK = 1;
  }

  return;
}
//=========================================================================
// finds NM nearest neighbors in 3D case
//=========================================================================
void LIST(Cell *C, int O)
{
  int i,j,k,q,q1,j1;
  int ic[3],N[3];
  double r,x,*R,a[3],RC;

  RC = C->Rc*C->Rc;

  if( C->ND>0 )
    JAR(C);
  N[0] = N[1] = N[2] = 0;
  Reciprocal(C);
  for(q=0;q<C->ND;q++)
    N[q] = ceil( C->Rc*VectorLen(C->R[q],3) );
  #pragma omp parallel private(j1,R,j,ic,q,r,q1,x,k,a) num_threads(C->NP)
  {
    R  = make_d1D(C->NM);
    #pragma omp for schedule(dynamic,C->NB)  
    for(i=0;i<C->N;i++)
    {
      C->Nn[i] = C->nn[i] = 0;
      for(j1=0;j1<C->NM;j1++)
      {
        R[j1] = 10000000000.0;
        C->Ni[i][j1] = 0;
        C->NDX[i][j1] = 0.0;
      }
      for(j=0;j<C->N;j++)
        for(ic[0]=-N[0];ic[0]<=N[0];ic[0]++)
        for(ic[1]=-N[1];ic[1]<=N[1];ic[1]++)
        for(ic[2]=-N[2];ic[2]<=N[2];ic[2]++)
          if(! (i==j&&ic[0]==0&&ic[1]==0&&ic[2]==0) )
          {
            for(q=0,r=0.0;q<D3;q++)
            {
              for(q1=0,x=0.0;q1<D3;q1++) 
                x += (double)(ic[q1]) *C->L[q1][q]; 
              x += C->X[j][q]-C->X[i][q]; 
              r += x*x; 
            }

            if(r<RC)
            {
              R[C->Nn[i]] = r;
              C->Ni[i][C->Nn[i]] = j;
              for(q=0;q<D3;q++)
              {
                for(q1=0,x=0.0;q1<D3;q1++)
                  x += (double)( ic[q1] ) *C->L[q1][q];
                C->S[i][C->Nn[i]][q] = x;
              }
	      C->Nn[i]++;
              if( C->Nn[i]==C->NM )
              {
	        if(C->JOBT/10==0)
		  fprintf(stderr,"Increase 'Max number of nearest neighbors' % lf % lf %3d\n",sqrt(r),C->Rc,C->NM);
	        else
		  fprintf(stderr,"Increase 'Max number of nearest neighbors' in setup % lf % lf\n",sqrt(r),C->Rc);
                SAVE_CELL(C,"CONTCAR",0);
                exit(1);
              }
            }
          }
      //===== to order all nearest neighbors, N^2 cost  =====
      if(O==1)
      {
        for(j=0;  j<C->Nn[i];j++)
	  for(k=j+1;k<C->Nn[i];k++)
	    if(R[k]<R[j])
	    {                                                                                                                                                                         
	      atomSwap(C,i,j,k);
	      dSwap(&R[j],&R[k]);
	    }                                
      //===== to put the nearest neighbor as C->Ni[i][0] =====
      }
      else
      {
        for(j=1;  j<C->Nn[i];j++)
	  if(R[j]<R[0])
	  {
	    atomSwap(C,i,j,0);
	    dSwap(&R[j],&R[0]);
	  }
      }
      C->nn[i] = C->Nn[i];

      for(j=0;j<C->Nn[i];j++)
        C->NDX[i][j] = sqrt(R[j]);

      if(C->JOBT/10>1&&C->MODT==1)
      for(j=0;j<C->Nn[i];j++)
        for(k=0;k<C->Nn[i];k++)
        {
          for(q=0,r=0.0;q<D3;q++)
            r += DX(C,i,j,q)*DX(C,i,k,q);
          C->cos[i][j][k] = r/(C->NDX[i][j]*C->NDX[i][k]);
          for(q=0;q<3;q++)
            a[q] = DX(C,i,j,q)-DX(C,i,k,q);
          C->ndx[i][j][k] = sqrt(DotProd(a,a,3));
        }
    }
    free_d1D(R);
  }
}
//========================================================================= 
// Gives distance between atoms i and j; no boudnary conditions !!!!
//=========================================================================
double NDR(Cell *C, int i, int j)
{
  int q;
  double x,r;

  for(r=0.0,q=0;q<D3;q++)
  {
    x = C->X[j][q] - C->X[i][q];
    r += x*x;
  }
  return sqrt(r);

}
//=========================================================================
// Gives q-distance between atom i and its n'th nearest neighbor: X[q][.n.] - X[q][i]
//=========================================================================
double DX(Cell *C, int i, int n, int q)
{
  return ( C->X[C->Ni[i][n]][q] - C->X[i][q] + C->S[i][n][q]); 
}
//=========================================================================                                                       
// Gives distance between atom i's neighbors j and k                                                                              
//=========================================================================                                                       
double ndx(Cell *C, int i, int j, int k)
{
  int q;
  double r[3];

  for(q=0;q<3;q++)
    r[q] = DX(C,i,j,q)-DX(C,i,k,q);
  return sqrt(DotProd(r,r,3));
}
//=========================================================================
// Gives distance between atom i and its .n.'th nearest neighbor: Ri.n.
//=========================================================================
double NDX(Cell *C, int i, int n)
{
  int q;
  double r,x;

  for(r=0.0,x=0.0,q=0;q<D3;q++)
  {
    x  = C->X[C->Ni[i][n]][q] - C->X[i][q] + C->S[i][n][q];
    r += x*x;
  }
  return sqrt(r);
}
//=========================================================================
//  Returns either periodic or non-periodic distance
//=========================================================================
double NDD(Cell *C, int i, int n)
{
  if( C->ND==0 )
    return NDR(C,i,n);
  return NDX(C,i,n);
}
//=========================================================================
// Gives cos of the angle formed by atom i and its nearest neighbors .j. and .k.
//=========================================================================
double Cos(Cell *C, int i, int j, int k)
{
  double t;
  int q;

  t = 0;
  for(q=0;q<D3;q++)
    t += DX(C,i,j,q)*DX(C,i,k,q);

  t /= ( NDX(C,i,j)*NDX(C,i,k) );

  return t;
}
//=========================================================================
// Return the neighbor number of atom j on the neighbor list of atom i
//========================================================================
int NN(Cell *C, int i, int j)
{
  int i1;

  for(i1=0;i1<C->Nn[i];i1++)
    if(C->Ni[i][i1]==j)
      return i1;
  return C->Nn[i];

}
//=========================================================================
void SCALE_Cell(Cell *C, double a)
{
  int i,q;
  for(i=0;i<C->N;i++)
    for(q=0;q<C->ND;q++)
      C->X[i][q] *= a;
  for(i=0;i<C->ND;i++)
    for(q=0;q<C->ND;q++)
      C->L[i][q] *= a;
  for(i=0;i<C->ND;i++)
    for(q=0;q<C->ND;q++)
      C->R[i][q] /= a;
}
//=========================================================================
void SCALE_LATTICE(Cell *C, double a)
{
  int i,q;
  for(i=0;i<3;i++)
    for(q=0;q<C->ND;q++)
      C->L[i][q] *= a;
}
//=========================================================================
void MATCH_LATTICE(Cell *C, Cell *D)
{
  int    i,q;
  double a;

  Relative(C);
  for(i=0;i<3;i++)
  {
    a = sqrt( (D->L[i][0]*D->L[i][0]+D->L[i][1]*D->L[i][1]+D->L[i][2]*D->L[i][2])/(C->L[i][0]*C->L[i][0]+C->L[i][1]*C->L[i][1]+C->L[i][2]*C->L[i][2]) );
    for(q=0;q<3;q++)
      C->L[i][q] *= a;
  }
  Real(C);
}
//==================================================================
void Print_RDF(Cell *C, char *file)
{
  int k,m,n;
  double dr;
  FILE *out;

  dr   = C->Rmax/(double)(C->NRDF);

  out = fopen(file,"w");

  for(k=0;k<C->NRDF;k++)
  {
    fprintf(out,"% lf   ",dr*(double)k);
    for(m=0;m<C->NSPC;m++)
      for(n=m;n<C->NSPC;n++)
	fprintf(out,"% lf  ",C->RDF[k][m][n]);
    fprintf(out,"\n");
  }
  fclose(out);

}
//==================================================================
void Print_RDF_FULL(Cell *C, char *file)
{
  int k,m,n;
  double dr,s;
  FILE *out;

  dr   = C->Rmax/(double)(C->NRDF);

  out = fopen(file,"w");

  fprintf(out,"#    R         total ");
  for(m=0;m<C->NSPC;m++)
    for(n=m;n<C->NSPC;n++)
      fprintf(out,"      %c%c   ",(char)(m+65),(char)(n+65));
  fprintf(out,"\n");

  for(k=0;k<C->NRDF;k++)
    {
      fprintf(out,"% lf   ",dr*(double)k);

      s = 0.0;
      for(m=0;m<C->NSPC;m++)
	for(n=m;n<C->NSPC;n++)
	  s += C->RDF[k][m][n];
      fprintf(out,"% lf  ",s);

      for(m=0;m<C->NSPC;m++)
	for(n=m;n<C->NSPC;n++)
	  fprintf(out,"% lf  ",C->RDF[k][m][n]);
      fprintf(out,"\n");
    }
  fclose(out);

}
//==================================================================
void RDF(Cell *C, int J)
{
  int i,j,k,m,n,o,Do;
  double r, dr, x, DR, f;

  dr   = C->Rmax/(double)(C->NRDF);
  Do   = (int)(3.0*C->DR/dr);       // go 3 sigmas left and right
  DR   = 0.5/(C->DR*C->DR);

  for(k=0;k<C->NRDF;k++)
    for(m=0;m<C->NSPC;m++)
      for(n=0;n<C->NSPC;n++)
        C->RDF[k][m][n] = 0.0;
  for(i=0;i<C->N;i++)
  {
    for(j=0;j<C->Nn[i];j++)
      if(C->ATMN[C->Ni[i][j]]>=C->ATMN[i])
      {
	r = NDX(C,i,j);
	o = (int)(r/dr);
	for(k=o-Do;k<=o+Do;k++)
	  if(k>=0&&k<C->NRDF)
	  {
	    x = (double)k*dr;
	    if(x<C->Rmin)
	      f = 1.0;
	    else
	      f = cos( 0.5*Pi*(x-C->Rmin)/(C->Rmax-C->Rmin) );
	    C->RDF[k][C->ATMN[i]][C->ATMN[C->Ni[i][j]]] += f*exp(-(r-x)*(r-x)*DR);
	  }
      }
  }
  for(k=0;k<C->NRDF;k++)
    for(m=0;m<C->NSPC;m++)
      C->RDF[k][m][m] *= 0.5;
  if(J==0)
    return;

  r = 0.0;
  for(k=0;k<C->NRDF;k++)
    for(m=0;m<C->NSPC;m++)
      for(n=m;n<C->NSPC;n++)
        r += C->RDF[k][m][n]*C->RDF[k][m][n];
  r = 1.0/sqrt(r);
  for(k=0;k<C->NRDF;k++)
    for(m=0;m<C->NSPC;m++)
      for(n=m;n<C->NSPC;n++)
        C->RDF[k][m][n] *= r;
}
//==================================================================
double CxC(Cell *C, Cell *D)
{
  int k,m,n;
  double r;

  r = 0.0;
  for(k=0;k<C->NRDF;k++)
    for(m=0;m<C->MNT;m++)
      for(n=m;n<C->MNT;n++)
	r += C->RDF[k][m][n]*D->RDF[k][m][n];
  return r;
}
//==================================================================
void Reciprocal(Cell *C)
{
  int i,q;
  double t;

  VectorProd(C->L[0],C->L[1],C->R[2]);
  VectorProd(C->L[1],C->L[2],C->R[0]);
  VectorProd(C->L[2],C->L[0],C->R[1]);
  t = 1.0/CrossProd(C->L[0],C->L[1],C->L[2]);
  for(i=0;i<3;i++)
    for(q=0;q<3;q++)
      C->R[i][q] *= t;
}
//==================================================================
void Relative(Cell *C)
{
  int i,j,q;
  double x[3];

  if(C->XT==0)
  {
    fprintf(stderr,"Cell is already relative\n");
    exit(0);
  }
  Reciprocal(C);  

  for(i=0;i<C->N;i++)
  {
    for(q=0;q<D3;q++)
      for(j=0,x[q]=0.0;j<D3;j++)
	x[q] += C->X[i][j]*C->R[q][j];
    for(q=0;q<D3;q++)
      C->X[i][q] = x[q];
  }
  C->XT = 0;
}
//==================================================================
void Real(Cell *C)
{
  int i,j,q;
  double x[3];

  if(C->XT==1)
  {
    fprintf(stderr,"Cell is already real\n");
    exit(0);
  }
  for(i=0;i<C->N;i++)
  {
    for(q=0;q<D3;q++)
      for(j=0,x[q]=0.0;j<D3;j++)
	x[q] += C->X[i][j]*C->L[j][q];
    for(q=0;q<D3;q++)
      C->X[i][q] = x[q];
  }
  C->XT = 1;
}
//======================================================================
void JAR(Cell *C)
{
  int i,q;
 
  if(C->ND==0)
    return;

  Relative(C);

  for(i=0;i<C->N;i++)
    for(q=0;q<D3;q++)
    {
      while(C->X[i][q] <  0.0)
	C->X[i][q] += 1.0;
      while(C->X[i][q] >= 1.0)
        C->X[i][q] -= 1.0;
    }
  Real(C);
}
//======================================================================
// Transforms Lattice Vectors to make them shorter as in JPCM, 20, 064210
// there is a typo in the paper, should be ceil(|a*b|/b^2)
// so far implemented only projections on a,b,c, not along the diagonals
//======================================================================
int SHRT_LV(Cell *C)
{
  int i,q,j,J;
  double P[3],a;

  abc(C);

  J = 1;  
  for(i=0;i<3;i++)
    for(j=0;j<2;j++)
    {
      if(J==1)
      {
	abc(C);
	for(q=0;q<3;q++)
	  P[q] = DotProd(C->L[(q+1)%3],C->L[(q+2)%3],3);
	J = 0;
      }
      if( fabs(P[(i+2)%3])/C->LAT[(i+1-j)%3] > C->LAT[(i+1-j)%3]*0.5 )
      {
	a = ceil(fabs(P[(i+2)%3])/(C->LAT[(i+1-j)%3]*C->LAT[(i+1-j)%3]));
	for(q=0;q<3;q++)
	  C->L[(i+j)%3][q] -= a*sign(P[(i+2)%3])*C->L[(i+1-j)%3][q];
	J = 1;
      }
    }
  if(J==1)
    abc(C);
  if( C->LAT[0]<C->R0 || C->LAT[1]<C->R0 || C->LAT[2]<C->R0 )
    return 0;

  JAR(C);

  return 1;
}
//======================================================================
// Orders the atoms by z vector
//======================================================================
void ORDER_Z(Cell *C)
{
  int i,j,q,k,K;

  K = C->ATMN[C->N-1]+1;

  Relative(C);
  for(k=i=0;k<K;k++)
    for(;i<C->N&&C->ATMN[i]==k;i++)
    {
      for(j=i+1;j<C->N&&C->ATMN[j]==k;j++)
        if(C->X[j][2]<C->X[i][2])
          for(q=0;q<3;q++)
          {
            dSwap(&C->X[i][q],&C->X[j][q]);
            iSwap(&C->FF[i][q],&C->FF[j][q]);
          }
   }
   Real(C);
}
//======================================================================
// Orders the atoms type
//======================================================================
void ORDER(Cell *C)
{
  int i,j,q;

  for(i=0;i<C->N-1;i++)
    for(j=i+1;j<C->N;j++)
      if(C->ATMN[j]<C->ATMN[i])
      {
        iSwap(&C->ATMN[j],&C->ATMN[i]);
        iSwap(&C->ATMZ[j],&C->ATMZ[i]);
        for(q=0;q<D3;q++)
        {
          dSwap(&C->X[i][q],&C->X[j][q]);
          iSwap(&C->FF[i][q],&C->FF[j][q]);
        }
      }

  for(j=0;j<C->NSPC;j++)
    C->SPCN[j]=0;

  for(i=0;i<C->N;i++)
    for(j=0;j<C->NSPC;j++)
      if(C->ATMN[i]==j)
        C->SPCN[j]++;
}
//======================================================================
void ADD(Cell *C, Cell *D, double x, double y, double z)
{
  int i,q;
  for(i=0;i<D->N;i++)
  {
    C->X[i+C->N][0] = D->X[i][0] + x;
    C->X[i+C->N][1] = D->X[i][1] + y;
    C->X[i+C->N][2] = D->X[i][2] + z;
    C->ATMN[i+C->N] = D->ATMN[i];
    C->Nn[i+C->N]   = D->Nn[i];
    for(q=0;q<3;q++)
      C->FF[i+C->N][q] = D->FF[i][q];
  }
  C->N += D->N;
}
//======================================================================
