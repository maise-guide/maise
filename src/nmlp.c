#include "nmlp.h"

PAR *OOO;

//================================================================
// ANN functions
//================================================================
double GPc(int G, double x)
{
  if(G==0) return 0.0;
  if(G==1) return 2.0*x*(x*x-1.0);
  return 0.0;
}
//================================================================
double GPp(int G, double x)
{
  if(G==0) return 1.0;
  if(G==1) return 1.0-x*x;
  return 0.0;
}
//================================================================
double GP(int G, double x)
{
  if(G==0) return x;
  if(G==1) return tanh(x);
  return 0.0;
}
//======================================================
double FRC_ANN_PARA(ANN *R, LNK *L, double ***e, double ***d)
{
  int i,j,k,n,m,l,q,nij,sij,spc;
  double tm,tn,tl;

  for(i=0,L->e=0.0;i<L->N;i++)
  {
    spc = L->ATMN[i];
    for(n=0;n<R->NU[0];n++)
      e[i][0][n] = (L->Cn[i][n]-R->Rmin[spc][n]*R->DR[spc][n])-0.0;
    for(k=0;k<R->NL-1;k++)
      for(m=0;m<R->NU[k+1];m++)
      {
        for( n=0, e[i][k+1][m]=R->B[spc][k][m]; n<R->NU[k]; n++)
          e[i][k+1][m] += e[i][k][n]*R->W[spc][k][m][n];
        e[i][k+1][m] =  GP(R->GT[k+1],e[i][k+1][m]);
        d[i][k+1][m] = GPp(R->GT[k+1],e[i][k+1][m]);
      }
    L->EA[i] = e[i][R->NL-1][0];
    L->e += e[i][R->NL-1][0];
  }
  
  for(q=0;q<6;q++)
    L->s[q] = 0.0;
  for(i=0;i<L->N;i++)    
  if(L->MRK[i]==1)
  {
    spc = L->ATMN[i];
    for(q=0;q<3;q++)
      L->f[i][q] = 0.0;
    // for MLP with 2 hidden layers
    // =======  dGni/dxi  ========
    for(m=0;m<R->NU[2];m++)
    {
      tm = R->W[spc][2][0][m]*d[i][2][m];
      for(n=0;n<R->NU[1];n++)
      {
        tn = tm*R->W[spc][1][m][n]*d[i][1][n];
        for(l=0;l<R->NU[0];l++)
        {
	  tl = tn*R->W[spc][0][n][l];
	  for(q=0;q<3;q++)
	    L->f[i][q] -= L->Fn[i][L->DNn[i]][l][q]*tl;
	  // ======= will not work if not all atoms are marked =======
	  if(R->EFS>1)
	    for(q=0;q<6;q++)
	      L->s[q] -= L->Sn[i][l][q]*tn*R->W[spc][0][n][l];
        }
      }
    }
    // =======  dGnj/dxi  ========
    for(j=0;j<L->DNn[i];j++)
    {
      sij = L->DNs[i][j];
      nij = L->DNi[i][j];
      for(m=0;m<R->NU[2];m++)
      {
	tm = R->W[sij][2][0][m]*d[nij][2][m];
	for(n=0;n<R->NU[1];n++)
	{
	  tn = tm*R->W[sij][1][m][n]*d[nij][1][n];
	  for(l=0;l<R->NU[0];l++)
	  {
	    tl = tn*R->W[sij][0][n][l];
	    for(q=0;q<3;q++)
	      L->f[i][q] -= L->Fn[i][j][l][q]*tl;
	  }
	}
      }
    }
  }

  return L->e;
}
//================================================================
double DF_ANN_PARA(ANN *R, LNK *L, double ***Bp, double ****Wp,  double ***e,  double ***d,  double ***c, double ***xn, double ***xm)
{
  int i,k,n,m,q,l,j,nij,sij,spc;
  double E,tm,tn,tl,tmn,rmn,tnl,rnl[3],dfq[3],dnl[3];
  
  for(i=0,E=0.0;i<L->N;i++)
  {
    spc = L->ATMN[i];
    for(n=0;n<R->NU[0];n++)
      e[i][0][n] = (L->Cn[i][n]-R->Rmin[spc][n]*R->DR[spc][n])-0.0;
    for(k=0;k<R->NL-1;k++)
      for(m=0;m<R->NU[k+1];m++)
      {
        for( n=0, e[i][k+1][m]=R->B[spc][k][m]; n<R->NU[k]; n++)
          e[i][k+1][m] += e[i][k][n]*R->W[spc][k][m][n];
        e[i][k+1][m] =  GP(R->GT[k+1],e[i][k+1][m]);
        d[i][k+1][m] = GPp(R->GT[k+1],e[i][k+1][m]);
        c[i][k+1][m] = GPc(R->GT[k+1],e[i][k+1][m]);
      }
    E += e[i][R->NL-1][0];
  }

  //===== calculate forces here first =====
  for(i=0;i<L->N;i++)
  if(L->MRK[i]==1)
  {
    spc = L->ATMN[i];
    for(q=0;q<3;q++)
      L->f[i][q] = 0.0;
    // ======= dGni/dxi ========
    for(m=0;m<R->NU[2];m++)
    {
      tm = R->W[spc][2][0][m]*d[i][2][m];
      for(n=0;n<R->NU[1];n++)
      {
	tn = tm*R->W[spc][1][m][n]*d[i][1][n];
	for(l=0;l<R->NU[0];l++)
	{
	  tl = tn*R->W[spc][0][n][l];
	  for(q=0;q<3;q++)
	    L->f[i][q] -= L->Fn[i][L->DNn[i]][l][q]*tl;
	}
      }
    }
    // ======= dGnj/dxi ========
    for(j=0;j<L->DNn[i];j++)
    {
      sij = L->DNs[i][j];
      nij = L->DNi[i][j];
      for(m=0;m<R->NU[2];m++)
      {
	tm = R->W[sij][2][0][m]*d[nij][2][m];
	for(n=0;n<R->NU[1];n++)
	{
	  tn = tm*R->W[sij][1][m][n]*d[nij][1][n];
	  for(l=0;l<R->NU[0];l++)
	  {
	    tl = tn*R->W[sij][0][n][l];
	    for(q=0;q<3;q++)
	      L->f[i][q] -= L->Fn[i][j][l][q]*tl;
	  }
	}
      }
    }
  }

  //===== now calculate weight derivatives =====
  for(i=0;i<L->N;i++)
  if(L->MRK[i]==1)
  {
    spc = L->ATMN[i];
    for(q=0;q<3;q++)
      dfq[q] = 2.0*(L->f[i][q]-L->F[i][q])*R->WF*L->W;

    k = R->NL-2;
    //===== calculate xn to store for backpropagation =====
    for(n=0;n<R->NU[k-1];n++)
      for(j=0;j<=L->DNn[i];j++)
      {
	sij = L->DNs[i][j];
	for(q=0;q<3;q++)
	  xn[n][j][q] = 0.0;
	for(l=0;l<R->NU[k-2];l++)
	  for(q=0;q<3;q++)
	    xn[n][j][q] += R->W[sij][k-2][n][l]*L->Fn[i][j][l][q]*dfq[q];
      }      
    //===== calculate xm to store for backpropagation ===== 
    for(m=0;m<R->NU[k];m++)
      for(j=0;j<=L->DNn[i];j++)
      {
	sij = L->DNs[i][j];
	nij = L->DNi[i][j];
	for(q=0;q<3;q++)
	  xm[m][j][q] = 0.0;
	for(n=0;n<R->NU[k-1];n++)
	{
	  tmn = R->W[sij][k-1][m][n]*d[nij][k-1][n];
	  for(q=0;q<3;q++)
	    xm[m][j][q] += tmn*xn[n][j][q];
	}
      }
    
    if(R->MIX==0)
    {
    for(m=0;m<R->NU[k];m++)
      for(j=0;j<=L->DNn[i];j++)
      {
        sij = L->DNs[i][j];
        nij = L->DNi[i][j];
        tm = R->W[sij][k][0][m]*c[nij][k][m];
        for(n=0;n<R->NU[k-1];n++)
        {
          tmn = tm*R->W[sij][k-1][m][n]*d[nij][k-1][n];
          for(q=0;q<3;q++)
            Bp[sij][k-1][m] -= tmn*xn[n][j][q];
        }
      }

    for(n=0;n<R->NU[k-1];n++)
      for(j=0;j<=L->DNn[i];j++)
      {
        sij = L->DNs[i][j];
        nij = L->DNi[i][j];
        for(m=0;m<R->NU[k];m++)
        {
          tn = R->W[sij][k][0][m]*d[nij][k][m]*R->W[sij][k-1][m][n]*c[nij][k-1][n];
          tm = R->W[sij][k][0][m]*c[nij][k][m]*R->W[sij][k-1][m][n]*d[nij][k-1][n];
          for(q=0;q<3;q++)
            Bp[sij][k-2][n] -= tn*xn[n][j][q] + tm*xm[m][j][q];
        }
      }

    for(m=0;m<R->NU[k];m++)
      for(j=0;j<=L->DNn[i];j++)
      {
        sij = L->DNs[i][j];
        nij = L->DNi[i][j];
        for(n=0;n<R->NU[k-1];n++)
	{
	  tn = d[nij][k][m]*R->W[sij][k-1][m][n]*d[nij][k-1][n];
	  for(q=0;q<3;q++)
	    Wp[sij][k][0][m] -= tn*xn[n][j][q];
	}
      }

    for(m=0;m<R->NU[k];m++)
      for(n=0;n<R->NU[k-1];n++)
        for(j=0;j<=L->DNn[i];j++)
        {
	  sij = L->DNs[i][j];
	  nij = L->DNi[i][j];
          tmn = R->W[sij][k][0][m]*d[nij][k][m]*d[nij][k-1][n];
	  for(q=0;q<3;q++)
	    Wp[sij][k-1][m][n] -= tmn*xn[n][j][q];
	  tmn = R->W[sij][k][0][m]*c[nij][k][m]*e[nij][k-1][n];
	  for(q=0;q<3;q++)
	    Wp[sij][k-1][m][n] -= tmn*xm[m][j][q];
        }
    }
    for(j=0;j<=L->DNn[i];j++)
      for(l=R->O;l<R->NU[0];l++)
      {
	sij = L->DNs[i][j];
	for(n=0;n<R->NU[1];n++)
        {
	  nij = L->DNi[i][j];
	  tnl = d[nij][1][n]*e[nij][0][l];
	  for(q=0;q<3;q++)
	  {
	    dnl[q] = L->Fn[i][j][l][q]*d[nij][1][n]*dfq[q];
	    rnl[q] = c[nij][1][n]*e[nij][0][l]*xn[n][j][q];
	  }
	  for(m=0;m<R->NU[2];m++)
	  {
	    tmn = R->W[sij][2][0][m]*d[nij][2][m]*R->W[sij][1][m][n];
	    rmn = R->W[sij][2][0][m]*c[nij][2][m]*R->W[sij][1][m][n]*tnl;
	    for(q=0;q<3;q++)
	      Wp[sij][0][n][l] -= tmn*(dnl[q]+rnl[q])+rmn*xm[m][j][q];
	  }
	}
      }
    
  }
  return E;
}
//================================================================
//
//================================================================
double DE_ANN_PARA(ANN *R, LNK *L, double ***Bp, double ****Wp,  double ***e)
{
  int i,k,n,m,l,spc;
  double E,c,tm,tn;

  for(i=0,E=0.0;i<L->N;i++)
  {
    spc = L->ATMN[i];
    for(n=0;n<R->NU[0];n++)
      e[i][0][n] = (L->Cn[i][n]-R->Rmin[spc][n]*R->DR[spc][n])-0.0;

    for(k=0;k<R->NL-1;k++)
      for(m=0;m<R->NU[k+1];m++)
      {
	for( n=0, e[i][k+1][m]=R->B[spc][k][m]; n<R->NU[k]; n++)
	  e[i][k+1][m] += e[i][k][n]*R->W[spc][k][m][n];
	e[i][k+1][m] =  GP(R->GT[k+1],e[i][k+1][m]);
      }
    E += e[i][R->NL-1][0];
  }
  
  c = 2.0*(E-L->E)/(double)(L->N*L->N)*R->WE*L->W;

  for(i=0;i<L->N;i++)
  {
    spc = L->ATMN[i];
    k = R->NL-2;
    if(R->MIX==0)
    {
      Bp[spc][k][0] += c;
      for(m=0;m<R->NU[k];m++)
	Wp[spc][k][0][m] += c*e[i][k][m];
      k--;
      
      for(m=0;m<R->NU[k+1];m++)
      {
	tm = c*R->W[spc][k+1][0][m]*GPp(R->GT[k+1],e[i][k+1][m]);
	Bp[spc][k][m] += tm;
	for(n=0;n<R->NU[k];n++)
	  Wp[spc][k][m][n] += tm*e[i][k][n];
      }
    }
    else
      k--;

    if(k)
    {
      k--;      
      for(m=0;m<R->NU[k+2];m++)
      {
	tm = c*R->W[spc][k+2][0][m]*GPp(R->GT[k+2],e[i][k+2][m]);
	for(n=0;n<R->NU[k+1];n++)
	{
	  tn = tm*R->W[spc][k+1][m][n]*GPp(R->GT[k+1],e[i][k+1][n]);
	  Bp[spc][k][n] += tn;
          for(l=R->O;l<R->NU[k];l++)
	    Wp[spc][k][n][l] += tn*e[i][k][l];
	}
      }
    }
  }  
  return E;
}
//================================================================
//
//================================================================
double TOT_ERR(ANN *R, LNK *L)
{
  int n;
  for(n=0,R->RE=0.0;n<R->N;n++)
    R->RE += pow( (L[n].E-ENE_ANN(R,&L[n]))/(double)L[n].N,2.0 );
  R->RE  = (R->RE/(double)R->N);
  return R->RE;  
}
//================================================================
double ENE_ANN(ANN *R, LNK *L)
{
  int i,k,n,m,spc;
  double E,a[500],b[500]; ///max number of vectors hardcoded!!!!

  E = 0.0;
  #pragma omp parallel private(spc,n,k,m,a,b) num_threads(R->NP)
  { 
    #pragma omp for reduction (+:E) schedule(dynamic,R->NB)
    for(i=0;i<L->N;i++)
    {
      spc = L->ATMN[i];
      for(n=0;n<R->NU[0];n++)
        a[n] = (L->Cn[i][n]-R->Rmin[spc][n]*R->DR[spc][n])-0.0;

      for(k=0;k<R->NL-1;k++)
      {
        for(m=0;m<R->NU[k+1];m++)
        {
	  for(n=0,b[m]=R->B[spc][k][m];n<R->NU[k];n++)
	    b[m] += a[n]*R->W[spc][k][m][n];
	  b[m] = GP(R->GT[k+1],b[m]);
        }
	
        for(m=0;m<R->NU[k+1];m++)
	  a[m] = b[m];
      }
      L->EA[i] = b[0];
      E += b[0];
    }
  }
  return E;
}
//======================================================
double FRC_ANN(ANN *R, LNK *L)
{
  int i,j,k,n,m,l,q,nij,sij,spc,nth;
  double tm,tn,**s,e;

  e = 0.0;
  s = make_d2D(R->NP,6);
  for(nth=0;nth<R->NP;nth++)
    for(q=0;q<6;q++)
      s[nth][q] = 0.0;
  
  #pragma omp parallel private(spc,n,k,m) num_threads(R->NP)
  {
    #pragma omp for reduction (+:e) schedule(dynamic,R->NB)
    for(i=0;i<L->N;i++)
    {
      spc = L->ATMN[i];
      for(n=0;n<R->NU[0];n++)
        R->e[spc][i][0][n] = (L->Cn[i][n]-R->Rmin[spc][n]*R->DR[spc][n])-0.0;
      for(k=0;k<R->NL-1;k++)
        for(m=0;m<R->NU[k+1];m++)
        {
	      for( n=0, R->e[spc][i][k+1][m]=R->B[spc][k][m]; n<R->NU[k]; n++)
	        R->e[spc][i][k+1][m] += R->e[spc][i][k][n]*R->W[spc][k][m][n];
	      R->e[spc][i][k+1][m] =  GP(R->GT[k+1],R->e[spc][i][k+1][m]);
	      R->d[spc][i][k+1][m] = GPp(R->GT[k+1],R->e[spc][i][k+1][m]);
        }
      L->EA[i] = R->e[spc][i][R->NL-1][0];
      e += R->e[spc][i][R->NL-1][0];
    }
  }

  #pragma omp parallel private(nth,spc,q,m,tm,n,tn,l,j,nij,sij) num_threads(R->NP)
  {
    nth = omp_get_thread_num();
    #pragma omp for schedule(dynamic,R->NB)
    for(i=0;i<L->N;i++)
    {
      spc = L->ATMN[i];
      for(q=0;q<3;q++)
        L->f[i][q] = 0.0;
      // for 3-layer MLP
      // =======  dGni/dxi  ========
      for(m=0;m<R->NU[2];m++)
      {
        tm = R->W[spc][2][0][m]*R->d[spc][i][2][m];
        for(n=0;n<R->NU[1];n++)
        {
	      tn = tm*R->W[spc][1][m][n]*R->d[spc][i][1][n];
	      for(l=0;l<R->NU[0];l++)
	      {
	        if(L->MRK[i]==1)
	          for(q=0;q<3;q++)
	            L->f[i][q] -= L->Fn[i][L->DNn[i]][l][q]*tn*R->W[spc][0][n][l];
	        for(q=0;q<6;q++)
	          s[nth][q] -= L->Sn[i][l][q]*tn*R->W[spc][0][n][l];
	      }
        }
      }      
      // =======  dGnj/dxi  ========
      if(L->MRK[i]==1)
        for(j=0;j<L->DNn[i];j++)
        {
	      nij = L->DNi[i][j];
	      if(j==L->DNn[i]) sij = spc; else sij = L->DNs[i][j];

	      for(m=0;m<R->NU[2];m++)
	      {
	        tm = R->W[sij][2][0][m]*R->d[L->DNs[i][j]][nij][2][m];
	        for(n=0;n<R->NU[1];n++)
  	        {
	          tn = tm*R->W[sij][1][m][n]*R->d[sij][nij][1][n];
	          for(l=0;l<R->NU[0];l++)
	          for(q=0;q<3;q++)
		        L->f[i][q] -= L->Fn[i][j][l][q]*tn*R->W[sij][0][n][l];;
	        }      
	      }
        }
    }  
  }
  L->e = e;
  for(q=0;q<6;q++)  
    for(nth=0,L->s[q]=0.0;nth<R->NP;nth++)
      L->s[q] += s[nth][q];

  free_d2D(s,R->NP);
  return L->e;
}
//================================================================
void INIT_MLP(ANN *R)
{
  int  k,n,m,i,spc;
  FILE *in;
  char s[400];

  sprintf(s,"mkdir -p %s/",R->otpt);
  system(s);

  sprintf(s,"%s/model",R->otpt);

  in = fopen(s,"r");
  //=====  read the NN if 'model' file exists  =====
  if( in!= NULL )
  {
    fclose(in);
    READ_ANN(R);
  }
  //=====  randomly initialize the NN otherwise  =====
  else
  {
    for(spc=0;spc<R->NSPC;spc++)
    {
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  for(n=0;n<R->NU[k];n++,i++)
	    R->W[spc][k][m][n] = 2.0*(Random()-0.5)/(double)R->NU[k];
      
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++,i++)
	  R->B[spc][k][m] = 2.0*(Random()-0.5)/(double)R->NU[k];
    }
  }
  R->O = 0;
  //=====  for stratified training load sub-species NNs  =====
  if(R->MIX==1) 
  {
    if(R->NSPC > 1 && R->NSPC < 4) 
    {
      READ_STRAT(R);
      return;
    }
    fprintf(stderr,"Stratified training is currently available only for binaries and ternaries!\n");
    exit(1);    
  }
}
//================================================================
// copy scp MLP weights into a 1D vector 
//================================================================
int W4V(ANN *R, double *V, int *SPC)
{
  int i=0,k,n,m,spc,ii,jj,N1,N2;

  N1 =  8;
  N2 = 22;
  for(spc=0;spc<R->NSPC;spc++)
    if( SPC[spc]==1 )
    {
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  V[i++] = R->B[spc][k][m];  
      for(k=0;k<1;k++)
        for(m=0;m<R->NU[1];m++)
	{
	  for(ii=0;ii<N1;ii++)
	    V[i++] = R->W[spc][k][m][ii*2+spc];
          for(jj=0;jj<N2;jj++)
            V[i++] = R->W[spc][k][m][ii*2+jj*3+spc];
	}
      for(k=1;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  for(n=0;n<R->NU[k];n++)
	    V[i++] = R->W[spc][k][m][n];
    }
  return i;
}
//================================================================
//  copy all MLP weights into a 1D vector of length R->NW
//================================================================
void W2V(ANN *R, double *V)
{
  int i=0,k,n,m,spc;


  for(spc=0;spc<R->NSPC;spc++)
    {
      // biases first
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  V[i++] = R->B[spc][k][m];
      
      // rest of the neurons next
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  for(n=0;n<R->NU[k];n++)
	    V[i++] = R->W[spc][k][m][n];
    }
}
//================================================================
//  copy weights from a 1D vector into the MLP's weights (W & B)
//================================================================
void V2W(ANN *R, double *V)
{
  int i=0,k,n,m,spc;

  for(spc=0;spc<R->NSPC;spc++)
    {
      // biases first
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  R->B[spc][k][m] = V[i++];
      
      // rest of the neurons next
      for(k=0;k<R->NL-1;k++)
	for(m=0;m<R->NU[k+1];m++)
	  for(n=0;n<R->NU[k];n++)
	    R->W[spc][k][m][n] = V[i++];
    }
}
//================================================================
//
//================================================================
void TRAN_MLP(ANN *R, Cell *C, LNK *L)
{
  int n,i;

  FILE *out;
  char s[400];

  OOO = (PAR *)malloc(R->NP*sizeof(PAR));
  for(n=0;n<R->NP;n++)
  {
    OOO[n].e  = make_d3D(R->A,R->NL,R->N0);
    OOO[n].Wp = make_d4D(R->NSPC,R->NL,R->N1,R->N0);
    OOO[n].Bp = make_d3D(R->NSPC+1,R->NL-1,R->N1);
    if(R->EFS==1||R->EFS==3)
    {
      OOO[n].d  = make_d3D(R->A,  R->NL,R->N0);
      OOO[n].c  = make_d3D(R->A,  R->NL,R->N0);
      OOO[n].xn = make_d3D(R->NU[1],R->DNm+1,3);
      OOO[n].xm = make_d3D(R->NU[2],R->DNm+1,3);
    }
  }

  sprintf(s,"%s/err-out.dat",R->otpt);

  //=====  Train the MLP here =====
  for(i=0;i<1;i++)
  {
    Random();
    
    printf("Total number of parameters: %3d\n",R->NW);
    out=fopen(s,"w");
    FPRNT_HEAD(out);
    fprintf(out,"|                            Model training                           |\n");
    fprintf(out,"=======================================================================\n\n");
    fprintf(out,"Total number of parameters: %3d\n",R->NW);
    fclose(out);
    
    if(R->MINT>=0 && R->MINT<4)
      MLP_MIN(R,L);
    else
    {
      fprintf(stderr,"Error: Enter valid value for MINT (0-3)\n");
      exit(1);
    }
     
  }
  return;
}
//================================================================
