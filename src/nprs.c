#include "nprs.h"

//==================================================================
//   0 in output for improper 'basis' file and 1 for good file
//==================================================================
void GNcheck(char *s, int coor)
{
  FILE   *in;
  int    n,index=0;
  char   buf[400];
  double a,r;
  
  if( (in=fopen(s,"r")) == 0 )
  {
    fprintf(stderr,"ERROR opening %s\n",s);
    exit(1);
  }

  while(fgets(buf,200,in))
    if(strncmp(buf,"  B2A",5)==0)
      break;

  fgets(buf,200,in);
  sscanf(buf,"%lf %lf\n",&a,&r);
  if(0)
  if( fabs(a-0.529177249) > 1e-8 )
  {
    fclose(in);
    fprintf(stderr,"Error in 'basis' file format(1)!\n"); 
    exit(1);
  }
  fgets(buf,200,in);
  fgets(buf,200,in);
  fgets(buf,200,in);
  if( strlen(buf) != 1 )  
  {
    fclose(in);
    fprintf(stderr,"Error in 'basis' file format(2)!\n"); 
    exit(1);
  }
  for(n=0;n<coor;n++)
  {
    fgets(buf,200,in);
    if( strlen(buf) == 1 ) 
      index++;
  }
  
  if( index > 0 ) 
  {
    fclose(in);
    fprintf(stderr,"Error: incorrect number of symmetry functions in 'basis' file!\n"); 
    exit(1);
  }

  fgets(buf,200,in);
  if( strlen(buf) != 1 ) 
  {
    fclose(in);
    fprintf(stderr,"Error in 'basis' file format(3) %s!\n",buf); 
    exit(1);
  }
  fgets(buf,200,in);
  if( strlen(buf) == 1 ) 
  {
    fclose(in);
    fprintf(stderr,"Error in 'basis' file format(4)!\n"); 
    exit(1);
  }
  
  fclose(in);
}
//==================================================================
//
//==================================================================
void Build_PRS(PRS *P, PRS *W, int J)
{
  FILE   *in;
  int    i,p,l,n,m,w0,w1,w00,w01,w11;
  char   buf[400],s[400],basis[400];
  double a,r;
  int    c=0,d,ver;
  
  P->Pl   = make_d2D(P->LM+1,P->LM+1);
  P->GW   = make_d2D(P->GM,P->GM);
  
  P->GF   = make_i2D(P->D,9);
  P->GP   = make_d2D(7,10);               // 6 PB types of parameters x 10 
                                                     
  if(J==0)
  {
    sprintf(basis,"model");
    ver=check_ver(basis);
    if( ver < 2400 )
    {
      printf("Error: the model file is in old format; convert the model first!\n");
      exit(1);
    }
  }
  else
    sprintf(basis,"basis");

  if(P->DSCR==1)
  {
    in = fopen("INI/w.dat","r");
    for(n=0;n<P->GN-1;n++)
    {
      for(m=0;m<=n;m++)
	fgets(buf,200,in);
      fgets(buf,200,in);
    }
    for(n=0;n<P->GN;n++)
    {
      fgets(buf,200,in);
      sscanf(buf,"%lf %lf %lf %lf %lf %lf",&P->GW[n][0],&P->GW[n][1],&P->GW[n][2],&P->GW[n][3],&P->GW[n][4],&P->GW[n][5]);
    }
    for(n=0;n<P->GN;n++)
      for(i=0;i<P->GN;i++)
	P->GW[n][i] *= sqrt((2.0*(double)(i+1.0)+5.0)/pow(P->RC,2.0*(double)(i+1.0)+5.0));
    fclose(in);
    in = fopen("INI/pl.dat","r");
    for(l=0;l<=P->LM;l++)
    {
      fscanf(in,"%s",s);
      for(p=l%2;p<=l;p+=2)
      {
	fscanf(in,"%lf\n",&P->Pl[l][p]);
	P->Pl[l][p] *= (double)(2*l+1)/(4.0*Pi);
      }
      fscanf(in,"\n");
    }
    fclose(in);
  }
  if(P->DSCR==2)
  {
    GNcheck(basis,P->NSYM);
    for(w0=0;w0<P->NSPC;w0++)
      for(w1=0;w1<P->NSPC;w1++)
      {
	w01         = w0*w0+w1*w1;
	W[w01].GF   = make_i2D(P->D,9);
	W[w01].GP   = make_d2D(7,10);
	W[w01].GT   = make_i1D(7);
	W[w01].NSPC = P->NSPC;
	W[w01].NSYM = P->NSYM;
	W[w01].D    = P->D;
	W[w01].IO   = P->IO;  
        W[w01].EFS  = P->EFS; 
        W[w01].FMRK = P->FMRK;
        // W[w01].SPCZ = make_i1D(P->NSPC);
	for(i=0;i<P->NSPC;i++) 
	  W[w01].SPCZ[i] = P->SPCZ[i];
      }

    in = fopen(basis,"r");
    while(fgets(s,200,in))
      if(strncmp(s,"  B2A",5)==0)
	break;

    for(w0=0;w0<P->NSPC;w0++)
    {
      w00 = 2*w0*w0;
      fgets(buf,200,in);
      sscanf(buf,"%lf %lf\n",&a,&r);
      fgets(buf,200,in);
      fgets(buf,200,in);
      fgets(buf,200,in);
      for(n=0;n<W[w00].NSYM;n++)
      {
	fgets(buf,200,in);
	sscanf(buf,"%d %d %d %d %d %d %d %d %d",&i,
	       &W[w00].GF[n][0],&W[w00].GF[n][1],&W[w00].GF[n][2],&W[w00].GF[n][3],
	       &W[w00].GF[n][4],&W[w00].GF[n][5],&W[w00].GF[n][6],&W[w00].GF[n][7]);
      }
    
      fgets(buf,200,in);
      for(n=0;n<7;n++)
      {
	fgets(buf,200,in);
	sscanf(buf,"%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",s,&W[w00].GT[n],
	       &W[w00].GP[n][0],&W[w00].GP[n][1],&W[w00].GP[n][2],&W[w00].GP[n][3],&W[w00].GP[n][4],
	       &W[w00].GP[n][5],&W[w00].GP[n][6],&W[w00].GP[n][7],&W[w00].GP[n][8],&W[w00].GP[n][9]);
      }
      for(n=0;n<2;n++)
	for(i=0;i<10;i++)
	  W[w00].GP[n][i] *= a*r;
    
      for(n=2;n<4;n++)
	for(i=0;i<10;i++)
	  W[w00].GP[n][i] /= a*a*r*r;

      for(n=0;n<P->D;n++)
        for(i=0;i<9;i++)
          P->GF[n][i] = W[0].GF[n][i];

      for(n=0;n<7;n++)
        for(i=0;i<10;i++)
	  P->GP[n][i] = W[0].GP[n][i];

      P->RC = W[w00].RC = W[w00].GP[0][0];
      for(n=0;n<W[w00].NSYM;n++)
	if(W[w00].GP[0][W[w00].GF[n][1]]>W[w00].RC)
	  W[w00].RC = W[w00].GP[0][W[w00].GF[n][1]];
      //counting actual number of G2 and G4 functions from 'basis' file 
      W[w00].NG2=0;
      W[w00].NG4=0;
      for(n=0;n<W[w00].NSYM;n++)
      {
	if(W[w00].GF[n][0]==2) W[w00].NG2++;
	if(W[w00].GF[n][0]==4) W[w00].NG4++;
      }
    
      //counting components number from 'basis' file and matching that with NCMP in setup file             
      for(c=0,n=W[w00].NSPC;n>0;c+=n,n--);

      d = (W[w00].NG2*W[w00].NSPC)+(W[w00].NG4*c);
      if(d!=W[w00].D) 
      { 
	fprintf(stderr,"Error: NCMP is not consistent with 'basis' file, suggested is: %d!\n",d); 
	exit(1); 
      }
      fgets(buf,200,in);
      //=====  parsers could be species-specific; if only one provided used it for all species  =====
      if( w0<P->NSPC-1 && fgets(buf,200,in)==0 )
      {
	fclose(in);
	in = fopen(basis,"r");
	while(fgets(s,200,in))
	  if(strncmp(s,"  B2A",5)==0)
	    break;
      }
    }
    fclose(in);
    for(w0=0;w0<P->NSPC;w0++)
      for(w1=0;w1<P->NSPC;w1++)
	if(w0!=w1)
	{
	  w00 = 2*w0*w0;
	  w01 = w0*w0+w1*w1;
	  w11 = 2*w1*w1;
          W[w01].RC = 0.5*(W[w00].RC+W[w11].RC);
	  // for parameters linear in distance
 	  for(n=0;n<2;n++)
	    for(i=0;i<10;i++)
	      W[w01].GP[n][i] = 0.5*(W[w00].GP[n][i]+W[w11].GP[n][i]);
	  // for parameters inverse square in distance (n1 and n2)
	  for(n=2;n<4;n++)
	    for(i=0;i<10;i++)
	      W[w01].GP[n][i] = pow( 0.5*(pow(W[w00].GP[n][i],-0.5)+ pow(W[w11].GP[n][i],-0.5)),-2.0 );
          for(n=0;n<P->D;n++)
            for(i=0;i<9;i++)
              W[w01].GF[n][i] = W[0].GF[n][i];
	}
    P->NG2 = W[0].NG2;
    P->NG4 = W[0].NG4;
  }
  return;
}
//==================================================================                                              
double dfc(double Rc, double r)
{
  if(r>Rc)
    return 0.0;
  return -0.5*Pi/Rc*sin(Pi*r/Rc);
}
//==================================================================                                                              
double fc(double Rc, double r)
{
  if(r>Rc)
    return 0.0;
  return 0.5*(cos(Pi*r/Rc)+1.0);
}
//==================================================================
double dGNR(PRS *P, int n, double Rc, double r, double fcij, double dfcij)
{
  return exp(-P->GP[2][P->GF[n][3]]*pow(r-P->GP[1][P->GF[n][2]],2.0))*(-2.0*P->GP[2][P->GF[n][3]]*(r-P->GP[1][P->GF[n][2]])*fcij+dfcij);
}
//==================================================================
double dgnr(PRS *P, int n, double Rc, double r)
{
  if(P->GF[n][0]==1)
    return dfc(Rc,r);
  if(P->GF[n][0]==2)
    return exp(-P->GP[2][P->GF[n][3]]*pow(r-P->GP[1][P->GF[n][2]],2.0))*(-2.0*P->GP[2][P->GF[n][3]]*(r-P->GP[1][P->GF[n][2]])*fc(Rc,r)+dfc(Rc,r));
  if(P->GF[n][0]==3)
    return dfc(Rc,r)*cos(P->GP[3][P->GF[n][4]]*r)-fc(Rc,r)*P->GP[3][P->GF[n][4]]*sin(P->GP[3][P->GF[n][4]]*r);
  return 0.0;
}
//==================================================================                                                              
double GNR(PRS *P, int n, double Rc, double r)
{
//  printf("%3d % lf % lf\n",n,P->GP[1][P->GF[n][2]],P->GP[2][P->GF[n][3]]);
  if(P->GF[n][0]==1)
    return fc(Rc,r);
  if(P->GF[n][0]==2)
    return fc(Rc,r)*exp(-P->GP[2][P->GF[n][3]]*pow(r-P->GP[1][P->GF[n][2]],2.0));
  if(P->GF[n][0]==3)
    return fc(Rc,r)*cos(P->GP[3][P->GF[n][4]]*r);
  if(P->GF[n][0]==5||P->GF[n][0]==6)
    return fc(Rc,r)*pow( P->GP[1][P->GF[n][2]]/r, P->GP[2][P->GF[n][3]] );      
  return 0.0;
}
//==================================================================                                                              
double GNA(PRS *P, int n, double R0, double n0, double r0, double R1, double n1, double r1, double R2, double n2, double r2, double a)
{
  return fc(R0,r0)*fc(R1,r1)*fc(R2,r2) * pow(1.0+P->GP[6][P->GF[n][7]]*a,P->GP[5][P->GF[n][6]]) * exp(-(r0*r0*n0+r1*r1*n1+r2*r2*n2));
}
//==================================================================
void LNK_IN(LNK *L, int o, char *path)
{
  int i,j,n,tmp;
  int N,D,EFS,DNn;
  char buf[400];
  FILE *in;
  
  sprintf(buf,"%s/e%06d",path,o);
  if( (in=fopen(buf,"r")) == 0 )
  {
    fprintf(stderr,"ERROR opening %s\n",buf);
    exit(1);
  }
  L->p = 0.0;
  fscanf(in,"%d %d %d %d %d %lf %lf %lf\n",&EFS,&D,&N,&L->NF,&DNn,&L->E,&L->p,&L->DE);
  L->N = N;

  Build_LNK(L,N,DNn,D,EFS);
  
  if(EFS==1||EFS==3)
  {
    for(i=0;i<L->NF;i++)
      fscanf(in,"%3d ",&L->Fi[i]);
    fscanf(in,"\n");
    for(i=0;i<L->N;i++)
      L->MRK[i] = 0;
    for(i=0;i<L->NF;i++)
      L->MRK[L->Fi[i]] = 1;
  }

  if(EFS==2||EFS==3)
    fscanf(in,"%lf %lf %lf %lf %lf %lf\n",&L->S[0],&L->S[1],&L->S[2],&L->S[3],&L->S[4],&L->S[5]);

  fscanf(in,"%s\n",L->path);
  //=====  for Energy  =====
  for(i=0;i<N;i++)
    {
      fscanf(in,"%d\n",&L->ATMN[i]);
      for(n=0;n<D;n++)
	fscanf(in,"%lf\n",&L->Cn[i][n]);
    }

  //=====  for Forces  =====
  if(EFS==1||EFS==3)
    for(i=0;i<N;i++)
      if(L->MRK[i]==1) 
	{
	  fscanf(in,"%d %d",&tmp,&L->DNn[i]);
	  fscanf(in,"%lf %lf %lf",&L->F[i][0],&L->F[i][1],&L->F[i][2]);
	  for(n=0;n<D;n++)
	    fscanf(in,"%lf %lf %lf\n",&L->Fn[i][L->DNn[i]][n][0],&L->Fn[i][L->DNn[i]][n][1],&L->Fn[i][L->DNn[i]][n][2]);
	  for(j=0;j<L->DNn[i];j++)
	    {
	      fscanf(in,"%4d\n",&L->DNi[i][j]);
	      L->DNs[i][j] = L->ATMN[L->DNi[i][j]];
	      for(n=0;n<D;n++)
		fscanf(in,"%lf %lf %lf\n",&L->Fn[i][j][n][0],&L->Fn[i][j][n][1],&L->Fn[i][j][n][2]);
	    }      
	  L->DNs[i][L->DNn[i]] = L->ATMN[i];
	  L->DNi[i][L->DNn[i]] = i;
	  n = 0;
	}
  
  //=====  for Stresses  =====
  if(EFS==2||EFS==3)
    for(i=0;i<L->N;i++)
      for(n=0;n<D;n++)
	fscanf(in,"%lf %lf %lf %lf %lf %lf\n",&L->Sn[i][n][0],&L->Sn[i][n][1],&L->Sn[i][n][2],&L->Sn[i][n][3],&L->Sn[i][n][4],&L->Sn[i][n][5]);

  fclose(in);
}
//==================================================================
void LNK_OUT(LNK *L, int o, char *path, int EFS, int D)
{
  int i,j,n,q,DNn;
  char buf[400];
  FILE *out;
  
  sprintf(buf,"%s/e%06d",path,o);
  out = fopen(buf,"w");
  
  for(i=DNn=0;i<L->N;i++)
    if(L->MRK[i]&&DNn<L->DNn[i])
      DNn = L->DNn[i];
  
  fprintf(out,"%d %3d %3d %3d %3d   % lf  % lf  % lf\n",EFS,D,L->N,L->NF,DNn,L->E,L->p,L->DE);
  if(EFS==1||EFS==3)
    {
      for(i=0;i<L->N;i++)
	if(L->MRK[i])
	  fprintf(out,"%3d ",i);
      fprintf(out,"\n");
    }
  if(EFS==2||EFS==3)
  {
    for(q=0;q<6;q++)
      fprintf(out,"% 28.16lf ",L->S[q]);
    fprintf(out,"\n");
  }
  fprintf(out,"%s\n",L->path);

  //=====  for Energy  =====
  for(i=0;i<L->N;i++)
  {
    fprintf(out,"%d\n",L->ATMN[i]);
    for(n=0;n<D;n++)
      fprintf(out,"% 28.16lf\n",L->Cn[i][n]);
  }

  //=====  for Forces  =====
  if(EFS==1||EFS==3)
  for(i=0;i<L->N;i++)
    if(L->MRK[i]==1)
    {
      fprintf(out,"%4d %4d ",i,L->DNn[i]);
      for(q=0;q<3;q++)
	fprintf(out,"% lf ",L->F[i][q]);
      fprintf(out,"\n");
      for(n=0;n<D;n++)
	fprintf(out,"% 28.16lf % 28.16lf % 28.16lf\n",L->Fn[i][L->DNn[i]][n][0],L->Fn[i][L->DNn[i]][n][1],L->Fn[i][L->DNn[i]][n][2]);
      for(j=0;j<L->DNn[i];j++)
      {
	fprintf(out,"%4d\n",L->DNi[i][j]);
	for(n=0;n<D;n++)
	  fprintf(out,"% 28.16lf % 28.16lf % 28.16lf\n",L->Fn[i][j][n][0],L->Fn[i][j][n][1],L->Fn[i][j][n][2]);
      }    
    }
  //=====  for Stresses  =====
  if(EFS==2||EFS==3)
    for(i=0;i<L->N;i++)
      for(n=0;n<D;n++,fprintf(out,"\n"))
	for(q=0;q<6;q++)
	  fprintf(out,"% 28.16lf ",L->Sn[i][n][q]);
  fclose(out);
}

//==================================================================                                                              
void PRS_BP(PRS *P, PRS *W, Cell *C, LNK *L, int o, char *path)
{
  int    i,j,k,n,q,*I,ii,jj;
  double t,A0,A1,A2,A3,A4,A5,A6,A7,B1,B2,B3,T0,F0,F1,F2,R01,R02,R12,n01,n02,n12;
  double fcijk,lzijk[10][10],a2ijk[10][10],a3ijk[10][10],n2ijk[10],a52,a53,a62,a63,ndxjk,fcij,dfcij,powz[10];
  double dxij[3],dxik[3],dxji[3],dxjk[3],dqij[3],dqik[3];
  int    l,z,n2,n4,n6,n7; 

  int    ind2[10],ind4[10];
  int    m,s0,s1,s2,s01,s02,s12,ds1,ds2,ds,g2,g4;

  L->N   = C->N;
  C->Rc  = C->rc = P->RC;

  for(s0=0;s0<P->NSPC;s0++)
    for(n=0;n<P->NSYM;n++)
      if(C->rc > W[2*s0*s0].GP[0][W[2*s0*s0].GF[n][1]])
	C->Rc = C->rc = W[2*s0*s0].GP[0][W[0].GF[n][1]];
  for(z=0;z<W[0].GT[5];z++)
    powz[z] = pow(2.0,1.0-P->GP[5][z]);
  LIST(C,0);  

  //===== symmetry functions by species  =====
  //-----ij-ijk-------------------------------
  //    |AA|AAA|
  //------------------------------------------
  //    |AA|AAA|AB|AAB|ABB|
  //    |BB|BBB|BA|BBA|BAA|
  //------------------------------------------
  //    |AA|AAA|AB|AAB|ABB|AC|AAC|ACC|ABC|
  //    |BB|BBB|BC|BBC|BCC|BA|BBA|BAA|BCA|
  //    |CC|CCC|CA|CCA|CAA|CB|CCB|CBB|CAB|
  //==========================================

  ind2[0] = 0;
  ind2[1] = 1*P->NG2 + 1*P->NG4;
  ind2[2] = 2*P->NG2 + 3*P->NG4;

  ind4[0] = 1*P->NG2;
  ind4[1] = 2*P->NG2 + 1*P->NG4;
  ind4[2] = 2*P->NG2 + 2*P->NG4;
  ind4[4] = 3*P->NG2 + 3*P->NG4;
  ind4[8] = 3*P->NG2 + 4*P->NG4;
  ind4[5] = 3*P->NG2 + 5*P->NG4;

  for(i=0;i<C->N;i++) 
    for(n=0;n<P->D;n++)   
      L->Cn[i][n] = 0.0;
  // For EFS
  #pragma omp parallel private(ds1,ds2,ds,jj,k,n,m,g2,g4,s0,s1,s01,R01,s2,s02,s12,R02,R12,fcijk,n2,n01,n02,n12,n2ijk,lzijk,z,l) num_threads(C->NP)
  {
    #pragma omp for schedule(dynamic,C->NB)
    for(i=0;i<C->N;i++)
    {
      s0 = C->ATMN[i];
      for(jj=0;jj<C->Nn[i];jj++)
      {
        n   = 0;
        s1  = C->ATMN[C->Ni[i][jj]];
        s01 = s0*s0+s1*s1;
        R01 = W[s01].GP[0][P->GF[n][1]];
	ds1 = (s1-s0+P->NSPC)%P->NSPC;
        
        if(C->NDX[i][jj]<R01)
        {
	  for(n=0,g2=0;n<P->NSYM;n++)
	    if(P->GF[n][0]==2)
	    {
	      m   = g2 + ind2[ds1];
	      L->Cn[i][m] += GNR(&W[s01],n,R01,C->NDX[i][jj]);
	      g2++;
	    }
	  // for reduced radius for triplets	
	  if( g2<P->D && C->NDX[i][jj]<W[s01].GP[0][P->GF[g2][1]] )
	  for(k=jj+1;k<C->Nn[i];k++)
          {
	    s2  = C->ATMN[C->Ni[i][k]];
	    s02 = s0*s0+s2*s2;
	    s12 = s1*s1+s2*s2;
	    R02 = W[s02].GP[0][P->GF[g2][1]];
	    R12 = W[s12].GP[0][P->GF[g2][1]];
	    ds2 = (s2-s0+P->NSPC)%P->NSPC;
	    if( C->NDX[i][k]<R02 && C->ndx[i][jj][k]<R12 )
	    {
	      // assuming that all G4 Rc are the same
	      fcijk = 2.0*fc(R01,C->NDX[i][jj])*fc(R02,C->NDX[i][k])*fc(R12,C->ndx[i][jj][k]);
	      for(n2=0;n2<W[0].GT[3];n2++)
	      {
	        n01 = W[s01].GP[3][n2];
	        n02 = W[s02].GP[3][n2];
	        n12 = W[s12].GP[3][n2];
	        n2ijk[n2] = fcijk*exp( -C->NDX[i][jj]*C->NDX[i][jj]*n01 - C->NDX[i][ k]*C->NDX[i][ k]*n02 - C->ndx[i][jj][k]*C->ndx[i][jj][k]*n12 );
	      }

	      for(z=0;z<W[0].GT[5];z++)
	      for(l=0;l<W[0].GT[6];l++)
	        lzijk[z][l] =  powz[z]*pow(1.0+P->GP[6][l]*C->cos[i][jj][k],P->GP[5][z]);
	      
	      ds = ds1*ds1 + ds2*ds2;
	      for(n=0,g4=0;n<P->NSYM;n++)
	        if(P->GF[n][0]==4)
	        {
		  m = g4 + ind4[ds];
		  L->Cn[i][m] += lzijk[W[0].GF[n][6]][W[0].GF[n][7]] * n2ijk[W[0].GF[n][4]];
		  g4++;
	        }
	    }
	  }
        }
      }
    }
  }
  // For FS
  if(P->EFS>0) 
  {
  #pragma omp parallel private(j,n,q,I,jj,ii,fcij,dfcij,dxji,g2,t,s0,m,k,A0,dxjk,n2,n2ijk,z,l,a2ijk,a3ijk,s1,s2,ds1,ds2,ds,a52,a53,a63,ndxjk,g4,n4,n6,n7,A2,A3,A4,A5,A6,A7,dxij,dqij,dxik,dqik,T0,F0,F1,F2,a62,B1,A1,B2,B3) num_threads(C->NP)
    {
      I = make_i1D(C->N);
      #pragma omp for schedule(dynamic,C->NB)
      for(i=0;i<C->N;i++)
      {      
  	s0 = C->ATMN[i];
        for(j=0;j<=L->DNn[i];j++)
	  for(n=0;n<P->D;n++)
	    for(q=0;q<3;q++)
	      L->Fn[i][j][n][q] = 0.0;
        // this dG/dxj needs to be done before dG/dxi so that dG/dxi are stored in Fn[i][L->DNn[i]]...
        if(P->EFS!=2 && L->MRK[i]==1)
        {
	  L->DNn[i]=0;	
	  for(j=0;j<C->N;j++)
	    I[j] = 0;
	  for(jj=0;jj<C->Nn[i];jj++)
	  {
	    j  = C->Ni[i][jj];
  	    s1 = C->ATMN[j];
	    if(I[j]==0)
	    {
	      I[j] = 1;
	      for(ii=0;ii<C->Nn[j];ii++)
	        if(C->Ni[j][ii]==i)
	        {
		  fcij  =  fc(P->GP[0][P->GF[0][1]],C->NDX[j][ii]);
		  dfcij = dfc(P->GP[0][P->GF[0][1]],C->NDX[j][ii]);
                  for(q=0;q<3;q++)
                    dxji[q] = DX(C,j,ii,q);
		  for(n=0,g2=0;n<P->NSYM;n++)
		    if(P->GF[n][0]==2)
		    {
		      t = dGNR(P,n,P->GP[0][P->GF[n][1]],C->NDX[j][ii],fcij,dfcij)/C->NDX[j][ii];
                      L->DNs[i][L->DNn[i]] = C->ATMN[j];
		      ds1 = (s0-s1+P->NSPC)%P->NSPC;
		      m  = g2 + ind2[ds1];
		      for(q=0;q<3;q++)
		        L->Fn[i][L->DNn[i]][m][q] += t*dxji[q]; // the sign is correct!
		      // for the pair function in stress calculation only the part on atom i is needed
		      g2++;
		    }
		  if( g2<P->D && C->NDX[j][ii]<P->GP[0][P->GF[g2][1]] )
		  for(k=0;k<C->Nn[j];k++)
		    if( k!=ii && C->ndx[j][ii][k]<P->GP[0][P->GF[g2][1]] && C->NDX[j][k]<P->GP[0][P->GF[g2][1]] )
		    {
		      A0 = fc(P->GP[0][P->GF[g2][1]],C->NDX[j][ii])*fc(P->GP[0][P->GF[g2][1]],C->ndx[j][ii][k]);
		      for(q=0;q<3;q++)
		        dxjk[q] = DX(C,j,k,q);
		      // works only for the same Rc for all species
		      for(n2=0;n2<W[0].GT[3];n2++)
		        n2ijk[n2] = 2.0*exp(-P->GP[3][n2]*(pow(C->NDX[j][ii],2.0)+pow(C->NDX[j][k],2.0)+pow(C->ndx[j][ii][k],2.0)))*fc(P->GP[0][P->GF[g2][1]],C->NDX[j][k]);
		      for(z=0;z<W[0].GT[5];z++)
		      for(l=0;l<W[0].GT[6];l++)
		      {
		        a2ijk[z][l] = powz[z]* A0*P->GP[6][l]*P->GP[5][z]*pow(1.0+P->GP[6][l]*C->cos[j][ii][k],P->GP[5][z]-1.0)/C->NDX[j][ii];
		        a3ijk[z][l] = powz[z]*pow(1.0+P->GP[6][l]*C->cos[j][ii][k],P->GP[5][z]);
		      }
		      s2    = C->ATMN[C->Ni[j][k]];
                      ds1   = (s0-s1+P->NSPC)%P->NSPC;
                      ds2   = (s2-s1+P->NSPC)%P->NSPC;
                      ds    = ds1*ds1 + ds2*ds2;

		      a52   = -C->cos[j][ii][k]/C->NDX[j][ii];
		      a53   = fc(P->GP[0][P->GF[g2][1]],C->ndx[j][ii][k])*dfc(P->GP[0][P->GF[g2][1]],C->NDX[j][ii])   /C->NDX[j][ii];
		      a63   = fc(P->GP[0][P->GF[g2][1]],C->NDX[j][ii])   *dfc(P->GP[0][P->GF[g2][1]],C->ndx[j][ii][k])/C->ndx[j][ii][k];
		      ndxjk = 1.0/C->NDX[j][k];

		      for(n=0,g4=0;n<P->NSYM;n++)
		        if(P->GF[n][0]==4)
		        {
                          m  = g4 + ind4[ds];
			  n4 = W[0].GF[n][4];
			  n6 = W[0].GF[n][6];
			  n7 = W[0].GF[n][7];
			  A2 = a2ijk[n6][n7];
			  A3 = a3ijk[n6][n7];
			  A4 = -2.0*P->GP[3][P->GF[n][4]]*A0;
			  A5 =  A2*a52 + A3*(a53+A4);
			  A6 =           A3*(a63+A4);
			  A7 =  A2*ndxjk;
			  for(q=0;q<3;q++)
			    L->Fn[i][L->DNn[i]][m][q] += n2ijk[n4]*((A5+A6)*dxji[q] + (A7-A6)*dxjk[q]);
			  g4++;
		        }
		    }
	        }
	      L->DNi[i][L->DNn[i]++]=j;
	    }
	  }
        }
        for(n=0;n<P->D;n++)
	  for(q=0;q<6;q++)
	  L->Sn[i][n][q] = 0.0;
        for(j=0;j<C->Nn[i];j++)
        {
          for(q=0;q<3;q++)
            dxij[q] = DX(C,i,j,q);
          for(q=0;q<3;q++)              
            dqij[q] = DX(C,i,j,(q+1)%3);
	  fcij  =  fc(P->GP[0][P->GF[0][1]],C->NDX[i][j]);
	  dfcij = dfc(P->GP[0][P->GF[0][1]],C->NDX[i][j]);
	  s1 = C->ATMN[C->Ni[i][j]];
	  
	  for(n=g2=0;n<P->NSYM;n++)
	    if(P->GF[n][0]==2)
	    {
	      t  = dGNR(P,n,P->GP[0][P->GF[n][1]],C->NDX[i][j],fcij,dfcij)/C->NDX[i][j];
	      ds1 = (s1-s0+P->NSPC)%P->NSPC;
	      m  = g2 + ind2[ds1];

	      if(P->EFS!=2 && L->MRK[i]==1)
	        for(q=0;q<3;q++)
		  L->Fn[i][L->DNn[i]][m][q] -= t*dxij[q];  // note that DX = -( Xi - Xj )
	      // for the pair function in stress calculation only the part on atom i is needed
	      if(P->EFS!=1)
	      {
	        for(q=0;q<3;q++)
		  L->Sn[i][m][q]   += t*dxij[q]*dxij[q];
	        for(q=0;q<3;q++)
		  L->Sn[i][m][q+3] += t*dxij[q]*dqij[q];
	      }
	      g2++;
	    }
          if( g2<P->D && C->NDX[i][j]<P->GP[0][P->GF[g2][1]] )
	  for(k=j+1;k<C->Nn[i];k++)
	    if( C->ndx[i][j][k]<P->GP[0][P->GF[g2][1]] && C->NDX[i][k]<P->GP[0][P->GF[g2][1]] )
	    {
	      for(q=0;q<3;q++)
	        dxik[q] = DX(C,i,k,q);
              for(q=0;q<3;q++)
                dqik[q] = DX(C,i,k,(q+1)%3);

	      s2    = C->ATMN[C->Ni[i][k]];
	      ds1   = (s1-s0+P->NSPC)%P->NSPC;
	      ds2   = (s2-s0+P->NSPC)%P->NSPC;
	      ds    = ds1*ds1 + ds2*ds2;

	      A0   = fc(P->GP[0][P->GF[g2][1]],C->NDX[i][j])*fc(P->GP[0][P->GF[g2][1]],C->NDX[i][k]);
	      T0   = fc(P->GP[0][P->GF[g2][1]],C->ndx[i][j][k]);
	      F0   = dfc(P->GP[0][P->GF[g2][1]],C->NDX[i][j])* fc(P->GP[0][P->GF[g2][1]],C->NDX[i][k])* fc(P->GP[0][P->GF[g2][1]],C->ndx[i][j][k]);
	      F1   =  fc(P->GP[0][P->GF[g2][1]],C->NDX[i][j])*dfc(P->GP[0][P->GF[g2][1]],C->NDX[i][k])* fc(P->GP[0][P->GF[g2][1]],C->ndx[i][j][k]);
	      F2   =  fc(P->GP[0][P->GF[g2][1]],C->NDX[i][j])* fc(P->GP[0][P->GF[g2][1]],C->NDX[i][k])*dfc(P->GP[0][P->GF[g2][1]],C->ndx[i][j][k]);
	      for(n2=0;n2<W[0].GT[3];n2++)
	        n2ijk[n2] = exp(-P->GP[3][n2]*(pow(C->NDX[i][j],2.0)+pow(C->NDX[i][k],2.0)+pow(C->ndx[i][j][k],2.0)));
	      for(z=0;z<W[0].GT[5];z++)
	      for(l=0;l<W[0].GT[6];l++)
	      {
	        a2ijk[z][l] = powz[z]*P->GP[6][l]*P->GP[5][z]*pow(1.0+P->GP[6][l]*C->cos[i][j][k],P->GP[5][z]-1.0)*A0;
	        a3ijk[z][l] = powz[z]*pow(1.0+P->GP[6][l]*C->cos[i][j][k],P->GP[5][z]);
	      }
	      a52 = (1.0/C->NDX[i][k]-C->cos[i][j][k]/C->NDX[i][j])/C->NDX[i][j];
	      a62 = (1.0/C->NDX[i][j]-C->cos[i][j][k]/C->NDX[i][k])/C->NDX[i][k];
	      a53 = fc(P->GP[0][P->GF[g2][1]],C->NDX[i][k])*dfc(P->GP[0][P->GF[g2][1]],C->NDX[i][j])/C->NDX[i][j];
	      a63 = fc(P->GP[0][P->GF[g2][1]],C->NDX[i][j])*dfc(P->GP[0][P->GF[g2][1]],C->NDX[i][k])/C->NDX[i][k];

              for(n=g4=0;n<P->NSYM;n++)
                if(P->GF[n][0]==4)
	        {
                  m  = g4 + ind4[ds];
		  B1 = n2ijk[W[0].GF[n][4]];
		  A1 = B1*T0;
		  n4 = W[0].GF[n][4];
		  n6 = W[0].GF[n][6];
		  n7 = W[0].GF[n][7];
		  A2 = a2ijk[n6][n7];
		  A3 = a3ijk[n6][n7];
		  A4 = -2.0*P->GP[3][P->GF[n][4]]*A0;
		  B2 = A3*A4*T0;
		  B3 = A2*T0;
		  A5 = A1*( A2*a52 + A3*(a53+A4) );
		  A6 = A1*( A2*a62 + A3*(a63+A4) );		
		  if(P->EFS!=2 && L->MRK[i]==1)
		    for(q=0;q<3;q++)
		      L->Fn[i][L->DNn[i]][m][q] -= 2.0*(A5*dxij[q] + A6*dxik[q]);
		  if(P->EFS!=1)
		  {
		    for(q=0;q<3;q++)
		    {
		      L->Sn[i][m][q] -= -2.0*B1*(B2-C->cos[i][j][k]*B3/(C->NDX[i][j]*C->NDX[i][j])+A3*F0/C->NDX[i][j]   )* dxij[q]*dxij[q];
		      L->Sn[i][m][q] -= -2.0*B1*(B2-C->cos[i][j][k]*B3/(C->NDX[i][k]*C->NDX[i][k])+A3*F1/C->NDX[i][k]   )* dxik[q]*dxik[q];
		      L->Sn[i][m][q] -= -2.0*B1*(B2                                               +A3*F2/C->ndx[i][j][k])*(dxij[q]-dxik[q])*(dxij[q]-dxik[q]);
		      L->Sn[i][m][q] -= -2.0*B1*(B3/(C->NDX[i][j]*C->NDX[i][k]))*2.0*(dxij[q]*dxik[q]);
		    }
		    for(q=0;q<3;q++)
		    {
		      L->Sn[i][m][q+3] -= -2.0*B1*(B2-C->cos[i][j][k]*B3/(C->NDX[i][j]*C->NDX[i][j])+A3*F0/C->NDX[i][j]   )* dxij[q]*dqij[q];
		      L->Sn[i][m][q+3] -= -2.0*B1*(B2-C->cos[i][j][k]*B3/(C->NDX[i][k]*C->NDX[i][k])+A3*F1/C->NDX[i][k]   )* dxik[q]*dqik[q];
		      L->Sn[i][m][q+3] -= -2.0*B1*(B2                                               +A3*F2/C->ndx[i][j][k])*(dxij[q]-dxik[q])*(dqij[q]-dqik[q]);
		      L->Sn[i][m][q+3] -= -2.0*B1*(B3/(C->NDX[i][j]*C->NDX[i][k]))*(dxij[q]*dqik[q]+dxik[q]*dqij[q]);
		    }
		  }
		  g4++;
	      }
	    }
        }
      }
      free_i1D(I);
    }
  }
  if(P->IO==1)
    LNK_OUT(L,o,path,P->EFS,P->D);
}
//==================================================================
double pl(PRS *P, int l, double x)
{
  int i;
  double t;
  
  for(i=l%2,t=0.0;i<=l;i+=2)
    t += P->Pl[l][i]*pow(x,(double)i);
  
  return t;
}
//==================================================================                                                                         
double gn(PRS *P, int n, double r)
{
  double t;
  int i;
  if(r > P->RC)
    return 0.0;
  for(i=0,t=0.0;i<P->GN;i++)
    t += pow(P->RC-r,(double)(i+3.0))*P->GW[n][i];
  return t;
}
//==================================================================
//  TO BE CHECKED AND OPTIMIZED
//==================================================================                                                              
void PRS_PS(PRS *P, Cell *C, LNK *L, int o, char *path)
{
  int i,j,k,n,p,l;
  double *gnij;
  
  L->N  = C->N;
  C->Rc = P->RC;
  gnij  = make_d1D(C->NM);
  LIST(C,0);
  for(i=0;i<C->N;i++)
    for(p=0;p<P->GN*P->LN;p++)
      L->Cn[i][p] = 0.0;
  
  for(i=0;i<C->N;i++)    
    for(n=0;n<P->GN;n++)
    {
      for(j=0;j<C->Nn[i];j++)
	gnij[j] = gn(P,n,C->NDX[i][j]);
      for(l=0;l<P->LN;l++)
	for(j=0;j<C->Nn[i];j++)
	  L->Cn[i][l+P->LN*n] += gnij[j]*gnij[j]*(double)(2*l+1)/(4.0*Pi);
      for(l=0;l<P->LN;l++)
	for(j=0;j<C->Nn[i];j++)
	  for(k=j+1;k<C->Nn[i];k++)
	    L->Cn[i][l+P->LN*n] += gnij[j]*gnij[k]*pl(P,l,C->cos[i][j][k]);
    }
  
  free_d1D(gnij);
  if(P->IO==1)
    LNK_OUT(L,o,path,P->EFS,P->D);
}
//==================================================================
// Simply find Rij for optimization of classical potentials
//==================================================================
void PRS_IJ(PRS *P, Cell *C, LNK *L, int o, char *path)
{
  int i,n,p,N,M[2];

  L->N  = C->N;
  C->Rc = P->RC;

  LIST(C,0);  
  N = P->D/C->NSPC;
  for(i=0;i<C->N;i++)
    for(p=0,M[0]=M[1]=0;p<N;p++)
    {
      n = (C->ATMN[i]+C->ATMN[C->Ni[i][p]])%2;  // 0 for A-A & B-B, 1 for A-B & B-A
      L->Cn[i][M[n]+n*N] = C->NDX[i][p];
      M[n]++;
    }
  if(P->IO==1)
    LNK_OUT(L,o,path,P->EFS,P->D);
}
//==================================================================
// Redirect to proper PARSER
//==================================================================
void PARS_STR(PRS *P, PRS *W, Cell *C, LNK *L, int o, char *path)
{
  int i;

  for(i=0;i<C->N;i++) L->ATMN[i]=C->ATMN[i];
  // Power Spectrum   (PS)
  if(P->DSCR==1)
    PRS_PS(P,C,L,o,path);
  // Parinello-Behler (PB)
  if(P->DSCR==2)
    PRS_BP(P,W,C,L,o,path);    
  // Classical potential (IJ)
  if(P->DSCR==3)
    PRS_IJ(P,C,L,o,path);
}
//==================================================================
//  TO BE WRITTEN: ANALYZE WHAT DATA TO INCLUDE
//==================================================================
double CHCK_FRCS(LNK *L, double FMAX)
{
  // returns: (-1) at least one component exceeds the R->FMAX limit; 
  // or abs(max atomic force component) in the structure
  int    i,q;
  double max = -1.0;

  for(i=0;i<L->N;i++)
    for(q=0;q<3;q++)
      if( fabs(L->F[i][q]) > FMAX )
	return -1.0;
      else if(fabs(L->F[i][q]) > max) 
	max = fabs(L->F[i][q]);
  return max;
}
//==================================================================
void SORT_FIT(int *NFIT, double *EFIT, int *RFIT, double *FFIT, double *EDIF, int n, int N, int TAG, double ECUT, double EMAX, double FMAX, char *add)
{
  int i,k,m;
  double Emin,Emax,El,Em,Eh;
  int M;

  M = (int)((double)N*ECUT);

  for(i=n,Emin = 10^6,Emax =-10^6;i<N+n;i++)
  {
    NFIT[i] = TAG; // first: directory's default TAG assigned to all!
    if(EFIT[i]<Emin)
      Emin = EFIT[i];
    if(EFIT[i]>Emax)
      Emax = EFIT[i];
  }
  for(i=n;i<N+n;i++)
    EDIF[i] = EFIT[i] - Emin;

  if (TAG <  0) // tag will be used; but all other conditions are relaxed!
  {
    Em = 1.0 * Emax;
    for(i=n;i<N+n;i++)
    {
      NFIT[i] = abs(TAG);
      FFIT[i] = 1000000.0;
      RFIT[i] = 1;
    }
    printf("%6d %6d %6d % 18.12lf % 18.12lf % 18.12lf   %s\n",N,N,M,Emin,Emax,Em,add);
    return;    
  }

  if (TAG == 4) // if directory is marked for discarding!
  {
    Em = ECUT * Emax;
    for(i=n;i<N+n;i++)
      NFIT[i] = 4;
    printf("%6d %6d %6d % 18.12lf % 18.12lf % 18.12lf   %s\n",N,0,M,Emin,Emax,Em,add);
    return;
  }

  if( fabs(ECUT-1.0)<1e-14 ) // if ECUT = 1.0
  {
    Em = ECUT * Emax;
    for(i=n;i<N+n;i++)
    {
      NFIT[i] = TAG;
      FFIT[i] = FMAX;
      if(EMAX == 0.0)
	RFIT[i] = 1;
      else if((EFIT[i]-Emin-EMAX)<=0.0)
	RFIT[i] = 1;
      else
	RFIT[i] = 0;
    }
    printf("%6d %6d %6d % 18.12lf % 18.12lf % 18.12lf   %s\n",N,N,M,Emin,Emax,Em,add);
    return;
  }

  El = Emin;
  Eh = Emax;
  
  for(m=0,k=0,Em=0.0;m<N;m++)
  {
    Em = 0.5*(El+Eh);
    for(i=n,k=0;i<N+n;i++)
      if(EFIT[i]<Em)
	k++;
    if(k==M)
      break;
    if(k<M)
      El = Em;
    if(k>M)
      Eh = Em;
  }
  printf("%6d %6d %6d % 18.12lf % 18.12lf % 18.12lf   %s\n",N,k,M,Emin,Emax,Em,add);  

  for(i=n;i<N+n;i++)
  {
    FFIT[i] = FMAX;
    if(EMAX == 0.0)
      RFIT[i] = 1;
    else if((EFIT[i]-Emin-EMAX)<=0.0)
      RFIT[i] = 1;
    else
      RFIT[i] = 0;
    if(EFIT[i]>=Em)
      NFIT[i] = 4; // means discarding!
  }
  return;
}
//==================================================================
// Parse all structures in directory specified as DEPO in setup
//==================================================================
void PARS_DAT(ANN *R, PRS *P, Cell *C, LNK *L)
{
  int    i,j,n,m,nn,k,Nmax,x,ND,*RFIT,*NFIT,NRDF,spc1,spc2,TAG,N,totf;
  double *EFIT,***H,EMAX,t,FMAX,ECUT,*FFIT,fmax,*EDIF;
  char    buf[500],buf2[400],s[700],d[600],dn[2000][400],kw[400];
  FILE    *stamp,*in,*dir,*rtable,*nd,*ve;
  PRS     W[9];   // 2*NSPC+1, so for maximum 3 species one needs 9

  int     rindex,poscars;
  int     *tag;
  time_t  rawtime;
  struct  tm * timeinfo;

  totf    = 0;
  poscars = 0;
  int     MAX=100;  //number of columns in progress bar
  double  prog;
  int     y;
  int     NP[]={0,0,0,0,0}; //defaults for TAG: (1) train+test;(2) train;(3) test; (4) discard
  int     *tag_2,*corr;

  printf("|                           Dataset parsing                           |\n");
  printf("=======================================================================\n\n");

  H = make_d3D(C->NRDF,C->MNT,C->MNT);
  NRDF = 0;
  C->NSPC = R->NSPC;
  for(k=0;k<C->NRDF;k++)
    for(spc1=0;spc1<C->NSPC;spc1++)
      for(spc2=0;spc2<C->NSPC;spc2++)
	H[k][spc1][spc2] = 0.0;

  P->EFS = R->EFS;
  Nmax = 0;
  sprintf(s,"mkdir -p %s/",R->data);
  system(s);

  ve = fopen("ve.dat","w");
  //===== TRAINING POSCARS PARSING =====
  // Creates a list of all POSCAR.0 in ANN/data/*
  sprintf(s,"ls %s/*/ -d | sort",R->depo);
  nd = popen(s,"r");
  
  ND = 0;
  n = 0;
  while(fgets(buf,200,nd))
  {
    sscanf(buf,"%s\n",dn[n++]);
    ND++;
  }
  
  pclose(nd);

  struct Node *tmpnd;
  struct Node *tmp[ND];
  struct Node *tindex;
  struct Node *temp;
  
  int p[ND];
  
  tmpnd = (struct Node*) calloc(1,sizeof(struct Node));
  tmpnd->next = NULL;
  tindex = (struct Node*) calloc(1,sizeof(struct Node));
  tindex->next = NULL;
  
  for(n=0;n<ND;n++)
  {
    tmp[n] = (struct Node*) calloc(1,sizeof(struct Node));
    tmp[n]->next = NULL;
  }
  
  N=0;
  for(n=0;n<ND;n++)
  {
    sprintf(s,"find -L %s -name POSCAR.0 | sort",dn[n]);
    in=popen(s,"r");
    p[n]=0;
    while(fgets(buf,200,in))
    {
      tprintf(&(tmp[n]),buf);
      tprintf(&tmpnd,buf);
      p[n]++;
      N++;
    }
    pclose(in);
  }
  
  RFIT = make_i1D(N);
  NFIT = make_i1D(N);
  EFIT = make_d1D(N);
  FFIT = make_d1D(N);
  EDIF = make_d1D(N);

  printf(" dir         poscars               Emin               Emax               Ecut                  path\n");
  for(n=m=0;n<ND;n++)
  {
    k = 0;      
        temp=tmp[n];
    tgets(buf,200,&temp);
    while(tgets(buf,200,&temp))
    {
      sprintf(d,"%s",buf);
      d[strlen(d)-strlen("POSCAR.0")-1] = 0;
      sprintf(s,"%s/POSCAR.0",d);
      READ_CELL(C,s);
      sprintf(s,"%s/dat.dat",d);
      EFIT[m+k] = 10000000.0;
      if( (in=fopen(s,"r")) != 0 )
      {
	fgets(d,200,in);
	sscanf(d,"%lf\n",&EFIT[m+k]);
	EFIT[m+k] /= (double)C->N;
	fclose(in);
	if( EFIT[m+k]>100000.0 )
	{
	  fprintf(stderr,"ERROR reading %s\n",s);
	  exit(1);
	}
	k++;
      }
      else
      {
	fprintf(stderr,"ERROR opening %s\n",s);
	exit(1);
      }
    }
    printf("%3d  ",n);
    sprintf(buf,"%s/tag",dn[n]);
    ECUT = 0.0;
    EMAX = 0.0;
    FMAX = 0.0;
    TAG  = 0;
    if( (in=fopen(buf,"r")) != 0 )
    {
      for(i=0;i<4;i++)
	if(fgets(s,200,in))
	{
	  if(i==0)
	    sscanf(s,"%d\n",&TAG);
	  else if(i==1)
	    sscanf(s,"%lf\n",&ECUT);
	  else if(i==2)
	    sscanf(s,"%lf\n",&EMAX);
	  else if(i==3)
	    sscanf(s,"%lf\n",&FMAX);
	}
      fclose(in);
    }
    if (TAG  == 0  ) TAG  = 1;
    if (ECUT == 0.0) ECUT = R->ECUT;
    if (EMAX == 0.0) EMAX = R->EMAX;
    if (FMAX == 0.0) FMAX = R->FMAX;
    
    SORT_FIT(NFIT,EFIT,RFIT,FFIT,EDIF,m,p[n],TAG,ECUT,EMAX,FMAX,dn[n]);
    m += p[n];      
  }

  printf("Total %6d POSCAR.0 files are found in %s/*\n\n",N,R->depo);
  
  Build_PRS(P,W,1);    
  C->Rmax = P->RC;
  C->Rmin = 0.99*C->Rmax;
  
  tag = make_i1D(N);    
  
  if(R->seed2<0)             //sequential tags
  {
    for(nn=0,n=0;nn<N;nn++) 
      if(NFIT[nn]>0) 
	tag[nn]=n++;
  }
  else                       //fixed or variable sequence of tags?
  {
    if(R->seed2==0) R->seed2=time(NULL);
    PlantSeeds(R->seed2); 
    for(nn=0;nn<N;nn++) 
      NP[NFIT[nn]]++;
    printf("Structures marked for:    TRAIN=% d    TEST=% d    TRAIN+TEST=% d    DISCARD=% d\n",NP[2],NP[3],NP[1],NP[4]); 
    corr=make_i1D(5);
    corr[1]=NP[2];
    corr[2]=0;
    corr[3]=NP[1]+NP[2];
    for(j=1;j<=3;j++)	
      if(NP[j]>0) 
      {
	tag_2=make_i1D(NP[j]);    
	for(nn=0;nn<NP[j];nn++)
	  do
	  {
	    rindex=0;
	    tag_2[nn]=(int)(NP[j]*Random());
	    for(n=0;n<nn;n++) if(tag_2[nn]==tag_2[n]) rindex=1;
	  }while(rindex==1);
	for(nn=0,n=0;nn<N;nn++) if(NFIT[nn]==j) {tag[nn]=corr[j]+tag_2[n++];}
	free_i1D(tag_2);
      }
    free_i1D(corr);
    PlantSeeds(R->seed);
  }
  
  sprintf(buf,"%s/index.dat",R->data);
  dir=fopen(buf,"w");
  fclose(dir);
  
  rindex=0;
  //=====  progress bar  =====
  printf("PARSED %6.2lf%% [", (double)rindex );
  for (x=0; x<MAX; x++) // previously x=i but i was not initialized
    printf(" ");
  printf("]\n\033[F\033[J");
  //===== File processing =====
  temp=tmpnd;
  tgets(L->path,200,&temp);
  for(nn=n=0;nn<N;nn++,n++)
  {
    tgets(L->path,200,&temp);
    L->path[strlen(L->path)-strlen("POSCAR.0")-1] = 0;
    sprintf(buf,"%s/POSCAR.0",L->path);
    sprintf(buf2,"%s/dat.dat",L->path);
    //sprintf(buf3,"%s/frc.dat",L->path);
    if(access(buf,F_OK)!=-1&&access(buf2,F_OK)!=-1)
    {
      READ_CELL(C,buf);
      if(C->N==0) 
      { 
	fprintf(stderr,"ERROR in parsing %s\n",buf); 
	exit(0); 
      };
      
      if(C->N>Nmax)
	Nmax = C->N;
      L->DE = EDIF[n];
      in=fopen(buf2,"r");
      L->p = 0.0;
      fgets(buf,200,in);
      sscanf(buf,"%lf %lf\n",&L->E,&L->p);
      fgets(buf,200,in);
      sscanf(buf,"%s",kw);
      if(strncmp(kw,"in",2)==0)
      {
	sscanf(buf,"%*s %*s %lf %lf %lf %lf %lf %lf",&L->S[0],&L->S[1],&L->S[2],&L->S[3],&L->S[4],&L->S[5]);
	fgets(buf,200,in);
      }
      sscanf(buf,"%s",kw);
      if(strncmp(kw,"PO",2)==0)
      {
	fgets(buf,200,in);
	for(i=0;i<C->N;i++)
	{
	  fgets(buf,200,in);
	  sscanf(buf,"%lf %lf %lf %lf %lf %lf\n",&t,&t,&t,&L->F[i][0],&L->F[i][1],&L->F[i][2]);
	}
	fgets(buf,200,in);
	fgets(buf,200,in);
	L->N = C->N;
      }
      sscanf(buf,"%s",kw);
      if(strncmp(kw,"PS",2)==0)  // PSTRESS in kB
      {
	sscanf(buf,"%s %lf\n",s,&L->p);
	L->p *= 0.1;                         
      }
      fclose(in);

      fmax = CHCK_FRCS(L,FFIT[nn]);
      if( NFIT[nn]<4 && fmax>-1.0 && RFIT[nn]>0 ) // NFIT = 4 means discarding!
      {
	P->IO = 1;
	if(R->EFS==1||R->EFS==3)
	{
	  L->NF = (int)(P->FMRK*(double)C->N);
	  totf += MARK_CL(C,L->NF,R->seed2);
	  for(i=0;i<C->N;i++)
	    if( C->FRC[i]==1 && VectorLen(L->F[i],3)<R->FMIN)
	      {
		L->MRK[i] = 0;
		L->NF--;
	      }
	    else
	      L->MRK[i] = C->FRC[i];
	}
	PARS_STR(P,W,C,L,tag[nn],R->data);
	fprintf(ve,"% 12.6lf   % 12.6lf   % 12.6lf   % 3d   %s\n",CELL_VOL(C)/(double)C->N,L->E/(double)C->N,fmax,C->N,L->path);
	for(i=0;i<C->N;i++)
	  if(NDX(C,i,0)>C->Rmax)
	    break;
	if(i==C->N)
	{
	  NRDF++;
	  RDF(C,1);
	  for(k=0;k<C->NRDF;k++)
	    for(spc1=0;spc1<C->NSPC;spc1++)
	      for(spc2=0;spc2<C->NSPC;spc2++)
		H[k][spc1][spc2] += C->RDF[k][spc1][spc2];
	}
	sprintf(buf,"%s/index.dat",R->data);
	rtable=fopen(buf,"a");       //Updates index.dat file for list of corresponding directories and r**** files
	fprintf(rtable,"e%06d   %s\n",tag[nn],L->path);
	rindex++;
	fclose(rtable);
	
	//===== updating PROGRESS BAR =====
	prog=100.0*(double)rindex/((double)(N));
	printf("PARSED %6.2lf%% [", prog );
	for (x=0; x<(int)prog; x++)
	  printf("=");
	
	for (x=prog; x<MAX; x++)
	  printf(" ");
	
	printf("]\n\033[F\033[J");
	//=================================
      }
    }
    if(nn<N) poscars=rindex;
  }
  
  sprintf(buf,"sort %s/index.dat",R->data);
  in=popen(buf,"r");
  while(fgets(buf,200,in))
    tprintf(&tindex,buf);
  pclose(in);
  
  sprintf(buf2,"%s/index.dat",R->data);
  dir=fopen(buf2,"w");
  temp=tindex;
  tgets(buf,200,&temp);
  while(tgets(buf,200,&temp))
    fprintf(dir,"%s",buf);
  fclose(dir);
  
  if(rindex==0){fprintf(stderr,"\nError: none of %d  structures is parsed! check energy and volume window\n",N);exit(1);} 
  printf("\nSuccessfully parsed %d POSCAR.0 out of total %d structures!\n",poscars,N);
  
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  sprintf(s,"%s/stamp.dat",R->data);
  stamp=fopen(s,"w");
  fprintf(stamp,"VERS   %s\n" ,VERSION);
  fprintf(stamp,"TIME   %s",asctime (timeinfo));
  fprintf(stamp,"DEPO   %s\n",R->depo);
  fprintf(stamp,"STRC   %d\n",N);
  fprintf(stamp,"PRSE   %d\n",rindex);
  fprintf(stamp,"PEFS   %d\n" ,P->EFS);
  fprintf(stamp,"FMRK   %lf\n",P->FMRK);
  fprintf(stamp,"PRSF   %d\n" ,totf);
  fprintf(stamp,"COMP   %d\n",R->D);
  fprintf(stamp,"NSYM   %d\n",R->NSYM);
  fprintf(stamp,"NSPC   %d\n",R->NSPC);
  strcpy(s,"");
  for(y=0;y<R->NSPC;y++)
  {
    sprintf(d,"%d ",R->SPCZ[y]);
    strcat(s,d);
  }
  fprintf(stamp,"TSPC   %s\n" ,s);
  fprintf(stamp,"NMAX   %d\n" ,Nmax);
  fprintf(stamp,"RAND   %ld\n",R->seed2);
  fprintf(stamp,"ECUT   %lf\n",R->ECUT);
  fprintf(stamp,"EMAX   %lf\n",R->EMAX);
  fprintf(stamp,"FMAX   %lf\n",R->FMAX);
  fclose(stamp);
  
  for(k=0;k<C->NRDF;k++)
    for(spc1=0;spc1<C->NSPC;spc1++)
      for(spc2=0;spc2<C->NSPC;spc2++)
	C->RDF[k][spc1][spc2] = H[k][spc1][spc2]/(double)NRDF;
  Print_RDF_FULL(C,"RDFP.dat");

  free_i1D(NFIT);
  free_d1D(EFIT);
  free_i1D(RFIT);
  free_d1D(FFIT);
  free_d1D(EDIF);
  free_i1D(tag);
  
  fclose(ve);  
  sprintf(s,"mv RDFP.dat  %s/ 2>/dev/null",R->data);
  system(s);  
  sprintf(s,"mv ve.dat   %s/ 2>/dev/null",R->data);
  system(s);  
  sprintf(s,"cp basis    %s/ 2>/dev/null",R->data);
  system(s);
}
//==================================================================
