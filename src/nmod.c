#include "nmod.h"

//========================================================
// Neural Network Job
//========================================================
void NNET_MAIN(ANN *R, PRS *P, Cell *C)
{
  LNK L;

  if(R->NSPC==0)
  {
    fprintf(stderr,"Error: species are not defined in setup file!\n");
    exit(1);
  }
  Build_Cell(C,1);

  L.B = 0;
  Build_LNK(&L,C->N,C->NM,P->D,3);

  if ( R->JOBT/10==3 ) { PARS_DAT(R,P,C,&L); }
  if ( R->JOBT/10==4 ) { TRAN_ANN(R,C);      }
  if ( R->JOBT/10==5 ) { EVAL_ANN(R,P,C,&L); }

  exit(0);
}
//================================================================
// Train function (main routine)
//================================================================
void TRAN_ANN(ANN *R, Cell *C)
{
  char   t[500];
  LNK    *L;
  struct timeval t1, t2;
  struct rusage t3;
  double wall_time,user_time,syst_time,cpus_time;

  if( R->JOBT==41 )  
    R->MIX=1; 
  else 
    R->MIX=0;   ///to use pre-trained elements for multi-components

  printf("|                            Model training                           |\n");
  printf("=======================================================================\n\n");

  // initiate the timing
  gettimeofday(&t1, NULL);  

  ANA_STR(R); 

  Build_ANN(R);
  
  L = (LNK *)malloc(R->STR*sizeof(LNK));
  
  LOAD_LNK(R,C,L);

  INIT_MLP(R);
  ADJT_LNK(R,L);
  TRAN_MLP(R,C,L);
  
  CHCK_ERR(R,L);

  // compute wall and cpu time
  getrusage(RUSAGE_SELF, &t3); 
  gettimeofday(&t2, NULL);
  syst_time  = (double) t3.ru_stime.tv_sec + (double) t3.ru_stime.tv_usec / (double) 1e6; // = kernel time
  user_time  = (double) t3.ru_utime.tv_sec + (double) t3.ru_utime.tv_usec / (double) 1e6; // = user time
  wall_time  = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / (double) 1e6;
  cpus_time  = syst_time + user_time;

  OUT_ANN(R,L,wall_time,cpus_time,syst_time,"out","err-ene.dat","err-frc.dat");

  sprintf(t,"paste %s/basis >> %s/model",R->data,R->otpt);
  system(t);

}
//==================================================================
//
//==================================================================
void CHCK_ERR(ANN *R, LNK *L)
{
  int n,ii,i,NF,q,N1,N2,N3,N4;
  double **F;
  //FILE *out;
  //char s[400];
  
  for(n=NF=0;n<R->N+R->TN;n++)
    if(NF<L[n].N)
      NF = L[n].N;
  F = make_d2D(NF,3);

  R->RT=R->ET=R->RE=R->EE=R->RF=R->EF= 0.0;
  N1=N2=N3=N4=0;
  for(n=0;n<R->N+R->TN;n++)
    if( R->WENE>0 || ( R->WENE<0 && (L[n].W<1.0) ) )
    {
      R->TEeval[n] = ENE_ANN(R,&L[n])/(double)L[n].N;
      R->E[n]      = L[n].E          /(double)L[n].N;      
      if(n<R->N)
      {
	R->RE += pow( L[n].E/(double)L[n].N - R->TEeval[n] ,2.0)*L[n].W;
	R->RT += pow( L[n].E - R->TEeval[n]*(double)L[n].N ,2.0)*L[n].W;
	N1 += L[n].N*L[n].N;
      }
      else
      {
	R->EE += pow( L[n].E/(double)L[n].N - R->TEeval[n] ,2.0)*L[n].W;
	R->ET += pow( L[n].E - R->TEeval[n]*(double)L[n].N ,2.0)*L[n].W;
	N2 += L[n].N*L[n].N;
      }
    }
  R->RE  = sqrt(R->RE/(double)R->N);
  if(N2 > 0)
    R->EE  = sqrt(R->EE/(double)R->TN);
  R->RT *= R->WE;
  R->ET *= R->WE;
  
  if(R->EFS==1||R->EFS==3)
  {
    for(n=0;n<R->N+R->TN;n++)
      if( R->WENE>0 || ( R->WENE<0 && (L[n].W<1.0) ) )
      {
	for(ii=0;ii<L[n].NF;ii++)
	  for(q=0;q<3;q++)
	    F[L[n].Fi[ii]][q] = L[n].F[L[n].Fi[ii]][q];
	FRC_ANN(R,&L[n]);
	
	for(ii=0;ii<L[n].NF;ii++)
	{
	  i = L[n].Fi[ii];
	  for(q=0;q<3;q++)
	    if(n<R->N)
	    {
	      R->RF += pow( L[n].f[i][q] - F[i][q], 2.0)*L[n].W;
	      N3++;
	    }
	    else
	    {
	      R->EF += pow( L[n].f[i][q] - F[i][q], 2.0)*L[n].W;
	      N4++;
	    }
	}
      }
    R->RT += R->RF*R->WF*R->WF;
    R->ET += R->EF*R->WF*R->WF;
    R->RF  = sqrt( R->RF/(double)N3 );
    if (N4 > 0)
      R->EF  = sqrt( R->EF/(double)N4 );
  }
  
  R->RT = sqrt( R->RT /(double)(N1+N3)) ;
  if (N2+N4 > 0)
    R->ET = sqrt( R->ET /(double)(N2+N4)) ;
  
  free_d2D(F,NF);

}
//==================================================================
//
//==================================================================
double CPU_TIME(double ti, char buf[200])
{
  double tf = cpu_time();  
  return tf;
}
//==================================================================
//
//==================================================================
void LOAD_LNK(ANN *R, Cell *C, LNK *L)
{
  int  n,l,i,spc;
  FILE *in;
  char s[400];
  
  sprintf(s,"%s/index.dat",R->data);   //Reading e-files index 
  printf("Loading list of parsed data from %s\n\n",s);
  in=fopen(s,"r");
  for(i=0;i<R->N;i++)  {fscanf(in,"e%d %s\n",&R->train[i],s);}
  for(i=0;i<R->TN;i++) {fscanf(in,"e%d %s\n",&R->test[i],s);}
  fclose(in);

  for(n=0;n<R->STR;n++)
    L[n].B = 0;

  for(n=0;n<R->N;n++)
    LNK_IN(&L[n],R->train[n],R->data);
  for(n=0;n<R->TN;n++)
    LNK_IN(&L[n+R->N],R->test[n],R->data);

  if(R->EFS==1||R->EFS==3)
  {
    for(n=R->NF =0;n<R->N ;n++)
      R->NF  += L[n].NF*3;
    for(n=R->TNF=0;n<R->TN;n++)
      R->TNF += L[n+R->N].NF*3;
  }
  for(n=0,R->Eavg=R->Edev=0.0;n<R->STR;n++)
  {
    R->Edev += pow(L[n].E/(double)L[n].N,2.0);
    R->Eavg += L[n].E/(double)L[n].N;
  }

  if(R->MODT==1)
  {
    for(spc=0;spc<R->NSPC;spc++)
      for(l=0;l<R->D;l++)
	R->Rmax[spc][l] = R->Rmin[spc][l] = L[0].Cn[0][l];

    for(n=0;n<R->STR;n++)
      for(i=0;i<L[n].N;i++)
	for(l=0;l<R->D;l++)
	  {
	    if(L[n].Cn[i][l]>R->Rmax[L[n].ATMN[i]][l])
	      R->Rmax[L[n].ATMN[i]][l] = L[n].Cn[i][l];
	    if(L[n].Cn[i][l]<R->Rmin[L[n].ATMN[i]][l])
	      R->Rmin[L[n].ATMN[i]][l] = L[n].Cn[i][l];
	    R->Rmin[L[n].ATMN[i]][l] = 0.0;
	  }  
    for(spc=0;spc<R->NSPC;spc++)
      for(l=0;l<R->D;l++)
	if(R->Rmax[spc][l]-R->Rmin[spc][l] < 1e-14)
	  printf(" The range for the %d %d component is too small\n",spc,l);
	else
	  R->DR[spc][l] = 1.0/(R->Rmax[spc][l]-R->Rmin[spc][l]);
    for(n=R->DNm=0;n<R->N;n++)
      for(i=0;i<L[n].N;i++)
	if(R->DNm<L[n].DNn[i])
	  R->DNm = L[n].DNn[i];
  }

  R->Eavg /= (double)R->STR;
  R->Edev  = sqrt( R->Edev/(double)R->STR - R->Eavg*R->Eavg );

  for(n=0;n<R->N+R->TN;n++)
  {
    if( L[n].DE < fabs(R->WENE) )
      L[n].W = 1.0;
    else
      L[n].W = fabs(R->WENE)/L[n].DE;
  }
}
//=========================================================
// Rescale inputs once to avoid multiplications by R->DR
//=========================================================
void ADJT_LNK(ANN *R, LNK *L)
{
  int  q,n,l,i,j;

  for(n=0;n<R->STR;n++)
    for(i=0;i<L[n].N;i++)
    {
      for(l=0;l<R->D;l++)
        L[n].Cn[i][l] *= R->DR[L[n].ATMN[i]][l];
      if(R->EFS==1)
	for(j=0;j<=L[n].DNn[i];j++)
	  for(l=0;l<R->D;l++)
	    for(q=0;q<3;q++)
	      L[n].Fn[i][j][l][q] *= R->DR[L[n].DNs[i][j]][l];
    }
}
//=========================================================
// Creating output file after training and testing network
//=========================================================
void OUT_ANN(ANN *R, LNK *L,double w_time, double c_time,double s_time,char *s1,char *s2, char *s3)
{
  char   s[500],sout[300];
  int    i,k,q,n,NF;
  FILE*  out;
  double **F,errf;
  
  sprintf(sout,"%s/err-out.dat",R->otpt);
  printf("\n\n job times: real= %1.2lf   user= %1.2lf   sys= %1.2lf   total CPU= %1.2lf (sec)\n",w_time,c_time-s_time,s_time,c_time);
  out=fopen(sout,"a");
  fprintf(out,"\n\n job times: real= %1.2lf   user= %1.2lf   sys= %1.2lf   total CPU= %1.2lf (sec)\n",w_time,c_time-s_time,s_time,c_time);
  fclose(out);
  if(R->EFS==0)
  {
    printf("\nStand. dev. ENE  % lf %6d\n",R->Edev,R->STR);
    printf("Train error ENE  % lf %6d\n",R->RE,R->N);
    printf("Test  error ENE  % lf %6d\n",R->EE,R->TN);
    out=fopen(sout,"a");
    fprintf(out,"\nStand. dev. ENE  % lf %6d\n",R->Edev,R->STR);
    fprintf(out,"Train error ENE  % lf %6d\n",R->RE,R->N);
    fprintf(out,"Test  error ENE  % lf %6d\n",R->EE,R->TN);
    fclose(out);
  }
  if(R->EFS==1||R->EFS==3)
  {
    printf("\nStand. dev. ENE          % lf                     %6d\n",R->Edev,R->STR);
    printf("Train error ENE FRC TOT  % lf % lf % lf %6d %6d\n",R->RE,R->RF,R->RT,R->N,R->NF);
    printf("Test  error ENE FRC TOT  % lf % lf % lf %6d %6d\n",R->EE,R->EF,R->ET,R->TN,R->TNF);
    out=fopen(sout,"a");
    fprintf(out,"\nStand. dev. ENE          % lf                %6d\n",R->Edev,R->STR);
    fprintf(out,"Train error ENE FRC TOT  % lf % lf % lf %6d %6d\n",R->RE,R->RF,R->RT,R->N,R->NF);
    fprintf(out,"Test  error ENE FRC TOT  % lf % lf % lf %6d %6d\n",R->EE,R->EF,R->ET,R->TN,R->TNF);
    fclose(out);
  }
  
  //=====  Writing Test datasets Energy file  =====
  sprintf(s,"%s/%s",R->otpt,s2);
  out = fopen(s,"w");
  fprintf(out,"(eV/atm)   strct.         DFT              NN               NN-DFT    path\n");
  fprintf(out,"-----------------------------------------------------------------------------------------\n");
  for(i=0;i<R->N;i++)
    fprintf(out,"trn%06d e%06d  %14.6lf   %14.6lf   %14.6lf  %s\n",i+1,(int)R->train[i],R->E[i],R->TEeval[i],R->TEeval[i]-R->E[i],L[i].path);
  fprintf(out,"\n");
  for(i=0;i<R->TN;i++)
    fprintf(out,"tst%06d e%06d  %14.6lf   %14.6lf   %14.6lf  %s\n",i+1,(int)R->test[i],R->E[i+R->N],R->TEeval[i+R->N],R->TEeval[i+R->N]-R->E[i+R->N],L[i+R->N].path);
  fclose(out);

    //=====  Writing Test datasets Force file  =====
  if(R->EFS==1||R->EFS==3)
  {
    sprintf(s,"%s/%s",R->otpt,s3);
    out = fopen(s,"w");
    fprintf(out,"(eV/Ang)   strct.        NN-DFT   path\n");
    fprintf(out,"----------------------------------------------------------\n");
    
    for(n=NF=0;n<R->N+R->TN;n++)
      if(NF<L[n].N)
	NF = L[n].N;
    F = make_d2D(NF,3);
    
    for(n=0;n<R->N+R->TN;n++)
    {
      for(k=0;k<L[n].NF;k++)
	for(q=0;q<3;q++)
	  F[L[n].Fi[k]][q] = L[n].F[L[n].Fi[k]][q];
      
      FRC_ANN(R,&L[n]);
      
      for(k=0,errf=0.0;k<L[n].NF;k++)
      {
	i = L[n].Fi[k];
	for(q=0;q<3;q++)
	  errf += pow( L[n].f[i][q] - F[i][q], 2.0);
      }
      if(n == R->N)
	fprintf(out,"\n");
      fprintf(out,"%s%06d e%06d  %12.5lf  %s\n",(n<R->N?"trn":"tst"),(n<R->N?n:n-R->N)+1,(n<R->N ? (int)R->train[n]:(int)R->test[n-R->N]),(L[n].NF > 0 ? sqrt(errf/(double)(3*L[n].NF)) : -1.0),L[n].path);
    }
    fclose(out);
  }

  SAVE_ANN(R,w_time,c_time);
  
  return;
}
//======================================================
//
//======================================================
void PLT_EVAL(char *path,char *s)
{
  FILE *out;

  out=popen("gnuplot","w");
  fprintf(out,"set terminal pngcairo\n");
  fprintf(out,"set output '%s.png'\n",s);
  fprintf(out,"set multiplot\n");
  fprintf(out,"set nokey\n");
  fprintf(out,"set title \"Test\"\n");
  fprintf(out,"set ylabel \"Energy diff.\"\n");
  fprintf(out,"set xlabel \"Test struc.\"\n");
  fprintf(out,"plot \"%s.dat\" using 1:2\n",s);
  pclose(out);
}
//======================================================
double EVAL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  double E;
  int    i;
  char   buf[200];

  if( R->MODT>1 )
    return ENE_POT(C)/(double)C->N;
  P->IO = 1;
  for(i=0;i<C->N;i++)
    L->MRK[i]=1;
  sprintf(L->path,"%s",buf);
  PARS_STR(P,W,C,L,0,".");
  E = ENE_ANN(R,L)/(double)L->N;

  return E;
}
//======================================================
// Evaluating energies for structures in ANN/eval
//======================================================
void EVAL_ANN(ANN *R, PRS *P, Cell *C, LNK *L)
{
  char   fname[700],s1[600],s[600],buf[200],buf2[200],nntp[20][2],str[200],CMP[10][200];
  FILE   *in,*out;
  FILE   *f1,*f2;
  int    q,i,N,rindex,j,Q,AB[100][10],k,I[1000];
  double r,averr,tenr,tenr2,stdev,dft[100],ann[1000],x[1000],Adft,Bdft,Aann,Bann;
  int    n;
  double *ener;
  int    counter;
  PRS    W[9];
  struct Node *temp;
  char   strtype[200];

  printf("|                           Model evaluation                          |\n");
  printf("=======================================================================\n\n");

  if(strlen(R->eval)==1)
  {
    printf("Error in \"setup\" file EVAL directory not specified!\n");
    exit(1);
  }

  sprintf(nntp[ 1],"M");
  sprintf(nntp[ 2],"G");
  sprintf(nntp[ 3],"S");
  sprintf(nntp[ 4],"L");
  
  if(R->NSPC==1) Q = 1;  // elemental systems
  if(R->NSPC==2) Q = 3;  // binary systems
  if(R->NSPC==3) Q = 7;  // ternary systems

  //===== CAUTION does not check setup with 'model' file for NCMP and ... ! =====
  R->STR=1;
  C->ND = 3;
  if(R->MODT==1)
  {
    Build_ANN(R);
    READ_ANN(R);
    Build_PRS(P,W,0);
  }
  else
    READ_POT(C,R->otpt);

  Adft = Bdft = Aann = Bann = 0.0;
  
  // ===== EVALUATE EOS (eval/EOS) =====
  for(q=0;q<Q;q++)
  {
    if(Q<7)
    {
      if(q==0) atom_symb(R->SPCZ[0],CMP[0]);
      if(q==1) atom_symb(R->SPCZ[1],CMP[1]);
      if(q==2) sprintf(CMP[2],"%s%s",CMP[0],CMP[1]);
    }
    else
    {
      if(q==0) atom_symb(R->SPCZ[0],CMP[0]);
      if(q==1) atom_symb(R->SPCZ[1],CMP[1]);
      if(q==2) atom_symb(R->SPCZ[2],CMP[2]);
      if(q==3) sprintf(CMP[3],"%s%s",CMP[0],CMP[1]);
      if(q==4) sprintf(CMP[4],"%s%s",CMP[0],CMP[2]);
      if(q==5) sprintf(CMP[5],"%s%s",CMP[1],CMP[2]);
      if(q==6) sprintf(CMP[6],"%s%s%s",CMP[0],CMP[1],CMP[2]);
    }
    
    struct Node *tmp1;
    tmp1 = (struct Node*) calloc(1,sizeof(struct Node));
    tmp1->next = NULL;

    sprintf(s,"ls -R 2>/dev/null  %s/EOS.%s.*",R->eval,CMP[q]);
    out=popen(s,"r");
    N=0;
    while(fgets(buf,200,out))
    {
      tprintf(&tmp1,buf);
      N++;
    }
    pclose(out);

    sprintf(s1,"mkdir -p %s",R->otpt);
    system(s1);

    if(N>0)
    {     
      printf(" =====  Evaluating EOS   =====\n");
      sprintf(s1,"%s/eos-%s.dat",R->otpt,CMP[q]);
      out=fopen(s1,"w");
      ener=make_d1D(21);
      printf("Total %d EOS files are found in (%s)!\n",N,R->eval); 
      printf(     "   Average DFT energy    Standard deviation     Average error    (eV/atom)\n");
      fprintf(out,"                             Average DFT energy, Standard deviation, Average error (eV/atom)\n");
      temp = tmp1;
      tgets(buf2,200,&temp);
      for(n=0;n<N;n++)
      {
	rindex = 0;
	tenr   = 0.0;
	tenr2  = 0.0;
	stdev  = 0.0;
	averr  = 0.0;
        tgets(buf2,200,&temp);
	      buf2[strlen(buf2)-1]=0;
        sprintf(s1,"%s",buf2);
			in = fopen(s1,"r");
	fgets(str,200,in);
	str[strlen(str)-1]=0;
	if(q<3)
	  sprintf(fname,"%s/%s%s%02d%s.dat",R->otpt,nntp[R->MODT],CMP[q],n,str);      
	if(q==3)
	  sprintf(fname,"%s/%s%s%02d%s.dat",R->otpt,nntp[R->MODT],CMP[q],n+20,str);
	f1 = fopen(fname,"w");
	fgets(buf,200,in);
	fgets(buf,200,in);
	while(fgets(buf,200,in))
	{
	  sscanf(buf,"%lf",&ann[0]);
	  f2 = fopen("tmp","w");
          fgets(buf,200,in);
          fprintf(f2,"%s",buf);
          while(fgets(buf,200,in))
	  {
	    sscanf(buf,"%s",s);
	    fprintf(f2,"%s",buf);
	    if(strncmp(s,"STR",3)==0)
	      break;
	  }
	  fclose(f2);
	  READ_CELL(C,"tmp");
	  ann[1] = EVAL_ENE(R,P,W,C,L);
	  LIST(C,0);
	  for(i=0,r=0.0;i<C->N;i++)
	    r += NDX(C,i,0)/(double)C->N;
	  if(strncmp(str,"nps",3)==0)
	    r = (double)C->N;
	  
	  ann[0]=ann[0]/(double)C->N;//DFT energy
	  ener[rindex]=ann[0];
	  averr+=pow(ann[1] -ann[0],2.0);
	  fprintf(f1,"% lf % lf % lf % lf\n",r,ann[0],ann[1] ,CELL_VOL(C)/(double)C->N);
	  rindex++;      	   	  
	}
	fclose(f1);
	fclose(in);
	for(i=0;i<rindex;i++) 
	  tenr+=ener[i]; 
	tenr=tenr/(double)rindex;
	for(i=0;i<rindex;i++) 
	  tenr2+=pow(ener[i]-tenr,2.0); 
	stdev=sqrt(tenr2/(double)rindex);
	printf(" %24.16lf   %24.16lf   %24.16lf   %s\n",tenr,stdev,sqrt(averr/(double)rindex),fname);
	fprintf(out," %24.16lf   %24.16lf   %24.16lf   %s\n",tenr,stdev,sqrt(averr/(double)rindex),fname);
	sprintf(s,"%s%02d%s.dat",nntp[R->MODT],n,str);
      }
      fclose(out);
      free_d1D(ener);
      system("rm -f e000000 tmp");

      if(q<3)
      {
	PLOT(R,0, 0, 6,21,"eos",CMP[q]);
	if(N>6&&N<9)
	  PLOT(R,0, 6, 7,21,"cls",CMP[q]);
	if(N==9)
	  PLOT(R,0, 6, 9,21,"cls",CMP[q]);
      }
      if(q==3)	
        PLOT(R,0,20,25,21,"eos",CMP[q]);
    }

    //=====  EVALUATE VACANCY FORMATION ENERGIES =====
    struct Node* tmp2;
    tmp2 = (struct Node*) calloc(1,sizeof(struct Node));
    tmp2->next = NULL;

    sprintf(s,"ls -R 2>/dev/null  %s/VAC.%s.*",R->eval,CMP[q]);
    out=popen(s,"r");
    N=0;
    while(fgets(buf,200,out))
    {
      tprintf(&tmp2,buf);
      N++;
    }
    pclose(out);

    if(N>0)
    {
      counter=0;
      printf(" =====  Evaluating VAC   =====\n");
      printf("Total %d VAC files are found in (%s)!\n",N,R->eval);
      temp=tmp2;
      tgets(buf2,200,&temp);
      sprintf(fname,"%s/vac-%s.dat",R->otpt,CMP[q]);
      f1 = fopen(fname,"w");
      
      for(n=0;n<N;n++)
      {
      	j = 0;
        tgets(buf2,200,&temp);
	buf2[strlen(buf2)-1]=0;
	sprintf(s1,"%s",buf2);
	in = fopen(s1,"r");
	fgets(str,200,in);
	sscanf(str,"%s\n",strtype);
	fgets(str,200,in);
	str[strlen(str)-1]=0;
	while(fgets(buf,200,in))
        {
	  sscanf(buf,"%lf",&dft[j]);
	  f2 = fopen("tmp","w");
          fgets(buf,200,in);
          fprintf(f2,"%s",buf);
          while(fgets(buf,200,in))
	  {
	    sscanf(buf,"%s",s);
	    fprintf(f2,"%s",buf);
	    if(strncmp(s,"STR",3)==0)
	      break;
	  }
	  fclose(f2);
	  READ_CELL(C,"tmp");
	  dft [j] /= (double)C->N;
	  ann[j] = EVAL_ENE(R,P,W,C,L);
	  j++;
	}	
	printf(    "% 3d  % lf % lf % lf   % lf % lf % lf   %s %3d\n",n,dft[1],dft[0],(dft[1]-dft[0])*(double)C->N,ann[1],ann[0],(ann[1]-ann[0])*(double)C->N,strtype,C->N);
	fprintf(f1,"% 3d  % lf % lf %s % 8.3lf\n",n,(dft[1]-dft[0])*(double)C->N,(ann[1]-ann[0])*(double)C->N,strtype,(ann[1]-ann[0])*(double)C->N-(dft[1]-dft[0])*(double)C->N);
	counter++;
	fclose(in);
      }
      system("rm -f e000000 tmp");
      fclose(f1);
      PLOT(R,1,0,counter,21,"vac",CMP[q]);
    }

    //===== EVALUATE SURFACE FORMATION ENERGIES =====
    struct Node *tmp3;
    tmp3 = (struct Node*) calloc(1,sizeof(struct Node));
    tmp3->next = NULL;

    sprintf(s,"ls -R 2>/dev/null  %s/SRF.%s*",R->eval,CMP[q]);
    out=popen(s,"r");
    N=0;
    while(fgets(buf,200,out))
    {
      tprintf(&tmp3,buf);
      N++;
    }
    pclose(out);

    if(N>0)
    {
      counter=0;
      printf(" =====  Evaluating SRF   =====\n");
      printf("Total %d SRF files are found in (%s)!\n",N,R->eval);
      temp=tmp3;
      tgets(buf2,200,&temp);
      sprintf(fname,"%s/srf-%s.dat",R->otpt,CMP[q]);
      f1 = fopen(fname,"w");
      for(n=0;n<N;n++)
      {
	j = 0;
        tgets(buf2,200,&temp);
	      buf2[strlen(buf2)-1]=0;
	      sprintf(s1,"%s",buf2);
	      in = fopen(s1,"r");
	      fgets(str,200,in);
	      sscanf(str,"%s\n",strtype);
	      fgets(str,200,in);
	      str[strlen(str)-1]=0;
	      while(fgets(buf,200,in))
	      {	  
	        sscanf(buf,"%lf",&dft[j]);
	        f2 = fopen("tmp","w");
	        fgets(buf,200,in);
	        fprintf(f2,"%s",buf);
	        while(fgets(buf,200,in))
	        {
	          sscanf(buf,"%s",s);
	          fprintf(f2,"%s",buf);
	          if(strncmp(s,"STR",3)==0)
	            break;
	        }
	        fclose(f2);
	        READ_CELL(C,"tmp");
	        x[j] = C->L[0][0]*C->L[1][1]; // assuming that C->L[0] is along x
	        dft [j] /= (double)C->N;
	        ann[j] = EVAL_ENE(R,P,W,C,L);
	        j++;
	      }
	      for(i=1;i<j;i++)
	{
	  printf(    "% 3d  % lf % lf % lf   % lf % lf % lf   %s %3d\n",counter,dft[i],dft[0],(dft[i]-dft[0])*(double)C->N/x[i]*0.5,ann[i],ann[0],(ann[i]-ann[0])*(double)C->N/x[i]*0.5,strtype,C->N);
	  fprintf(f1,"% 3d  % lf % lf %s % 8.3lf\n",counter,(dft[i]-dft[0])*(double)C->N/x[i]*0.5,(ann[i]-ann[0])*(double)C->N/x[i]*0.5,strtype,(ann[i]-ann[0])*(double)C->N/x[i]*0.5-(dft[i]-dft[0])*(double)C->N/x[i]*0.5);
	  counter++;
	}
	fclose(in);
      }
      system("rm -f e000000 tmp");
      fclose(f1);
      PLOT(R,4,0,counter,21,"srf",CMP[q]);
    }

  }
  if(R->NSPC>1)
  for(q=0;q<Q;q++)
  if( (Q==3&&q==2) || (Q==7&&q>2&&q<6) )
  {
    //===== EVALUATE SUBSTITUTIONAL DEFECTS E'S =====                      
    struct Node* tmp4;
    tmp4 = (struct Node*) calloc(1,sizeof(struct Node));
    tmp4->next = NULL;

    sprintf(s,"ls -R 2>/dev/null  %s/SUB.%s*",R->eval,CMP[q]);
    out=popen(s,"r");
    N=0;
    while(fgets(buf,200,out))
    {
      tprintf(&tmp4,buf);
      N++;
    }
    pclose(out);
   
    if(N>0)
    {
      counter=0;
      printf(" =====  Evaluating SUB   =====\n");
      printf("Total %d SUB files are found in (%s)!\n",N,R->eval);

      temp=tmp4;
      tgets(buf2,200,&temp);

      sprintf(fname,"%s/sub-%s.dat",R->otpt,CMP[q]);
      f1 = fopen(fname,"w");    
      for(n=0;n<N;n++)
      {
	j = 0;
        tgets(buf2,200,&temp);
	buf2[strlen(buf2)-1]=0;
	sprintf(s1,"%s",buf2);
	in = fopen(s1,"r");
	fgets(str,200,in);
	sscanf(str,"%s\n",strtype);
	fgets(str,200,in);
	str[strlen(str)-1]=0;
	while(fgets(buf,200,in))
        {
	  sscanf(buf,"%lf",&dft[j]);
	  f2 = fopen("tmp","w");
          fgets(buf,200,in);
          fprintf(f2,"%s",buf);
          while(fgets(buf,200,in))
	  {
	    sscanf(buf,"%s",s);
	    fprintf(f2,"%s",buf);
	    if(strncmp(s,"STR",3)==0)
	      break;
	  }
	  fclose(f2);
	  READ_CELL(C,"tmp");
	  for(i=0;i<C->NSPC;i++)	  
	    AB[j][i] = C->SPCN[i];
	  dft[j] /= (double)C->N;
	  ann[j] = EVAL_ENE(R,P,W,C,L);
	  j++;
	}
	for(k=2;k<4;k++)
        {
	  dft[k] = dft[k]*(double)(AB[k][0]+AB[k][1]) - dft[0]*(double)AB[k][0] - dft[1]*(double)AB[k][1];
	  ann[k] = ann[k]*(double)(AB[k][0]+AB[k][1]) - ann[0]*(double)AB[k][0] - ann[1]*(double)AB[k][1];
	}
	printf(    "% 3d    % lf % lf  % lf % lf %s  % 8.3lf % 8.3lf\n",counter,dft[2],ann[2],dft[3],ann[3],strtype,ann[2]-dft[2],ann[3]-dft[3]);
	fprintf(f1,"% 3d    % lf % lf  % lf % lf %s  % 8.3lf % 8.3lf\n",counter,dft[2],ann[2],dft[3],ann[3],strtype,ann[2]-dft[2],ann[3]-dft[3]);
	counter++;
	fclose(in);
      }
      system("rm -f e000000 tmp");
      fclose(f1);
      PLOT(R,2,0,counter,21,"sub",CMP[q]);
    }


    //===== EVALUATE PHASE DIAGRAMS (HFORM) =====
    struct Node* tmp5;
    tmp5 = (struct Node*) calloc(1,sizeof(struct Node));
    tmp5->next = NULL;

    sprintf(s,"ls -R 2>/dev/null  %s/PHD.%s*",R->eval,CMP[q]);
    out=popen(s,"r");
    N=0;
    while(fgets(buf,200,out))
    {
      tprintf(&tmp5,buf);
      N++;
    }
    pclose(out);

    if(N>0)
    {
      counter=0;
      printf(" =====  Evaluating PHD   =====\n");
      printf("Total %d PHD files are found in (%s)!\n",N,R->eval);
      temp=tmp5;
      tgets(buf2,200,&temp);
      for(n=0;n<N;n++)
      {
	printf("\n");
	sprintf(fname,"%s/phd-%s%02d.dat",R->otpt,CMP[q],n);
	f1 = fopen(fname,"w");
	j = 0;
        tgets(buf2,200,&temp);
	      buf2[strlen(buf2)-1]=0;
        sprintf(s1,"%s",buf2);
	in = fopen(s1,"r");
	fgets(str,200,in);
	fgets(str,200,in);
	str[strlen(str)-1]=0;
	while(fgets(buf,200,in))
	{
	  sscanf(buf,"%lf",&dft[j]);
	  f2 = fopen("tmp","w");	  
	  fgets(buf,200,in);
	  fprintf(f2,"%s",buf);
	  while(fgets(buf,200,in))
	  {
	    sscanf(buf,"%s",s);
	    fprintf(f2,"%s",buf);
	    if(strncmp(s,"STR",3)==0)
	      break;
	  }
	  fclose(f2);
	  READ_CELL(C,"tmp");
	  for(i=0;i<C->NSPC;i++)
	    AB[j][i] = C->SPCZ[i];
	  // to be generalized
	  if(j==1)
	  {
	    AB[j][1] = AB[j][0];
	    AB[j][0] = 0;
	  }
	  x[j] = (double)AB[j][1]/(double)(AB[j][0]+AB[j][1]);
	  dft[j] /= (double)C->N;
	  ann[j] = EVAL_ENE(R,P,W,C,L);
	  printf("%3d %3d %3d % lf % lf\n",j,AB[j][0],AB[j][1],dft[j],ann[j]);
	  Adft = dft[0]; Bdft = dft[1];
	  Aann = ann[0]; Bann = ann[1];
	  j++;
	}
	Sort(x,I,j);
	for(k=0;k<j;k++)
	{
	  dft[I[k]] -= (Adft*(double)AB[I[k]][0] + Bdft*(double)AB[I[k]][1])/(double)(AB[I[k]][0]+AB[I[k]][1]);
	  ann[I[k]] -= (Aann*(double)AB[I[k]][0] + Bann*(double)AB[I[k]][1])/(double)(AB[I[k]][0]+AB[I[k]][1]);
	  printf(    "% lf  % lf  % lf   %3d\n",x[I[k]],dft[I[k]],ann[I[k]],I[k]);
	  fprintf(f1,"% lf  % lf  % lf   %3d\n",x[I[k]],dft[I[k]],ann[I[k]],I[k]);
	  counter++;
	}
	fclose(in);
	fclose(f1);
      }
      PLOT(R,3,0,counter,21,"phd",CMP[q]);
    }
    
  }
  system("rm -f e000000 tmp");

}
//==================================================================
