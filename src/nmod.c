#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "cdef.h"
#include "edef.h"
#include "ndef.h"
#include "cutl.h"
#include "cell.h"
#include "cpot.h"
#include "cmod.h"
#include "nprs.h"
#include "nutl.h"
#include "util.h"
#include "nmlp.h"
#include "plot.h"
#include "nmod.h"
#include "nmin.h"

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
  double time1,time2;
  char   t[200];
  LNK    *L;

  if( R->JOBT==41 )  
    R->MIX=1; 
  else 
    R->MIX=0;   ///to use pre-trained elements for multi-components
  
  printf(" =====   Training ANN ...    =====\n");

  ANA_STR(R); 
  
  time1=cpu_time();
  Build_ANN(R);
  READ_ATM(R,C);
  
  L = (LNK *)malloc(R->STR*sizeof(LNK));
  
  LOAD_LNK(R,C,L);

  INIT_MLP(R);
  TRAN_MLP(R,C,L);
  
  CHCK_ERR(R,L);
  time2=cpu_time();

  OUT_ANN(R,time2-time1,"out","energy_errors.dat","error_results.png");

  sprintf(t,"paste %s/basis >> %s/model",R->data,R->otpt);
  system(t);

}
//==================================================================
//
//==================================================================
void CHCK_ERR(ANN *R, LNK *L)
{
  int n,ii,i,j,NF,q,N1,N2,N3,N4;
  double **F;
  FILE *out;
  char s[200];
  
  for(n=NF=0;n<R->N+R->TN;n++)
    if(NF<L[n].N)
      NF = L[n].N;
  F = make_d2D(NF,3);
  
  R->RT=R->ET=R->RE=R->EE=R->RF=R->EF= 0.0;
  N1=N2=N3=N4=0;
  for(n=0;n<R->N+R->TN;n++)
    {
      R->TEeval[n] = ENE_ANN(R,&L[n])/(double)L[n].N;
      R->E[n]      = L[n].E          /(double)L[n].N;      
      if(n<R->N)
	{
	  R->RE += pow( L[n].E/(double)L[n].N - R->TEeval[n] ,2.0);
	  R->RT += pow( L[n].E - R->TEeval[n]*(double)L[n].N ,2.0);
	  N1 += L[n].N*L[n].N;
	}
      else
	{
	  R->EE += pow( L[n].E/(double)L[n].N - R->TEeval[n] ,2.0);
	  R->ET += pow( L[n].E - R->TEeval[n]*(double)L[n].N ,2.0);
	  N2 += L[n].N*L[n].N;
	}
    }
  R->RE  = sqrt(R->RE/(double)R->N);
  R->EE  = sqrt(R->EE/(double)R->TN);
  R->RT *= R->WE;
  R->ET *= R->WE;
  
  if(R->EFS==1||R->EFS==3)
  {
    for(n=0;n<R->N+R->TN;n++)
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
	    R->RF += pow( L[n].f[i][q] - F[i][q], 2.0);
	    N3++;
	  }
	  else
	  {
	    R->EF += pow( L[n].f[i][q] - F[i][q], 2.0);
	    N4++;
	  }
      }
    }
    R->RT += R->RF*R->WF*R->WF;
    R->ET += R->EF*R->WF*R->WF;
    R->RF  = sqrt( R->RF/(double)N3 );
    R->EF  = sqrt( R->EF/(double)N4 );
  }
  
  R->RT = sqrt( R->RT /(double)(N1+N3)) ;
  R->ET = sqrt( R->ET /(double)(N2+N4)) ;
  
  free_d2D(F,NF);
  sprintf(s,"%s/error.dat",R->otpt);
  out = fopen(s,"w");
  if(R->EFS==1)
    fprintf(out,"% 5d % 5d   % lf % lf   % lf % lf  % lf % lf\n",R->M,R->N+N1,R->RE,R->EE,R->RF,R->EF,R->RT,R->ET);
  else
    fprintf(out,"% 5d % 5d   % lf % lf   % lf\n",R->M,R->N,R->RE,R->EE,R->time);
  fclose(out);
  
  if(R->MODT==1)
    for(i=0;i<R->N+R->TN;i++)    
      for(j=0;j<L[i].N;j++)
      {
	R->E[i]       += (R->E0[R->SPCZ[L[i].ATMN[j]]]/(double)L[i].N);
	R->TEeval[i]  += (R->E0[R->SPCZ[L[i].ATMN[j]]]/(double)L[i].N);
      }
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
  int n,l,i,spc;
  FILE *in;
  char s[200];
  
  sprintf(s,"%s/index.dat",R->data);   //Reading e-files index 
  printf("%s\n",s);
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
    for(i=0;i<L[n].N;i++)
      L[n].E -= R->E0[R->SPCZ[L[n].ATMN[i]]];
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
}
//=========================================================
// Creating output file after training and testing network
//=========================================================
void OUT_ANN(ANN *R, double time, char *s1,char *s2, char *s3)
{
  char   s[200];
  int    i;
  FILE*  out;
  
  if(R->EFS==0)
  {
    printf("Stand. dev. ENE  % lf %6d\n",R->Edev,R->STR);
    printf("Train error ENE  % lf %6d\n",R->RE,R->N);
    printf("Test  error ENE  % lf %6d\n",R->EE,R->TN);
  }
  if(R->EFS==1||R->EFS==3)
  {
    printf("Train error ENE FRC TOT  % lf % lf % lf %6d %6d\n",R->RE,R->RF,R->RT,R->N,R->NF);
    printf("Test  error ENE FRC TOT  % lf % lf % lf %6d %6d\n",R->EE,R->EF,R->ET,R->TN,R->TNF);
  }
  
  //=====  Writing Test datasets Energy file  =====
  sprintf(s,"%s/%s",R->otpt,s2);
  out = fopen(s,"w");
  fprintf(out,"Energy Bias=%24.16lf\n------------------------------------\n",R->Ebias);
  fprintf(out,"     Strct.     Target energy              Evaluated energy            Eval-Target energy\n");
  fprintf(out,"-----------------------------------------------------------------------------------------\n");
  for(i=0;i<R->N;i++)
    fprintf(out,"%06d e%06d  %24.16lf   %24.16lf   %24.16lf\n",i+1,(int)R->train[i],R->E[i]-R->Ebias,R->TEeval[i]-R->Ebias,R->TEeval[i]-R->E[i]);
  fprintf(out,"\n");
  for(i=0;i<R->TN;i++)
    fprintf(out,"%06d e%06d  %24.16lf   %24.16lf   %24.16lf\n",i+1,(int)R->test[i],R->E[i+R->N]-R->Ebias,R->TEeval[i+R->N]-R->Ebias,R->TEeval[i+R->N]-R->E[i+R->N]);
  fclose(out);
  
  SAVE_ANN(R,time);

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
  char   fname[200],s1[200],s[200],buf[200],buf2[200],nntp[20][2],str[200],CMP[10][200];
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
  if( R->MODT>1 )
    READ_POT(C,R->otpt);
  else
  {
    Build_ANN(R); 
    READ_ANN(R);
    Build_PRS(P,W,0);
  }
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
      sprintf(s1,"%s/eos_%s.dat",R->otpt,CMP[q]);
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
	  for(i=0,ann[1]=0.0;i<C->N;i++)
	    ann[1] += (R->E0[C->ATMZ[i]]/C->N);
	  ann[1] += EVAL_ENE(R,P,W,C,L);
	  LIST(C);
	  for(i=0,r=0.0;i<C->N;i++)
	    r += NDX(C,i,0)/(double)C->N;
	  if(strncmp(str,"nps",3)==0)
	    r = (double)C->N;
	  
	  ann[0]=ann[0]/(double)C->N;//DFT energy
	  ener[rindex]=ann[0];
	  averr+=pow(ann[1] -ann[0],2.0);
	  fprintf(f1,"% lf % lf % lf % lf\n",r,ann[0],ann[1] ,Cell_VOLUME(C)/(double)C->N);
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
      sprintf(fname,"%s/vac_%s.dat",R->otpt,CMP[q]);
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
	  for(i=0,ann[j]=0.0;i<C->N;i++)
	    ann[j] += (R->E0[C->ATMZ[i]]/C->N);
	  dft [j] /= (double)C->N;
	  ann[j] += EVAL_ENE(R,P,W,C,L);
	  j++;
	}	
	printf(    "% 3d  % lf % lf % lf   % lf % lf % lf   %s %3d\n",n,dft[1],dft[0],(dft[1]-dft[0])*(double)C->N,ann[1],ann[0],(ann[1]-ann[0])*(double)C->N,strtype,C->N);
	fprintf(f1,"% 3d  % lf % lf %s % 8.3lf\n",n,(dft[1]-dft[0])*(double)C->N,(ann[1]-ann[0])*(double)C->N,strtype,(ann[1]-ann[0])*(double)C->N-(dft[1]-dft[0])*(double)C->N);
	counter++;
	fclose(in);
      }
      system("rm -f e000000");
      system("rm tmp");
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
      sprintf(fname,"%s/srf_%s.dat",R->otpt,CMP[q]);
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
	        for(i=0,ann[j]=0.0;i<C->N;i++)
	          ann[j] += (R->E0[C->ATMZ[i]]/C->N);
	        dft [j] /= (double)C->N;
	        ann[j] += EVAL_ENE(R,P,W,C,L);
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
      system("rm -f e000000");
      system("rm tmp");
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

      sprintf(fname,"%s/sub_%s.dat",R->otpt,CMP[q]);
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
	  for(i=0,ann[j]=0.0;i<C->N;i++)
	    ann[j] += (R->E0[C->ATMZ[i]]/C->N);
	  dft[j] /= (double)C->N;
	  ann[j] += EVAL_ENE(R,P,W,C,L);
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
      system("rm -f e000000");
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
	sprintf(fname,"%s/phd_%s%02d.dat",R->otpt,CMP[q],n);
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
	  for(i=0,ann[j]=0.0;i<C->N;i++)
	    ann[j] += (R->E0[C->ATMZ[i]]/C->N);
	  dft[j] /= (double)C->N;
	  ann[j] += EVAL_ENE(R,P,W,C,L);
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
  system("rm -f e000000");

}
//==================================================================
