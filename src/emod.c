#include "emod.h"

extern const double Pi;

ANN *RR;
PRS *PP;
LNK *LL;
PRS  WW[9];

//==================================================================
//     Initialize and generate a new population
//==================================================================
void INIT_TR(Tribe *T)
{
  int i,j,k,p,N,SES[4],FES[4],n;
  double dista,scale;
  char buf[200];

  system("mkdir -p EVOS\n");

  for(n=0;n<4;n++)
  {
    SES[n] = T->SES[n];
    FES[n] = T->FES[n];
  }
  if(FES[3]==T->N)    
    T->FES[0] = T->SES[1] = T->FES[1] = T->SES[2] = T->FES[2] = T->SES[3] = T->FES[3] = 2*T->N; // generate random NPs with NANO_TETR by default
  else
  {
    T->FES[0] = T->SES[1] = T->N      + (int)((double)(FES[0]-SES[0])/(double)(FES[3]-SES[0])*(double)T->N);
    T->FES[1] = T->SES[2] = T->SES[1] + (int)((double)(FES[1]-SES[1])/(double)(FES[3]-SES[0])*(double)T->N);
    T->FES[2] = T->SES[3] = T->SES[2] + (int)((double)(FES[2]-SES[2])/(double)(FES[3]-SES[0])*(double)T->N);
    T->FES[3] =             T->SES[3] + (int)((double)(FES[3]-SES[3])/(double)(FES[3]-SES[0])*(double)T->N);
    if( T->FES[3] != 2*T->N )
    {
      printf("Please make sure that the sum of TETR, PLNT, PACK, and BLOB operations for nanoparticles is a proper fraction of NPOP\n");
      fprintf(stderr,"Please make sure that the sum of TETR, PLNT, PACK, and BLOB operations for nanoparticles is a proper fraction of NPOP\n");
      exit(1);
    }
  }

  for(n=0;n<4;n++)
    if( (T->FES[n]-T->SES[n]) > 0)
      printf("INIT %3d %s %3d %3d %3d\n",n,T->NES[n],T->FES[n]-T->SES[n],T->SES[n],T->FES[n]);

  N = 2*T->N;

  for(k=0;k<T->NSPC;k++)               //scan through array of types of elements
  {
    i = T->SPCZ[k];
    if (i > 95)
    {
      sprintf(buf,"Element %d is not in the library\n",T->SPCZ[k]);
      fprintf(stderr,"Element %d is not in the library\n",T->SPCZ[k]);
      Print_LOG(buf);
      exit(0);
    }
    if( T->p >= T->minpres[0] && T->p <= T->maxpres[0] )
    {
      atom_symb(i,T->SS[k]);
      scale = (T->p-T->minpres[0])/(T->maxpres[0]-T->minpres[0]);
      dista = (T->minpres[i]-T->maxpres[i]);
      T->Rm[k] = (T->minpres[i]-dista*scale)*T->Rhc;
      if(T->Rm[k]<0.00001)
      {
	sprintf(buf,"Equilibrium distance for element %d is not correct\n",T->SPCZ[k]);
	fprintf(stderr,"Equilibrium distance for element %d is not correct\n",T->SPCZ[k]);
	Print_LOG(buf);
	exit(0);
      }
    }
    else
    {
      sprintf(buf,"Pressure %lf for element %d is not in the defined range (% lf-% lf)\n",T->p,T->SPCZ[k],T->minpres[0],T->maxpres[0]);
      fprintf(stderr,"Pressure %lf for element %d is not in the defined range (% lf-% lf)\n",T->p,T->SPCZ[k],T->minpres[0],T->maxpres[0]);
      Print_LOG(buf);
      exit(0);
    }
  }

  for(k=0;k<T->NSPC;k++)
  {
    sprintf(buf,"%3d %3s %3d  % lf % lf\n",T->SPCZ[k],T->SS[k],T->SPCN[k],T->Rm[k]/T->Rhc,T->Rm[k]);
    Print_LOG(buf);
  }

  for(p=0;p<2*T->N+2;p++)
  {
    T->C[p].R0 = 1.0;                  //set equilibrium bond length of C
    T->C[p].NSPC=T->NSPC;
    for(k=0;k<T->NSPC;k++)               //seting the number of each type for C[p]
      T->C[p].SPCN[k] = T->SPCN[k];
    for(k=0;k<T->NSPC;k++)
      T->C[p].SPCZ[k] = T->SPCZ[k];
    for(k=i=0;k<T->NSPC;k++)             //set each atom in C[p] to the appropriate atom type
      for(j=0;j<T->SPCN[k];j++,i++)
	T->C[p].ATMN[i] = k;
    T->C[p].N = i;                     //defining number of atoms for C[p]
    //===== This seems redundant and does not work for CODE=0 =====
    if(0)
    for(k=0;k<T->C[p].N;k++)
      T->C[p].ATMZ[k] = T->C[p].ATMZ[i];
    T->C[p].MNT = T->C[p].NSPC = T->NSPC;  //Setting max number of atom types
  }

  if(T->ND==3)                         //for crystals: define hard sphere volume for all the atoms
    for(i=0,T->VOL=0.0;i<T->NSPC;i++)
      T->VOL += (4.0/3.0*Pi)*pow(T->Rm[i]/0.70,3.0)*(double)T->SPCN[i];
  if(T->ND==2)                         //for thin films: define hard sphere volume for all the atoms
    for(i=0,T->VOL=0.0;i<T->NSPC;i++)
      T->VOL += (1.0/1.0*Pi)*pow(T->Rm[i]/0.70,2.0)*(double)T->SPCN[i];
  
  sprintf(buf,"%lf\n",T->VOL);
  Print_LOG(buf);

  //If this is the first step of a EVOS run then generate structures
  if(T->n<=0)                   // if T->n < 0 the particles will be re-generated for REFL and INVS
    for(p=T->N;p<2*T->N;p++)    //fill the second half of the array of cells (first half is for energy ordering)
    {
      if(T->ND==3)              // generate crystals
      {
	if(T->JS==0)            // start from random structures
        {
	  RAND_CL(T,&T->C[p],&T->C[N],0);
	  SHRT_LV(&T->C[p]);
	}
	else if(T->JS==1)            // start from specified structures
        {
	  sprintf(buf,"INI/POSCAR%03d",0);
	  READ_CELL(&T->C[N],buf);
	  TEMP_CL(T,&T->C[N],p);
	}
	else if(T->JS==2 )            // start from random structures but keep LV constant - fixed atoms too
        {
	  sprintf(buf,"INI/POSCAR%03d",0);
	  READ_CELL(&T->C[p],buf);
	  RAND_CL(T,&T->C[p],&T->C[N],1);
	}
      }
      if(T->ND==2)              // generate films (beta version)
      {
	if(T->JS==0)            // start from random structures
	{
	  RAND_CL(T,&T->C[p],&T->C[N],0);
	  SHRT_LV(&T->C[p]);
	}
        else if(T->JS==1)            // start from specified structures
	{
	  sprintf(buf,"INI/POSCAR%03d",0);
	  READ_CELL(&T->C[N],buf);
	  TEMP_CL(T,&T->C[N],p);
	}
      }
      if(T->ND==0)              // generate particles, cubic lattice cell is always fixed
      {
	if( p>=T->SES[0] && p< T->FES[0] )
	  NANO_TETR(T,0,&T->C[N],p);
	if( p>=T->SES[1] && p< T->FES[1] )
	  NANO_PLNT(T,1,p);
	if( p>=T->SES[2] && p< T->FES[2] )
	  NANO_PACK(T,2,&T->C[N],p);
        if( p>=T->SES[3] && p< T->FES[3] )
          NANO_BLOB(T,3,&T->C[N],p);
      }
      //      LIST(&T->C[p]);
      T->C[p].P=0.0;
      T->C[p].dE = -0.02;
    }

  if(T->CODE!=0)
  {
    sprintf(buf,"cp INI/sub0 INI/sub");
    system(buf);
    sprintf(buf,"sed -i 's/XX:XX:XX/%02d:%02d:%02d/' INI/sub",T->time/3600,(T->time/60)%60,T->time%60);
    system(buf);
  }
  if(T->CODE==1) //===== for VASP =====
  {
    sprintf(buf,"cp INI/INCAR0 INI/INCAR");
    system(buf);
    sprintf(buf,"sed -i 's/P.PPP/%lf/' INI/INCAR",T->p*10.0);
    system(buf);
  }

  if(T->n>0)  //If this is not the first step in the EVOS then read in from last step
    for(p=0;p<2*T->N;p++)
    {
      sprintf(buf,"EVOS/G%03d/M%03d/CONTCAR.0",T->n-1,p);
      if(READ_CELL(&T->C[p],buf)==0)
      {
	printf("ERROR: specified generation to restart from does not exist\n");
	exit(0);
      }
      sprintf(buf,"EVOS/G%03d/M%03d/OSZICAR.0",T->n-1,p);
      T->C[p].P = Read_OSZI(buf);
      if( fabs(T->C[p].P-1.0)<1e-15 )
      {
	printf("ERROR: specified generation to restart from does not exist\n");
	exit(0);
      }
      SHRT_LV(&T->C[p]);
      LIST(&T->C[p],0);      
      T->P1[p] = T->P2[p] = T->C[p].P;
      T->C[p].dE = -0.02;
    }
  for(n=0;n<3;n++)
  {
    T->SES[n] = SES[n];
    T->FES[n] = FES[n];
  }
}
//==================================================================
//     Read setup and exit if needed
//==================================================================
void EXIT_TR(Tribe *T)
{
  int i,p;
  char buf[400],s[200];
  FILE *in;

  if( !(  in = fopen("setup","r")) )
  {
    sprintf(buf,"Failed to open setup\n");
    fprintf(stderr,"Failed to open setup\n");
    Print_LOG(buf);
    exit(1);
  }
  while(fgets(buf,200,in))
    if(strncmp(buf,"JOBT",4)==0) 
      sscanf(buf+4,"%d" ,&T->JOBT);
  fclose(in);

  if(T->JOBT!=12)
    return;

  for(p=T->N;p<T->N*2;p++)
  {
    sprintf(buf,"EVOS/G%03d/M%03d/jobid",T->n,p);
    if( (in = fopen(buf,"r"))!=0 )
    {
      fgets(buf,200,in);
      for(i=0;i<strlen(buf)&&(buf[i]<48||buf[i]>57);i++);
      if(i<strlen(buf))
      {
	sscanf(buf+i,"%s",s);
	if( T->QT==0 )
	  sprintf(buf,"qdel %s",s);
	if( T->QT==1 )
	  sprintf(buf,"scancel %s",s);
	system(buf);
	sprintf(buf,"deleting job %s\n",s);
	Print_LOG(buf);
      }
      fclose(in);
    }
  }
  sprintf(buf,"rm -r EVOS/G%03d/",T->n);
  system(buf);
  sprintf(buf,"Hard exit at iteration %3d \n",T->n);
  fprintf(stderr,"Hard exit at iteration %3d \n",T->n);
  Print_LOG(buf);
  exit(1);
}
//==================================================================
//     Rank Tribe: Sorting is N^2
//==================================================================
void RANK_TR(Tribe *T)
{
  int i,k,p,*I,M;
  double Pmin,Pmax,t;
  Pmin =  10000.0;
  Pmax = -10000.0;  
  I = make_i1D(2*T->N);
  char buf[200];

  for(p=0;p<2*T->N;p++)
  {
    T->E[p] = T->C[p].P;
    if(T->E[p]<Pmin)
      Pmin = T->E[p];
    if(T->E[p]>Pmax)
      Pmax = T->E[p];
  }

  if( Pmax-Pmin < 1e-12 )
  {
    sprintf(buf,"The energy diversity is too small % lf \n",Pmax-Pmin);
    fprintf(stderr,"The energy diversity is too small % lf \n",Pmax-Pmin);
    Print_LOG(buf);
    exit(1);
  }

  for(p=0;p<2*T->N;p++)
    T->f[p] = 0.5*( 1.0 - tanh(2.0*(T->E[p]-Pmin)/(Pmax-Pmin)-1.0) );
  Sort(T->E,I,2*T->N);
  for(p=0;p<2*T->N;p++)
    T->S[p] = I[p];
  for(p=0;p<2*T->N;p++)
  {
    Copy_C(&T->C[    p ],&T->C[2*T->N]);
    Copy_C(&T->C[  I[p]],&T->C[    p ]);
    Copy_C(&T->C[2*T->N],&T->C[  I[p]]);

    dSwap(&T->f[p], &T->f[I[p]]);
    dSwap(&T->E[p], &T->E[I[p]]);
    dSwap(&T->P1[p],&T->P1[I[p]]);
    dSwap(&T->P2[p],&T->P2[I[p]]);

    for(k=0;k<2*T->N;k++)
      if(I[k]==p)
	I[k] = I[p];
  }  

  for(p=0;p<2*T->N;p++)
  {
    LIST(&T->C[p],0);
    RDF(&T->C[p],1);
  }

  for(p=2*T->N-1;p>=1;p--)
    for(k=p-1;k>=0;k--)
      if( (t=COMP_CL(&T->C[p],&T->C[k])) > T->CUT)
      {
        sprintf(buf,"found a duplicate %3d %3d % 24.16lf\n",k,p,t);
	Print_LOG(buf);
	T->f[p] *= 0.001 + (1.0-t)*(1.0-t);
	break;
      }  
  for(p=M=0;p<2*T->N;p++)
    if( T->f[p]>1e-10 )
    {
      M++;
      if( M==2*T->N-(int)(2.0*T->HE*(double)T->N) )
	break;
    }
  sprintf(buf,"discarded all above %3d % 15.6lf\n",p,T->E[p]);
  Print_LOG(buf);
  for(;p<2*T->N;p++)
    T->f[p] = 1e-15;

  if(1)
    for(p=0;p<2*T->N;p++)
    {
      sprintf(buf,"%3d % 24.16lf % 24.16lf\n",p,T->E[p],T->f[p]);
      Print_LOG(buf);
    }
  for(i=0;i<T->C[0].N;i++)
    T->C[0].ATMN[i] = T->C[0].ATMN[i];
  SAVE_CELL(&T->C[0],"poscar",0);
  if(access("INI/",F_OK)!=-1)
    SAVE_CELL(&T->C[0],"INI/POSCAR000",0);
  free_i1D(I);
}
//==================================================================
//     Select Tribe
//==================================================================
void SLCT_TR(Tribe *T)
{
  int k,p;
  double *f,r;
  char buf[200];
  
  f = make_d1D(2*T->N+1);

  for(p=1,f[0]=0.0;p<2*T->N+1;p++)
    f[p] = f[p-1] + T->f[p-1];

  for(p=0;p<2*T->N+1;p++)
    f[p] /= f[2*T->N];
  for(p=0;p<2*T->N;p++)
    T->I[T->n][ p] = 0;
  for(p=0;p<T->NB;p++)
    T->I[T->n][ p] = 1;

  for(k=T->NB;k<T->N;k++)
  {
    r = Random();
    for(p=1;p<2*T->N&&f[p]<r;p++);
    if(T->I[T->n][ p-1]==1)
      k--;
    else
      T->I[T->n][ p-1] = 1;
  }    
  for(p=0;p<2*T->N;p++)
    printf("% 3d %3d % lf %3d % lf\n",p,T->I[T->n][p],f[p],T->S[p],T->E[p]);

  for(k=p=0;p<2*T->N;p++)
    if(T->I[T->n][ p]==1)
    {
      Copy_C(&T->C[p],&T->C[k]);
      T->S[ k] = T->S[ p];
      T->P1[k] = T->P1[p];
      T->P2[k] = T->P2[p];
      T->E[ k] = T->E[ p];
      T->f[ k] = T->f[ p];
      sprintf(buf,"%3d %3d %3d % 15.6lf\n",k,p,T->S[k],T->E[k]);
      Print_LOG(buf);
      k++;
    }
  free_d1D(f);
}
//==================================================================
//     Plot Tribe
//==================================================================
void PLOT_TR(Tribe *T)
{
  int p;

  FILE *out1, *out2, *out3, *out4;

  if(T->n==0)
    system("rm -f erank.dat elink.dat ebest.dat eplot.dat");
  out1 = fopen("erank.dat","a");
  out2 = fopen("elink.dat","a");
  out3 = fopen("ebest.dat","a");
  out4 = fopen("eplot.dat","a");

  for(p=0;p<T->N;p++)
    fprintf(out1,"%4d % 15.9lf\n",T->n,T->E[p]/(double)T->C[p].N);
  if(T->n>0)
    for(p=0;p<T->N;p++)
      fprintf(out2,"%4d % 15.9lf\n%4d % 15.9lf\n\n%4d % 15.9lf\n%4d % 15.9lf\n\n",T->n-1,T->P1[p]/(double)T->C[p].N,T->n,T->E[p]/(double)T->C[p].N,T->n-1,T->P2[p]/(double)T->C[p].N,T->n,T->E[p]/(double)T->C[p].N);
  fprintf(out3,"%4d % 15.9lf\n",T->n,T->E[0]/(double)T->C[0].N);
  for(p=0;p<T->N;p++)
    fprintf(out4,"% 15.9lf % 15.9lf %4d %4d\n",CELL_VOL(&T->C[p])/(double)T->C[p].N,T->E[p]/(double)T->C[p].N,T->n,p);

  for(p=0;p<2*T->N;p++)
    T->P1[p] = T->P2[p] = T->E[p];

  if(T->JOBT==14)
    printf("Analyzed  %3d structures in EVOS/G%03d\n",T->N,T->n);

  fclose(out1);
  fclose(out2);
  fclose(out3);
  fclose(out4);
}
//==================================================================
//     Relax Tribe
//==================================================================
void QSUB_TR(Tribe *T, int p)
{
  char buf[200];

  sprintf(buf,"mkdir -p EVOS/G%03d/M%03d",T->n,p);
  system(buf);
  if(T->CODE==1) //===== for VASP =====
  {
    sprintf(buf,"cp INI/POTCAR EVOS/G%03d/M%03d",T->n,p);
    system(buf);
    sprintf(buf,"cp INI/INCAR EVOS/G%03d/M%03d",T->n,p);
    system(buf);
    sprintf(buf,"EVOS/G%03d/M%03d/KPOINTS",T->n,p);
    KMESH(&T->C[p],T->KM,buf,T->ND);
  }
  sprintf(buf,"EVOS/G%03d/M%03d/POSCAR",T->n,p);
  SAVE_CELL(&T->C[p],buf,0);
  sprintf(buf,"cp INI/sub g%03d",p);
  system(buf);
  sprintf(buf,"sed -i 's/GGGG/G%03d/' g%03d",T->n,p);
  system(buf);
  sprintf(buf,"sed -i 's/MMMM/M%03d/' g%03d",p,p);
  system(buf);
  if( T->JOBT==14 )
    return;
  if( T->QT==0 )
    sprintf(buf,"qsub g%03d >> EVOS/G%03d/M%03d/jobid",p,T->n,p);
  if( T->QT==1 )
    sprintf(buf,"sbatch g%03d >> EVOS/G%03d/M%03d/jobid",p,T->n,p);
  if( T->QT==2 )
    sprintf(buf,"bsub < g%03d >> EVOS/G%03d/M%03d/jobid",p,T->n,p);
  system(buf);  

}
//==================================================================
//     Relax Tribe
//==================================================================
void RELX_INT(Tribe *T)
{
  int  p;
  char buf[200];

  EXIT_TR(T);
  for(p=T->N;p<T->N*2;p++)
  {
    sprintf(buf,"mkdir EVOS/G%03d/M%03d",T->n,p);
    system(buf);
    sprintf(buf,"EVOS/G%03d/M%03d/POSCAR.0",T->n,p);
    SAVE_CELL(&T->C[p],buf,0);

    Copy_C(&T->C[p],&T->C[2*T->N+1]);    

    CELL_RELX(RR,PP,WW,&T->C[2*T->N+1],LL);
//    CELL_PHON(RR,PP,WW,&T->C[2*T->N+1],LL);

    Copy_C(&T->C[2*T->N+1],&T->C[p]);

    sprintf(buf,"EVOS/G%03d/M%03d/CONTCAR.0",T->n,p);
    SAVE_CELL(&T->C[2*T->N+1],buf,0);
    sprintf(buf,"mv OUTCAR EVOS/G%03d/M%03d/OUTCAR.0",T->n,p);
    system(buf);
    sprintf(buf,"mv OSZICAR EVOS/G%03d/M%03d/OSZICAR.0",T->n,p);
    system(buf);
  }
}
//==================================================================
//     Relax Tribe
//==================================================================
void RELX_TR(Tribe *T)
{
  int i,p,Q,t1,t2,*I,m;
  char buf[400],s[200];
  FILE *in;

  I = make_i1D(2*T->N);

  if( T->JOBT!= 14 )
  {
    system("sleep 3");      // to allow completion of PLOT_TR output
    sprintf(buf,"EVOS/G%03d",T->n);
    if(chdir(buf)==0)
    {
      chdir("../../");
      sprintf(buf,"Directory EVOS/G%03d already exists!\n",T->n);
      Print_LOG(buf);
      fprintf(stderr,"Directory EVOS/G%03d already exists!\n",T->n);
      exit(1);
    }
  }
  sprintf(buf,"mkdir -p EVOS/G%03d",T->n);
  system(buf);
  if( T->MODE<=0 && T->n>0)
  {
    for(p=0;p<T->N;p++)
    {
      sprintf(buf,"mkdir EVOS/G%03d/M%03d",T->n,p);
      system(buf);
      sprintf(buf,"EVOS/G%03d/M%03d/CONTCAR.0",T->n,p);
      SAVE_CELL(&T->C[p],buf,0);
      sprintf(buf,"cp EVOS/G%03d/M%03d/OSZICAR.0 EVOS/G%03d/M%03d",T->n-1,T->S[p],T->n,p);
      system(buf);
      sprintf(buf,"cp EVOS/G%03d/M%03d/OUTCAR.0 EVOS/G%03d/M%03d",T->n-1,T->S[p],T->n,p);
      system(buf);
    }
  }
  //===== Relax everything internally =====
  if(T->CODE==0)
    RELX_INT(T);
  //===== Relax everything externally =====
  else
  {
    for(p=T->N;p<2*T->N;p++)
      I[p] = 0;

    if( T->MODE <= 0 )
      for(p=T->N;p<2*T->N;p++)
	QSUB_TR(T,p);
    if( T->MODE < 0 )
    {
      printf("Generated %3d structures in EVOS/G%03d\n",T->N,T->n);
      exit(0);
    }

    p = T->N;
    while(p<2*T->N)
    {
      system("sleep 5");
      system("date -d now +%s > stamp");
      t1 = TIME("stamp");
      
      EXIT_TR(T);               //Check setup and EXIT if needed
      for(p=T->N;p<2*T->N;p++)
	if(I[p]==0)
	{
	  sprintf(buf,"EVOS/G%03d/M%03d/OUTCAR.0",T->n,p);
	  Q = Check_OUTCAR(&T->C[p],buf);
	  sprintf(buf,"EVOS/G%03d/M%03d/stamp",T->n,p);
	  t2 = TIME(buf);
	  if( Q==1 )
	  {
	    I[p] = 1;
	    sprintf(buf,"EVOS/G%03d/M%03d/CONTCAR.0",T->n,p);
	    READ_CELL(&T->C[p],buf);
	  }
	  if( T->MODE > 0 )
	    t1 = t2;
	  if( Q==0 || (Q<0 && t1>=0 && t2>=0 && ( (int)(t1-t2) >T->time ) ) || (Q==1 && CHCK_Rm(&T->C[p],T->Rm,0.7)==0) )
	  {
	    sprintf(buf,"Failed to relax structure %d %d %8d sec %8d %8d % d\n",T->n,p,t1-t2,t1,t2,Q);
	    fprintf(stderr,"Failed to relax structure %d %d %8d sec %8d %8d % d\n",T->n,p,t1-t2,t1,t2,Q);
	    Print_LOG(buf);
	    I[p] = -1;
	    if( Q<0 )
	    {
	      sprintf(buf,"EVOS/G%03d/M%03d/jobid",T->n,p);
	      if( (in = fopen(buf,"r"))!=0 )
	      {
		fgets(buf,200,in);
		for(i=0;i<strlen(buf)&&(buf[i]<48||buf[i]>57);i++);
		if( i<strlen(buf) )
		{
		  printf("%3d %3d %s %s\n",T->n,p,buf,buf+i);
		  sscanf(buf+i,"%s",s);
		  if( T->QT==0 )
		    sprintf(buf,"qdel %s",s);
		  if( T->QT==1 )
		    sprintf(buf,"scancel %s",s);
		  if( T->QT==2 )
		    sprintf(buf,"bkill %s",s);
		  system(buf);
		  sprintf(buf,"deleting job %s\n",s);
		  Print_LOG(buf);
		}
		fclose(in);
	      }
	    }
	  }
	}
      for(p=T->N;p<2*T->N;p++)
	if(I[p]==0)
	  break;
    }
    for(m=0,p=T->N;p<2*T->N;p++)
      if( I[p]==1 )
	m++;
    if( m==0 )
    {
      sprintf(buf,"All relaxations failed\n");
      fprintf(stderr,"All relaxations failed\n");
      Print_LOG(buf);
      exit(1);
    }
    else
    {
      sprintf(buf,"Iteration %3d: Successfully relaxed %3d members\n",T->n,m);
      Print_LOG(buf);
    }

    for(p=T->N;p<2*T->N;p++)
      if(I[p]==-1)
      {
	//=== if the relaxation failed in generation > 0 copy from the previous generation ===
	i = p-T->N;
	//=== if the relaxation failed in generation = 0 copy from the     same generation ===
	if(T->n==0)
	  for(m=1;m<T->N;m++)
	  {
	    i = T->N + (p+m)%T->N;
	    if(I[i]==1)
	      break;
	  }
	sprintf(buf,"rm -r EVOS/G%03d/M%03d/*",T->n,p);
	system(buf);
	sprintf(buf,"cp    EVOS/G%03d/M%03d/* EVOS/G%03d/M%03d",T->n,i,T->n,p);
	system(buf);
      }
  }

  for(p=T->N;p<2*T->N;p++)
  {
    sprintf(buf,"EVOS/G%03d/M%03d/CONTCAR.0",T->n,p);
    READ_CELL(&T->C[p],buf);
    sprintf(buf,"EVOS/G%03d/M%03d/OSZICAR.0",T->n,p);
    T->C[p].P = Read_OSZI(buf);
    if(T->ND==0)
      NANO_ROT(&T->C[p],0);
  }
  // -----------------------------------------------------------------------
  // only T->N random structures are relaxed at n=0
  // the duplicates will be eliminated anyway
  // -----------------------------------------------------------------------
  if(T->n==0)
    for(p=0;p<T->N;p++)
    {
      sprintf(buf,"mkdir EVOS/G%03d/M%03d",T->n,p);
      system(buf);
      sprintf(buf,"cp EVOS/G%03d/M%03d/CONTCAR.0 EVOS/G%03d/M%03d",T->n,T->N+p,T->n,p);
      system(buf);
      sprintf(buf,"cp EVOS/G%03d/M%03d/OSZICAR.0 EVOS/G%03d/M%03d",T->n,T->N+p,T->n,p);
      system(buf);
      sprintf(buf,"cp EVOS/G%03d/M%03d/OUTCAR.0  EVOS/G%03d/M%03d",T->n,T->N+p,T->n,p);
      system(buf);
      Copy_C(&T->C[T->N+p],&T->C[p]);
    }

  if(T->CODE>1)
  {
    for(p=T->N;p<2*T->N;p++)
    {
      sprintf(buf,"rm -f g%03d*",p);
      system(buf);
    }
    sprintf(buf,"rm -f EVOS/G%03d/M*/OUT* EVOS/G%03d/M*/POTCAR* EVOS/G%03d/M*/INCAR* EVOS/G%03d/M*/KPOINTS*",T->n-2,T->n-2,T->n-2,T->n-2);
    if(T->CODE==1)
      system(buf);
  }
  sprintf(buf,"cp EVOS/G%03d/M000/CONTCAR.0 poscar",T->n);

  free_i1D(I);

}
//==================================================================
//     Evolve Tribe
//==================================================================
void EVLV_TR(Tribe *T, int J)
{
  char buf[200];
  int  N;

  N = 2*T->N;

  if(T->n>=0)
  {
    RANK_TR(T);
    SLCT_TR(T);
    if(J)
      PLOT_TR(T);
  }

  if(J==0||T->n<T->NI-1)
  {
    //=====  nanoparticle evolution operations  =====
    if(T->ND==0)
    {
      NANO_TETR(T, 0 ,&T->C[N],-1);  // tetris-type generation of random NPs
      NANO_PLNT(T, 1 ,-1);           // seeding NPs with INI/POSCAR000
      NANO_PACK(T, 2 ,&T->C[N],-1);  // seeding NPs with bcc, fcc, hcp
      NANO_BLOB(T, 3 ,&T->C[N],-1);  // random generation of NPs
      NANO_MATE(T, 4 );              // crossover using plane cuts
      NANO_SWAP(T, 5 );              // crossover using core-shells
      NANO_RUBE(T, 6 );              // Rubik's cube operation
      NANO_SYMM(T, 7 );              // symmetrization using reflection
      NANO_SYMM(T, 8 );              // symmetrization using inversion
      NANO_CHOP(T, 9 ,&T->C[N],-1);  // chopping for creating facets
      NANO_MUTE(T,10 );              // distortion using random shifts
    }
    //=====  bulk structure evolution operations  =====
    else
    {
      BULK_MATE(T, 4 );              // crossover using plance cuts
      BULK_MUTE(T,10 );              // distortion using random shifts
    }
  }
  if(T->JOBT==11)
  {
    sprintf(buf,"Soft exit at iteration %d\n",T->n);
    Print_LOG(buf);
  }
}
//==================================================================
void ANA_EVOS(Tribe *T, Cell *C, Cell *D)
{
  int i,n,N,m,M,g,G,*I,J;
  double E0,*E,*V,tol;
  char s[400],buf[400],dir[200],run[10][200];
  FILE *in,*out;

  tol = 0.1;

  J = 0;
  getcwd(dir,200);
  sprintf(run[0],"EVOS");
  while(chdir(run[J])==0)
  {
    chdir(dir);
    J++;
    sprintf(run[J],"POOL/RUN%d",J-1);    
  }
  if(J==0)
  {
    printf("Directory EVOS is not found\n");
    exit(0);
  }
  J--;

  sprintf(s,"find %s -name CONTCAR.0 | sort",run[J]);
  in = popen(s,"r");
  M = 0;
  while(fgets(buf,200,in))
    M++;
  pclose(in);

  E = make_d1D(M);
  V = make_d1D(M);
  I = make_i1D(M);
  char DIR[M][200];

  //===== for analysis of EVOS get the lowest enthalpy from ebest.dat =====
  E0 = 1e12;
  if(J==0)
  {
    sprintf(s,"tail -1 ebest.dat");
    in = popen(s,"r");
    fscanf(in,"%d %lf\n",&i,&E0);
    pclose(in);
    sprintf(s,"mkdir -p POOL");
    system(s);
  }
  else
  {
    sprintf(s,"find %s -name CONTCAR.0 | sort",run[J]);
    in = popen(s,"r");
    for(m=0;m<M;m++)
    {
      fscanf(in,"%s\n",dir);
      dir[strlen(dir)-strlen("CONTCAR.0")-1] = 0;
      sprintf(s,"%s/CONTCAR.0",dir);
      READ_CELL(C,s);
      sprintf(s,"%s/OSZICAR.0",dir);
      E[m] = Read_OSZI(s)/(double)C->N;
      sprintf(s,"%s/OUTCAR.0",dir);
      if( fabs((E[m]-1.0)*(double)C->N)>1e-14 && E[m]<E0 && Check_OUTCAR(C,s)==1 )
	E0 = E[m];
    }
    pclose(in);
  }
  sprintf(s,"find %s -name CONTCAR.0 | sort",run[J]);
  in = popen(s,"r");
  for(m=G=0;m<M;m++)
  {
    fscanf(in,"%s\n",dir);
    dir[strlen(dir)-strlen("CONTCAR.0")-1] = 0;
    sprintf(s,"%s/CONTCAR.0",dir);
    READ_CELL(C,s);
    sprintf(s,"%s/OSZICAR.0",dir);
    E[G] = Read_OSZI(s)/(double)C->N;
    V[G] = CELL_VOL(C) /(double)C->N;
    sprintf(s,"%s/OUTCAR.0",dir);
    if( fabs((E[G]-1.0)*(double)C->N)>1e-14 && E[G]<(E0+T->DE) && Check_OUTCAR(C,s)==1 )
      sprintf(DIR[G++],"%s",dir);
  }
  pclose(in);

  Sort(E,I,G);

  system("rm -f POOL/POSCAR*");

  for(N=0,g=G-1;g>=0;g--)
  {
    sprintf(s,"%s/CONTCAR.0",DIR[I[G-g-1]]);
    READ_CELL(C,s);
    LIST(C,0);
    RDF(C,1);
    for(n=0;n<N;n++)
    {
      sprintf(buf,"POOL/POSCAR%03d",n);
      READ_CELL(D,buf);
      LIST(D,0);
      RDF(D,1);
      if(CxC(C,D)>T->CUT)
	break;
    }
    if(n==N)
    {
      sprintf(buf,"POOL/POSCAR%03d",N);
      sprintf(C->TAG,"%s % lf",DIR[I[G-g-1]],E[I[G-g-1]]);
      if(C->ND==0)
	NANO_ROT(C,0);
      SAVE_CELL(C,buf,0);
      N++;
    }
  }        

  sprintf(buf,"mkdir -p POOL/SET%d; cp POOL/POSCAR* POOL/SET%d",J,J);
  system(buf);

  sprintf(buf,"ana%d.dat",J);
  out = fopen(buf,"w");
  for(n=0;n<N;n++)
  {
    sprintf(buf,"POOL/POSCAR%03d",n);
    READ_CELL(C,buf);
    READ_CELL(D,buf);
    sscanf(C->TAG,"%s %lf",buf,&E[n]);
    V[n] = CELL_VOL(C)/(double)C->N;
    if(C->ND==0)
    {
      CENTER(C,0.5);
      CENTER(D,0.5);
    }
    FIND_WYC(C,D,tol,1);
    READ_CIF(C,"str.cif",tol,C->NM,"POSCAR");
    FIND_PRS(C,D,tol);
    printf("%4d %3d %s%d % 12.6lf % 12.6lf % 12.6lf\n",n,C->SGN,C->PRS,C->N,V[n],E[n],E[n]-E0);
    fprintf(out,"%4d %3d %s%d % 12.6lf % 12.6lf % 12.6lf\n",n,C->SGN,C->PRS,C->N,V[n],E[n],E[n]-E0);
  }
  fclose(out);
  printf("Found %6d candidates in %s with lowest enthalpy % 12.6lf\n",M,run[J],E0);
  printf("Found %6d structures in %s with enthalpy range  % 12.6lf % 12.6lf\n",G,run[J],E0,E0+T->DE);
  printf("Found %6d structures in %s with enthalpy range  % 12.6lf % 12.6lf   different by CxC % lf\n",N,run[J],E0,E0+T->DE,T->CUT);

  sprintf(buf,"POOL/INI%d",J);
  if(chdir(buf)==0)
  {
    chdir("..");
    sprintf(buf,"mkdir -p RUN%d",J);
    system(buf);    
    for(n=0;n<N;n++)
    {
      sprintf(buf,"cp -r INI%d RUN%d/M%03d; cp POSCAR%03d RUN%d/M%03d/POSCAR",J,J,n,n,J,n);
      system(buf);
      sprintf(buf,"cd RUN%d/M%03d; mv sub g%03d; sbatch g%03d",J,n,n,n);
      system(buf);
    }
  }

  free_d1D(E);
  free_d1D(V);
  free_i1D(I);
}
//==================================================================
// build the tribe
//==================================================================
void INIT_EVOS(Tribe *T, Cell *C)
{
  int k,n,i,j;
  char buf[200];

  T->C  = (Cell *)malloc((2*T->N+2)*sizeof(Cell));

  C->N = C->A = C->A*2;
  Build_Cell(C,0);

  C->XT = 1;
  for(j=0;j<T->NSPC;j++)
    C->SPCZ[j] = T->SPCZ[j];
  for(j=i=0;j<C->NSPC;j++)
    for(k=0;k<C->SPCN[j];k++)
      C->ATMZ[i++] = C->SPCZ[j];

  C->NSPC = T->NSPC;
  C->N = C->A = C->A/2;

  //===== conventional cells in EVOS analysis could be 4 times larger =====
  if(T->JOBT==13)
    C->N = C->A = C->A*4;

  for(n=0;n<2*T->N+2;n++)
  {
    T->C[n].A    = C->A;  // for NANO_ROT we need +2 atoms
    T->C[n].N    = C->N;
    T->C[n].NM   = C->NM;
    T->C[n].ND   = C->ND;
    T->C[n].Rmax = C->Rmax;
    T->C[n].Rmin = C->Rmin;
    T->C[n].DR   = C->DR;
    T->C[n].JOBT = C->JOBT;
    T->C[n].XT   = 1;
    T->C[n].NSPC = C->NSPC;
    T->C[n].Rc   = C->Rc = C->Rmax;
    T->C[n].rc   = C->rc = C->Rmax;
    T->C[n].MODT = C->MODT;
    T->C[n].NP   = C->NP;
    T->C[n].NB   = C->NB;

    for(k=0;k<C->NSPC;k++)
    {
      T->C[n].SPCZ[k] = T->SPCZ[k];
      T->C[n].SPCN[k] = T->SPCN[k];
    }
    sprintf(T->C[n].WDIR,"%s",C->WDIR);
    Build_Cell(&T->C[n],0);
  }

  if(T->JOBT==13)
    C->N = C->A = C->A/4;
  else
    for(n=0;n<2*T->N+2;n++)
      Copy_C(C,&T->C[n]);

  if(T->JOBT!=10&&T->JOBT!=14)
    return;

  if(T->CODE==0)
    READ_POT(C,".");

  T->C[2*T->N+1] = *C;

  //===== for planting seed structures =====
  T->pos = 0;
  if(T->TES[1]>0)
    T->pos = 1;
  for(;T->pos<T->N;T->pos++)
  {
    sprintf(buf,"INI/POSCAR%03d",T->pos);
    if( READ_CELL(&T->C[2*T->N+1],buf)==0 )
      break;
  }
  sprintf(buf,"PLNT will plant %3d from %3d structures\n",T->TES[1],T->pos);
  Print_LOG(buf);

  T->I  = make_i2D(T->NI,2*T->N);
  T->S  = make_i1D(2*T->N);
  T->E  = make_d1D(2*T->N+1);
  T->P1 = make_d1D(2*T->N);
  T->P2 = make_d1D(2*T->N);
  T->T  = make_d2D(T->NI,2*T->N);
  T->f  = make_d1D(2*T->N+1);  
  for(n=0;n<2*T->N+2;n++)
    T->C[n].P = 1.0e9;
  
  return;
}
//==================================================================
// prepare EVOS for manual execution
//==================================================================
void EVOS_PREP(Tribe *T)
{
  char buf[200];
  FILE *in;

  // T->MODE = -1: generate population
  // T->MODE =  1: evaluate population

  sprintf(buf,"grep -c -s . ebest.dat");
  in = popen(buf,"r");
  fscanf(in,"%d",&T->n);
  pclose(in);
  sprintf(buf,"ls -l EVOS/G*/M%03d/OUTCAR.0 2>/dev/null | wc -l ",T->N);
  in = popen(buf,"r");
  fscanf(in,"%d",&T->MODE);
  pclose(in);
  sprintf(buf,"EVOS/G%03d/M%03d",T->n,T->N);
  if( T->MODE==T->n && chdir(buf)==0 )
  {
    chdir("../../../");
    printf("Population EVOS/G%03d already exists and should be run first.\n",T->n);
    exit(0);
  }
  T->MODE = 2*(T->MODE - T->n)-1;
  if( abs(T->MODE)!= 1 )
  {
    printf("The number of completed generations in EVOS is inconsistent with the number of analyzed generations in ebest.dat\n");
    exit(0);
  }
}
//==================================================================
void EVOS_MAIN(Tribe *T, ANN *R, PRS *P, Cell *C)
{
  char buf[200];

  printf("|                         Evolutionary search                         |\n");
  printf("=======================================================================\n\n");

  RR = R;
  PP = P;

  INIT_EVOS(T,C);
  R->MODT = C->MODT;

  if(T->seed==0)                          // Check if seed is by hand or auto if so:
    T->seed=time(NULL);                   // SET TO TIME
  sprintf(buf,"seed:  %ld\n\n",T->seed);
  Print_LOG(buf);
  PlantSeeds(T->seed);                    // Initialize random number generator

  //=====  Analyze completed EVOS  =====
  if(T->JOBT==13)
  {    
    ANA_EVOS(T,&T->C[0],&T->C[1]);
    exit(0);
  }

  //=====  Manual  EVOS  execution =====
  if(T->JOBT==14)
  {
    EVOS_PREP(T);
    INIT_TR(T);
    RELX_TR(T);
    EVLV_TR(T,1);
    T->n++;
    T->MODE *= -1;
    RELX_TR(T);
    exit(0);
  }

  //===== Automated EVOS execution =====
  T->MODE = 0;                            // automated EVOS with JOBT=10
  INIT_TR(T);                             // Populate the tribe T
  if(T->n!=0)                             // If this is a continuation run evolve first
    EVLV_TR(T,0);
  for(;T->n<T->NI&&T->JOBT==10;T->n++)    // for each genertation relax then evolve
  {
    RELX_TR(T);
    EVLV_TR(T,1);
  }
  if(T->JOBT==12)
  {
    sprintf(buf,"Hard exit at iteration %3d \n",T->n);
    Print_LOG(buf);
    exit(0);
  }

}
//==================================================================
