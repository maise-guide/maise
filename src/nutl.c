#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cdef.h"
#include "ndef.h"
#include "edef.h"
#include "util.h"
#include "nprs.h"
#include "nmlp.h"
#include "nutl.h"

//==============================================================================
void READ_STRAT3(ANN *R)
{
  int   i,j,k,l,m,n;
  char  s[200],s1[200],s2[200];
  int   NG[20],NC[3];
  FILE  *in,*mlp;
  int   t,a;

  int   tag[10]={4,2,1,0,0,0,0,0,0,0};        //variable cmp
  int   x,z,ind;
  int   spc[3][2]={{0,1},{0,2},{1,2}};
  int   bin[3][2]={{0,1},{0,2},{1,2}};
  int   ter[3][3]={{0,1,3},{0,2,5},{3,4,5}};

  //===== reading basis parameters from the first model file in the output directroy =====
  for(i=0; i < 20;i++) 
    NG[i]=0;
  atom_symb(R->SPCZ[0],s1);
  atom_symb(R->SPCZ[1],s2);
  sprintf(s,"%s/%s%s.dat",R->otpt,s1,s2);
  if( (in=fopen(s,"r")) == 0 )
  {
    printf("ERROR opening %s file to read basis parameters\n",s);         
    exit(1);
  }
  while( fgets(s,200,in) )
    if( strncmp(s,"  B2A",5) == 0 )
      break;
  for(i=0; i < 4;i++) 
    fgets(s,200,in);
  for(i=0; i < R->NSYM;i++)
  {
    fgets(s,200,in);
    sscanf(s,"%d %d %d %d %d %d %d %d %d\n",&t,&a,&t,&t,&t,&t,&t,&t,&t);
    NG[a]++;
    if( a > 2 ) 
      NG[0]++;                       // total number of triplet (NG4) functions
  }
  fclose(in);

  NC[0]=R->NSPC;                     // number of components/sym_function: pair
  for(t=R->NSPC,a=0;t > 0;t--) 
    a+=t;

  NC[1]=a;                           // number of components/sym_function: triplet

  R->nw = R->NSPC*(NG[0])*R->NU[1];  // total number of adjustable weights

  R->mask = make_i2D(R->NSPC+1,R->D);

  for(i=0; i < R->D;i++)                // mark all as NON-adjustable=0 at first
    for(j=0; j < R->NSPC;j++) 
      R->mask[j][i] = 0;

  for(i=0,l=1; i < R->NSPC-1;i++)
    for(j=i+1; j < R->NSPC;j++)
    {
      atom_symb(R->SPCZ[i],s );
      atom_symb(R->SPCZ[j],s2);
      sprintf(R->file_list[l],"%s/%s%s.dat",R->otpt,s,s2);
      l++;
    }

  for(i=0; i < R->NSPC;i++)
    for(n=l=0; n < R->D;n++) 
      if( n == NG[2]*NC[0]+l*NC[1]+tag[i])
      {
	R->mask[i][n] = 1;
	l++;
      }

  //===== reading 'model' file parameters =====
  for(ind=0; ind < 3;ind++)
  {
    mlp = fopen(R->file_list[ind+1],"r");
    while(fgets(s,200,mlp))
      if(strncmp(s,"        minimum",15)==0)
        break;
    for(z=0;z<2;z++)
    {
      for(i=0;i<NG[2];i++)
      {
        j=i*3;
        for(x=0;x<2;x++)
        {
          fgets(s,200,mlp);
          sscanf(s,"%lf %lf %lf\n",&R->Rmin[spc[ind][z]][j+bin[ind][x]],&R->Rmax[spc[ind][z]][j+bin[ind][x]],&R->DR[spc[ind][z]][j+bin[ind][x]]);
        }
      }
      for(i=0;i<NG[0];i++)
      {
        j=NG[2]*3+i*6;
        for(x=0;x<3;x++)
        {
          fgets(s,200,mlp);
          sscanf(s,"%lf %lf %lf\n",&R->Rmin[spc[ind][z]][j+ter[ind][x]],&R->Rmax[spc[ind][z]][j+ter[ind][x]],&R->DR[spc[ind][z]][j+ter[ind][x]]);
        }
      }
    }
    

    while(fgets(s,200,mlp))
      if(strncmp(s,"|                        neural",31)==0)
        break;
    fgets(s,200,mlp);
    
    for(z=0;z < 2;z++)
    {
      for(k=i=0; k < R->NL-1;k++)
	for(m=0; m < R->NU[k+1];m++)
	{
	  fgets(s,200,mlp);
	  sscanf(s,"%lf\n",&R->B[spc[ind][z]][k][m]);
	}

      for(m=k=0; m < R->NU[k+1];m++)
      {
	for(i=0; i < NG[2];i++)
	{
	  j=i*3;
	  for(x=0; x < 2;x++)
	  {
	    fgets(s,200,mlp); 
	    sscanf(s,"%lf\n",&R->W[spc[ind][z]][k][m][j+bin[ind][x]]);
	  }
	}

	for(i=0; i < NG[0];i++)
	{
	  j=NG[2]*3+i*6;
	  for(x=0; x < 3;x++)
	  {
	    fgets(s,200,mlp); 
	    sscanf(s,"%lf\n",&R->W[spc[ind][z]][k][m][j+ter[ind][x]]);
	  }
	}
      }

      for(k=1; k < R->NL-1;k++)
	for(m=0; m < R->NU[k+1];m++)
	  for(n=0; n < R->NU[k];n++)
	  {
	    fgets(s,200,mlp);
	    sscanf(s,"%lf\n",&R->W[spc[ind][z]][k][m][n]);
	  }
    }
    fclose(mlp);
  }
}
//==============================================================================
//                                            
//==============================================================================
void READ_STRAT2(ANN *R)
{
  int    i,spc,j,k,l,m,n;
  char   s[200],s2[200];
  FILE   *in,*mlp;
  int    t,a;
  int    NG[20],NC[3];
  int    tag[10]={0,2,0,0,0,0,0,0,0,0};           //fixed cmp
  double **w;

  //===== reading basis parameters from the first model file in the output directroy =====
  for(i=0; i < 20;i++) 
    NG[i]=0;
  atom_symb(R->SPCZ[0],s2);
  sprintf(s,"%s/%s.dat",R->otpt,s2);
  if( (in=fopen(s,"r")) == 0 )
  {
    printf("ERROR opening %s file to read basis parameters\n",s);         
    exit(1);
  }
  while( fgets(s,200,in) )
    if( strncmp(s,"  B2A",5) == 0 )
      break;

  for(i=0; i < 4;i++)
    fgets(s,200,in);

  for(i=0; i < R->NSYM;i++)
  {
    fgets(s,200,in);
    sscanf(s,"%d %d %d %d %d %d %d %d %d\n",&t,&a,&t,&t,&t,&t,&t,&t,&t);
    NG[a]++;
    if( a > 2) 
      NG[0]++;          // total number of triplet (NG4) functions, in case of combinations                             
  }
  fclose(in);



  NC[0]=R->NSPC;        // number of components/sym_function: pair
  for(t=R->NSPC,a=0; t > 0;t--) 
    a+=t;

  NC[1]=a;              // number of components/sym_function: triplet

  R->nw = R->NSPC*(NG[2]+2*NG[0])*R->NU[1];            // total number of adjustable weights  
  w     = make_d2D(R->NSPC+1,((R->NW-R->nw)/R->NSPC)); // vectors containing NON-adjustable weights per species  

  R->mask = make_i2D(R->NSPC+1,R->D);
  //creating mask list: mark all as adjustable=1
  for(i=0; i < R->D;i++) 
    for(j=0; j < R->NSPC;j++) 
      R->mask[j][i]=1;

  for(i=0,l=1; i < R->NSPC;i++)
  {
    atom_symb(R->SPCZ[i],s);
    sprintf(R->file_list[l],"%s/%s.dat",R->otpt,s);
    l++;
  }

  for(spc=0; spc < R->NSPC;spc++)
    for(n=j=l=0; n < R->D;n++) 
    {
      if( j<NG[2] && n == (j*NC[0])+spc ) 
      { 
	R->mask[spc][n]=0;
	j++;
      }
      else if( n == NG[2]*NC[0]+l*NC[1]+tag[spc] )
      {
	R->mask[spc][n]=0;
	l++;
      }
    }

  //===== reading 'model' file parameters =====
  for(spc=0; spc < R->NSPC;spc++)
  {
    mlp = fopen(R->file_list[spc+1],"r");     //file_list actual index starts from spc=1 to spc=NSPC
    
    while(fgets(s,200,mlp))
      if(strncmp(s,"|  number of weights",20)==0)
        break;
    sscanf(s,"%s %s %s %s %s %d %s \n",s2,s2,s2,s2,s2,&t,s2);
    if( t != (R->NW-R->nw)/R->NSPC )
    {
      fprintf(stderr,"Error: not consistent elemental and multi-component parameters\n");
      exit(1);
    }
    
    while(fgets(s,200,mlp))
      if(strncmp(s,"        minimum",15)==0)
        break;
    
    for(i=0; i < R->D;i++)
      if( R->mask[spc][i] == 0 )
      {
	fgets(s,200,mlp);
	sscanf(s,"%lf %lf %lf\n",&R->Rmin[spc][i],&R->Rmax[spc][i],&R->DR[spc][i]);
      }      

    while(fgets(s,200,mlp))
      if(strncmp(s,"|                        neural",31)==0)
        break;
    fgets(s,200,mlp);
 
    for(i=0; i < t;i++)
    {
      fgets(s,200,mlp);
      sscanf(s,"%lf\n",&w[spc][i]);
    }
    fclose(mlp);
  }

  for(spc=0; spc < R->NSPC;spc++)
  {      
    i=0;           
    for(k=0; k < R->NL-1;k++)
      for(m=0; m < R->NU[k+1];m++)
	R->B[spc][k][m] = w[spc][i++];

    k=0;
    for(m=0; m < R->NU[k+1];m++)
      for(n=0; n < R->NU[k];n++)
	if( R->mask[spc][n] == 0 )
	  R->W[spc][k][m][n] = w[spc][i++]; 

    for(k=1; k <R->NL-1;k++)
      for(m=0;m<R->NU[k+1];m++)
	for(n=0;n<R->NU[k];n++)
	  R->W[spc][k][m][n] = w[spc][i++];
  }
}
//==============================================================================
//                                            
//==============================================================================
void ANA_STR(ANN *R)
{
  char   s[200];
  FILE   *in;
  int    N,MAX,COMP,ns,spc[10],t[10];
  int    i,n,k,tag=0;

  sprintf(s,"%s/stamp.dat",R->data);
  if( (in=fopen(s,"r")) == 0 )
  {
    fprintf(stderr,"%s not found. Please parse the data in %s first\n",s,R->data);
    exit(0);
  }

  fgets(s,200,in);
  fgets(s,200,in);
  fgets(s,200,in);
  fscanf(in,"%s %d\n",s,&N);//# of parsed structures
  fscanf(in,"%s %d\n",s,&COMP);//# of components in stamp file
  fgets(s,200,in);
  fscanf(in,"%s %d\n",s,&ns);//# of species
  fgets(s,200,in);
  for(i=0,n=0,k=4; i < ns;i++,k+=n) sscanf(s+k,"%d%n",&t[i],&n);
  for(i=0; i < ns;i++) spc[i]=t[i];
  fscanf(in,"%s %d\n",s,&MAX);
  fclose(in);

  if( ns != R->NSPC )
  {fprintf(stderr,"Error: inconsistent number of species in stamp.dat and setup files!\n");exit(1);}

  for(i=0; i < ns;i++) 
    if( spc[i] != R->SPCZ[i] ) 
      tag++;

  if( tag != 0 ) 
  {
    fprintf(stderr,"Error: inconsistent types of species in stamp.dat and setup files!\n");
    exit(1);
  }

  if( COMP != R->D ) 
  {
    fprintf(stderr,"Error: inconsistent number of components in stamp.dat and setup files!\n");
    exit(1);
  }

  if( R->N < 0 )
    R->N  = -(int)((double)N*(double)R->N/100.0);

  if( R->TN <0 )
    R->TN = -(int)((double)N*(double)R->TN/100.0);

  if( N < (R->N+R->TN) ) 
  {
    fprintf(stderr,"Error: insufficient number of available structures!\n");
    exit(1);
  }

  R->STR = R->N+R->TN;
  R->A   = MAX;
}
//==============================================================================
//
//==============================================================================
void RNDNM2(int *NM, int N, int sd)
{
  int    i,j,*NN;

  NN = make_i1D(N);

  if( sd == 1 ) srand (time(NULL));

  for(i=0; i < N;i++)
    NN[i] = i;

  for(i=0; i < N;i++)
  {
    j = floor(((double)rand()/(double)RAND_MAX)*(double)(N-i));
    NM[i] = NN[j];
    NN[j] = NN[N-i-1];
  }

  if(0)
    for(i=0; i < N;i++)
      printf("%d\n",NM[i]);

  free_i1D(NN);
}
//==============================================================================
// Reading cpu time                                                                                      
//==============================================================================
double cpu_time( )
{
  double value;
  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//==============================================================================
// Mark atoms for force trainig
//==============================================================================
void MARK_CL(Cell *C, int M, long seed)
{
  int i,j,m,k;

  if( seed == 0 ) seed=time(NULL);
  PlantSeeds(seed);
  if( M == C->N )
  {
    for(i=0; i < C->N;i++)
      C->FRC[i] = 1;

    return;
  }

  for(i=0; i < C->N;i++)
    C->FRC[i] = 0;

  for(m=i=0; m < M;m++)
  {
    j = (int)(Random()*(double)(C->N-m));
    for(k=0; k < j;k++)
      if( C->FRC[(i+k+1)%C->N] == 1 )
	j++;

    i = (i+k)%C->N;
    C->FRC[i] = 1;
    while( C->FRC[i%C->N] == 1 )
      i++;

    i = i%C->N;
  }

  for(i=k=0; i < C->N;k+=C->FRC[i],i++);

  if( k != M )
  {
    fprintf(stderr,"MARK_CL FAIL\n");
    exit(1);
  }
}
//==============================================================================
// Building ANN object using configuration
//==============================================================================
void Build_ANN(ANN *R)
{
  int    n;

  if( (R->JOBT/10 == 3 ) || (R->JOBT/10 == 4) )
  {
    R->E    = make_d1D(R->STR);            // total energies as target values energy
    R->Fc   = make_d1D(R->STR);            // force
    R->train= make_i1D(R->STR);            // Holds file index of structures to be used for training
    R->test = make_i1D(R->STR);            // Holds file index of structures to be used for testing
    R->TEeval=make_d1D(R->STR);            // Holds evaluated energies for test structures
  }

  R->R1   = make_d1D(R->D+6);              // input vector
  R->G    = make_d1D(R->D);                // weight of input variables
  R->Rmin = make_d2D(R->NSPC,R->D);
  R->Rmax = make_d2D(R->NSPC,R->D);
  R->DR   = make_d2D(R->NSPC,R->D);

  //====== MLP elements ====== 
  if( R->MODT == 1 )
  {
    for(n=R->N0=0; n < R->NL;n++)
      if( R->N0 < R->NU[n] )
	R->N0 = R->NU[n];

    for(n=1,R->N1=0; n < R->NL;n++)
      if( R->N1 < R->NU[n] )
	R->N1 = R->NU[n];

    for(n=R->NW=0; n < R->NL-1;n++)
      R->NW += (R->NU[n]+1)*R->NU[n+1];

    R->NW*=R->NSPC;

    R->W  = make_d4D(R->NSPC,R->NL,R->N1,R->N0);  // could be improved
    R->B  = make_d3D(R->NSPC+1,R->NL-1,R->N1);    
    R->Wp = make_d4D(R->NSPC,R->NL,R->N1,R->N0);
    R->Bp = make_d3D(R->NSPC+1,R->NL-1,R->N1);
    R->e  = make_d4D(R->NSPC,R->A,R->NL,R->N0); // for E training - not needed for openmp
    R->d  = make_d4D(R->NSPC,R->A,R->NL,R->N0); // for E training or F/S 
    R->c  = make_d4D(R->NSPC,R->A,R->NL,R->N0); // for F/S training
  }
  //======  SC or GP elements ======
  if( R->MODT == 3 || R->MODT == 4 )
  {
    R->WW  = make_d2D(3,20);   // need just 6, to be fixed
    R->WWp = make_d2D(3,20);   // need just 6, to be fixed    
    R->We  = make_d2D(3,R->A); 
  }  

  return;
}
//==============================================================================
void Build_LNK(LNK *L, int N, int NM, int D, int EFS)
{
  int i;

  if( L->B == 1 )
    return;

  L->N   = N;
  L->ATMN = make_i1D(N);
  L->NM  = NM;
  L->DNn = make_i1D(N);
  L->MRK = make_i1D(N);
  L->Fi  = make_i1D(N);
  L->EA  = make_d1D(N);
  L->F   = make_d2D(N,3);  
  L->f   = make_d2D(N,3);
  L->S   = make_d1D(6);
  L->s   = make_d1D(6);
  L->Cn  = make_d2D(N,D);

  if( EFS > 0 )
  {
    L->Fn  = make_d4D(N,NM+1,D,3);
    L->Sn  = make_d3D(N,D,6);
    L->DNi = make_i2D(N,NM);
    L->DNs = make_i2D(N,NM);
    for(i=0;i<N;i++)
      L->MRK[i] = L->Fi[i] = L->DNn[i] = 0;
  }

  L->B = 1;
  L->p = 0.0;
  sprintf(L->path,"none");
}
//==============================================================================
void SAVE_ANN(ANN *R, double TIME)
{
  int    i,spc;
  char   s[200],b[200],t[200],e[200];
  double *V,sum;
  FILE*  out;
  time_t rawtime;
  struct tm * timeinfo;

  V = make_d1D(R->NW);
  W2V(R,V);
  for(i=0,sum=0.0; i < R->NW;i++)
    sum += fabs(V[i]);

  time(&rawtime);
  timeinfo = localtime ( &rawtime );
  sprintf(b,"-----------------------------------------------------------------------------\n");
  sprintf(e,"                                                                        |\n");

  sprintf(s,"%s/%s",R->otpt,"model");
  out=fopen(s,"w");

  fprintf(out,"%s|                    neural network general information                     |\n%s",b,b);
  sprintf(t," %014lX",(long)(sum*1e11));
  fprintf(out,"|  model unique ID      | %s%s",t,e+22+strlen(t));

  sprintf(t," %d",R->NSPC);
  fprintf(out,"|  number of species    | %s%s",t,e+22+strlen(t));
  for(i=0,t[0]=0; i < R->NSPC;i++)
    sprintf(t+strlen(t)," %2d",R->SPCZ[i]);
  fprintf(out,"|  species types        | %s%s",t,e+22+strlen(t));
  for(i=0,t[0]=0,sprintf(t+strlen(t)," "); i < R->NSPC;i++,sprintf(t+strlen(t)," "))
    atom_symb(R->SPCZ[i],t+strlen(t));
  fprintf(out,"|  species names        | %s%s",t,e+22+strlen(t));
  for(i=0,t[0]=0; i < R->NSPC;i++)
    sprintf(t+strlen(t),"%12.8lf",R->E0[R->SPCZ[i]]);
  fprintf(out,"|  species references   | %s%s",t,e+22+strlen(t));

  sprintf(t," %d",R->NL);  
  fprintf(out,"|  number of layers     | %s%s",t,e+22+strlen(t));
  for(i=0,t[0]=0; i < R->NL;i++)
    sprintf(t+strlen(t)," %d",R->NU[i]);

  fprintf(out,"|  architecture         | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->NW);
  fprintf(out,"|  number of weights    | %s%s",t,e+22+strlen(t));
  sprintf(t," doi");
  fprintf(out,"|  reference            | %s%s",t,e+22+strlen(t));

  fprintf(out,"%s|                            data information                               |\n%s",b,b);
  sprintf(t," %d",R->N);
  fprintf(out,"|  train energy data    | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->TN);
  fprintf(out,"|  test  Energy data    | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->NF);
  fprintf(out,"|  train force  data    | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->TNF);
  fprintf(out,"|  test  Force  data    | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->STR);
  fprintf(out,"|  total structures     | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->A);
  fprintf(out,"|  max atom number      | %s%s",t,e+22+strlen(t));
  sprintf(t," %lf   eV/atom",R->Edev);
  fprintf(out,"|  standard deviation   | %s%s",t,e+22+strlen(t));

  fprintf(out,"%s|                            training details                               |\n%s",b,b);
  sprintf(t," %s",asctime(timeinfo));
  t[strlen(t)-1] = 0;
  fprintf(out,"|  trained on           | %s%s",t,e+22+strlen(t));
  sprintf(t," %s",R->VER);
  fprintf(out,"|  maise version        | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->NP);
  fprintf(out,"|  parallelization      | %s%s",t,e+22+strlen(t));
  sprintf(t," %1.2lf sec",TIME);
  fprintf(out,"|  training time        | %s%s",t,e+22+strlen(t));
  sprintf(t," %1.1e",R->LREG);
  fprintf(out,"|  regularization       | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->ITER);
  fprintf(out,"|  number of epochs     | %s%s",t,e+22+strlen(t));

  fprintf(out,"%s|                              performance                                  |\n%s",b,b);
  sprintf(t," %lf  eV/atom",R->RE);
  fprintf(out,"|  train energy error   | %s%s",t,e+22+strlen(t));
  sprintf(t," %lf  eV/atom",R->EE);
  fprintf(out,"|  test  Energy error   | %s%s",t,e+22+strlen(t));
  sprintf(t," %lf  eV/Ang",R->RF);
  fprintf(out,"|  train force  error   | %s%s",t,e+22+strlen(t));
  sprintf(t," %lf  eV/Ang",R->EF);
  fprintf(out,"|  test  Force  error   | %s%s",t,e+22+strlen(t));

  fprintf(out,"%s|                          input component range                            |\n%s",b,b);
  fprintf(out,"        minimum                    maximum                     range\n");
  for(spc=0; spc < R->NSPC;spc++)
    for(i=0; i < R->D;i++) 
      fprintf(out,"% 21.16lf % 25.16lf  % 25.16lf\n",R->Rmin[spc][i],R->Rmax[spc][i],R->DR[spc][i]);

  fprintf(out,"%s|                        neural network weights                             |\n%s",b,b);

  for(i=0; i < R->NW;i++)
    fprintf(out,"% 21.16lf\n",V[i]);

  free_d1D(V);
  fprintf(out,"%s",b);
  fclose(out);
}
//==============================================================================
//
//==============================================================================
void READ_ANN(ANN *R)
{
  char   s[200],s1[200];
  int    i,n,spc,k;
  FILE   *in;
  double *V;

  sprintf(s,"model");
  if( R->JOBT/10 == 4 )
    sprintf(s,"%s/model",R->otpt);

  if( (in=fopen(s,"r")) == 0 )
  {
    fprintf(stderr,"ERROR opening %s\n",s);
    exit(1);
  }

  R->NSPC = 0;
  while( fgets(s,200,in) )
    if( strncmp(s,"|  model unique ID      |",25) == 0 )
    {
      sscanf(s+26,"%s" ,R->ID);
      break;
    }

  while( fgets(s,200,in) )  
    if( strncmp(s,"|  number of species    |",25) == 0 )   
    {
      sscanf(s+26,"%d" ,&R->NSPC); 
      break;
    }

  if( R->NSPC == 0 )
  {
    fprintf(stderr,"Error in reading number of species from %s\n",s);
    exit(1);
  }

  while( fgets(s,200,in) )
    if( strncmp(s,"|  species types        |",25) == 0 )
    {
      for(i=0,n=0,k=26; i < R->NSPC;i++,k+=n) 
	sscanf(s+k,"%d%n",&R->SPCZ[i],&n);
      break;
    }

  strcpy(R->compound,"");
  for(i=0; i < R->NSPC;i++)
  {
    atom_symb(R->SPCZ[i],s1);
    strcat(R->compound,s1);
  }

  while( fgets(s,200,in) )
    if( strncmp(s,"|  species references   |",25) == 0 )
    {
      for(i=0,n=0,k=26; i < R->NSPC;i++,k+=n)
        sscanf(s+k,"%lf%n",&R->E0[R->SPCZ[i]],&n);
      break;
    }

  while( fgets(s,200,in) )
    if( strncmp(s,"|  number of layers     |",25) == 0 )
    {
      sscanf(s+26,"%d",&R->NL);
      break;
    }

  while( fgets(s,200,in) )
    if( strncmp(s,"|  architecture         |",25) == 0 )
    {
      for(i=0,n=0,k=26; i < R->NL;i++,k+=n)
        sscanf(s+k,"%d%n",&R->NU[i],&n);
      break;
    }

  //Default type of the neurons in hiddern layers
  // =1 (tanh) and =0 (linear) for output layer
  if(R->NL==4) 
  {
    R->GT[0]=R->GT[1]=R->GT[2]=1;
    R->GT[3]=0;
  }
  else if (R->NL==2) 
  {
    R->GT[0]=R->GT[1]=1;
    R->GT[2]=0;
  }
  
  while( fgets(s,200,in) )
    if( strncmp(s,"|  number of weights    |",25) == 0 )
    {     
      sscanf(s+26,"%d",&R->NW);
      break;
    }

  while( fgets(s,200,in) )
    if( strncmp(s,"        minimum          ",25) == 0 )
      break;

  for(spc=0; spc < R->NSPC;spc++)
    for(n=0; n < R->NU[0];n++)
      fscanf(in,"%lf %lf %lf\n",&R->Rmin[spc][n],&R->Rmax[spc][n],&R->DR[spc][n]);

  fgets(s,200,in);
  fgets(s,200,in);
  fgets(s,200,in);

  V = make_d1D(R->NW);
  for(n=0; n < R->NW;n++)
    fscanf(in,"%lf\n",&V[n]);
  fclose(in);
  V2W(R,V);
  free_d1D(V);

  return;
}
//==============================================================================
int symb_atom(char *s)
{
  if ( strcmp(s,"H" ) == 0 || strcmp(s,"h" ) == 0 ) return 1   ;
  if ( strcmp(s,"He") == 0 || strcmp(s,"he") == 0 ) return 2   ;
  if ( strcmp(s,"Li") == 0 || strcmp(s,"li") == 0 ) return 3   ;
  if ( strcmp(s,"Be") == 0 || strcmp(s,"be") == 0 ) return 4   ;
  if ( strcmp(s,"B" ) == 0 || strcmp(s,"b" ) == 0 ) return 5   ;
  if ( strcmp(s,"C" ) == 0 || strcmp(s,"c" ) == 0 ) return 6   ;
  if ( strcmp(s,"N" ) == 0 || strcmp(s,"n" ) == 0 ) return 7   ;
  if ( strcmp(s,"O" ) == 0 || strcmp(s,"o" ) == 0 ) return 8   ;
  if ( strcmp(s,"F" ) == 0 || strcmp(s,"f" ) == 0 ) return 9   ;
  if ( strcmp(s,"Ne") == 0 || strcmp(s,"ne") == 0 ) return 10  ;
  if ( strcmp(s,"Na") == 0 || strcmp(s,"na") == 0 ) return 11  ;
  if ( strcmp(s,"Mg") == 0 || strcmp(s,"mg") == 0 ) return 12  ;
  if ( strcmp(s,"Al") == 0 || strcmp(s,"al") == 0 ) return 13  ;
  if ( strcmp(s,"Si") == 0 || strcmp(s,"si") == 0 ) return 14  ;
  if ( strcmp(s,"P" ) == 0 || strcmp(s,"p" ) == 0 ) return 15  ;
  if ( strcmp(s,"S" ) == 0 || strcmp(s,"s" ) == 0 ) return 16  ;
  if ( strcmp(s,"Cl") == 0 || strcmp(s,"cl") == 0 ) return 17  ;
  if ( strcmp(s,"Ar") == 0 || strcmp(s,"ar") == 0 ) return 18  ;
  if ( strcmp(s,"K" ) == 0 || strcmp(s,"k" ) == 0 ) return 19  ;
  if ( strcmp(s,"Ca") == 0 || strcmp(s,"ca") == 0 ) return 20  ;
  if ( strcmp(s,"Sc") == 0 || strcmp(s,"sc") == 0 ) return 21  ;
  if ( strcmp(s,"Ti") == 0 || strcmp(s,"ti") == 0 ) return 22  ;
  if ( strcmp(s,"V" ) == 0 || strcmp(s,"v" ) == 0 ) return 23  ;
  if ( strcmp(s,"Cr") == 0 || strcmp(s,"cr") == 0 ) return 24  ;
  if ( strcmp(s,"Mn") == 0 || strcmp(s,"mn") == 0 ) return 25  ;
  if ( strcmp(s,"Fe") == 0 || strcmp(s,"fe") == 0 ) return 26  ;
  if ( strcmp(s,"Co") == 0 || strcmp(s,"co") == 0 ) return 27  ;
  if ( strcmp(s,"Ni") == 0 || strcmp(s,"ni") == 0 ) return 28  ;
  if ( strcmp(s,"Cu") == 0 || strcmp(s,"cu") == 0 ) return 29  ;
  if ( strcmp(s,"Zn") == 0 || strcmp(s,"zn") == 0 ) return 30  ;
  if ( strcmp(s,"Ga") == 0 || strcmp(s,"ga") == 0 ) return 31  ;
  if ( strcmp(s,"Ge") == 0 || strcmp(s,"ge") == 0 ) return 32  ;
  if ( strcmp(s,"As") == 0 || strcmp(s,"as") == 0 ) return 33  ;
  if ( strcmp(s,"Se") == 0 || strcmp(s,"se") == 0 ) return 34  ;
  if ( strcmp(s,"Br") == 0 || strcmp(s,"br") == 0 ) return 35  ;
  if ( strcmp(s,"Kr") == 0 || strcmp(s,"kr") == 0 ) return 36  ;
  if ( strcmp(s,"Rb") == 0 || strcmp(s,"rb") == 0 ) return 37  ;
  if ( strcmp(s,"Sr") == 0 || strcmp(s,"sr") == 0 ) return 38  ;
  if ( strcmp(s,"Y" ) == 0 || strcmp(s,"y" ) == 0 ) return 39  ;
  if ( strcmp(s,"Zr") == 0 || strcmp(s,"zr") == 0 ) return 40  ;
  if ( strcmp(s,"Nb") == 0 || strcmp(s,"nb") == 0 ) return 41  ;
  if ( strcmp(s,"Mo") == 0 || strcmp(s,"mo") == 0 ) return 42  ;
  if ( strcmp(s,"Tc") == 0 || strcmp(s,"tc") == 0 ) return 43  ;
  if ( strcmp(s,"Ru") == 0 || strcmp(s,"ru") == 0 ) return 44  ;
  if ( strcmp(s,"Rh") == 0 || strcmp(s,"rh") == 0 ) return 45  ;
  if ( strcmp(s,"Pd") == 0 || strcmp(s,"pd") == 0 ) return 46  ;
  if ( strcmp(s,"Ag") == 0 || strcmp(s,"ag") == 0 ) return 47  ;
  if ( strcmp(s,"Cd") == 0 || strcmp(s,"cd") == 0 ) return 48  ;
  if ( strcmp(s,"In") == 0 || strcmp(s,"in") == 0 ) return 49  ;
  if ( strcmp(s,"Sn") == 0 || strcmp(s,"sn") == 0 ) return 50  ;
  if ( strcmp(s,"Sb") == 0 || strcmp(s,"sb") == 0 ) return 51  ;
  if ( strcmp(s,"Te") == 0 || strcmp(s,"te") == 0 ) return 52  ;
  if ( strcmp(s,"I" ) == 0 || strcmp(s,"i" ) == 0 ) return 53  ;
  if ( strcmp(s,"Xe") == 0 || strcmp(s,"xe") == 0 ) return 54  ;
  if ( strcmp(s,"Cs") == 0 || strcmp(s,"cs") == 0 ) return 55  ;
  if ( strcmp(s,"Ba") == 0 || strcmp(s,"ba") == 0 ) return 56  ;
  if ( strcmp(s,"La") == 0 || strcmp(s,"la") == 0 ) return 57  ;
  if ( strcmp(s,"Ce") == 0 || strcmp(s,"ce") == 0 ) return 58  ;
  if ( strcmp(s,"Pr") == 0 || strcmp(s,"pr") == 0 ) return 59  ;
  if ( strcmp(s,"Nd") == 0 || strcmp(s,"nd") == 0 ) return 60  ;
  if ( strcmp(s,"Pm") == 0 || strcmp(s,"pm") == 0 ) return 61  ;
  if ( strcmp(s,"Sm") == 0 || strcmp(s,"sm") == 0 ) return 62  ;
  if ( strcmp(s,"Eu") == 0 || strcmp(s,"eu") == 0 ) return 63  ;
  if ( strcmp(s,"Gd") == 0 || strcmp(s,"gd") == 0 ) return 64  ;
  if ( strcmp(s,"Tb") == 0 || strcmp(s,"tb") == 0 ) return 65  ;
  if ( strcmp(s,"Dy") == 0 || strcmp(s,"dy") == 0 ) return 66  ;
  if ( strcmp(s,"Ho") == 0 || strcmp(s,"ho") == 0 ) return 67  ;
  if ( strcmp(s,"Er") == 0 || strcmp(s,"er") == 0 ) return 68  ;
  if ( strcmp(s,"Tm") == 0 || strcmp(s,"tm") == 0 ) return 69  ;
  if ( strcmp(s,"Yb") == 0 || strcmp(s,"yb") == 0 ) return 70  ;
  if ( strcmp(s,"Lu") == 0 || strcmp(s,"lu") == 0 ) return 71  ;
  if ( strcmp(s,"Hf") == 0 || strcmp(s,"hf") == 0 ) return 72  ;
  if ( strcmp(s,"Ta") == 0 || strcmp(s,"ta") == 0 ) return 73  ;
  if ( strcmp(s,"W" ) == 0 || strcmp(s,"w" ) == 0 ) return 74  ;
  if ( strcmp(s,"Re") == 0 || strcmp(s,"re") == 0 ) return 75  ;
  if ( strcmp(s,"Or") == 0 || strcmp(s,"or") == 0 ) return 76  ;
  if ( strcmp(s,"Ir") == 0 || strcmp(s,"ir") == 0 ) return 77  ;
  if ( strcmp(s,"Pt") == 0 || strcmp(s,"pt") == 0 ) return 78  ;
  if ( strcmp(s,"Au") == 0 || strcmp(s,"au") == 0 ) return 79  ;
  if ( strcmp(s,"Hg") == 0 || strcmp(s,"hg") == 0 ) return 80  ;
  if ( strcmp(s,"Tl") == 0 || strcmp(s,"tl") == 0 ) return 81  ;
  if ( strcmp(s,"Pb") == 0 || strcmp(s,"pb") == 0 ) return 82  ;
  if ( strcmp(s,"Bi") == 0 || strcmp(s,"bi") == 0 ) return 83  ;
  if ( strcmp(s,"Po") == 0 || strcmp(s,"po") == 0 ) return 84  ;
  if ( strcmp(s,"At") == 0 || strcmp(s,"at") == 0 ) return 85  ;
  if ( strcmp(s,"Rn") == 0 || strcmp(s,"rn") == 0 ) return 86  ;
  if ( strcmp(s,"Fr") == 0 || strcmp(s,"fr") == 0 ) return 87  ;
  if ( strcmp(s,"Ra") == 0 || strcmp(s,"ra") == 0 ) return 88  ;
  if ( strcmp(s,"Ac") == 0 || strcmp(s,"ac") == 0 ) return 89  ;
  if ( strcmp(s,"Th") == 0 || strcmp(s,"th") == 0 ) return 90  ;
  if ( strcmp(s,"Pa") == 0 || strcmp(s,"pa") == 0 ) return 91  ;
  if ( strcmp(s,"Ul") == 0 || strcmp(s,"ul") == 0 ) return 92  ;
  if ( strcmp(s,"Np") == 0 || strcmp(s,"np") == 0 ) return 93  ;
  if ( strcmp(s,"Pu") == 0 || strcmp(s,"pu") == 0 ) return 94  ;
  if ( strcmp(s,"Am") == 0 || strcmp(s,"am") == 0 ) return 95  ;
  if ( strcmp(s,"Cm") == 0 || strcmp(s,"cm") == 0 ) return 96  ;
  if ( strcmp(s,"Bk") == 0 || strcmp(s,"bk") == 0 ) return 97  ;
  if ( strcmp(s,"Cf") == 0 || strcmp(s,"cf") == 0 ) return 98  ;
  if ( strcmp(s,"Es") == 0 || strcmp(s,"es") == 0 ) return 99  ;
  if ( strcmp(s,"Fm") == 0 || strcmp(s,"fm") == 0 ) return 100 ;

  return 0; //not recognized element  
} 
//==============================================================================
void  atom_symb(int i, char *s)      
{
  strcpy(s,"--");
  if ( i == 1  ) strcpy(s,"H" );
  if ( i == 2  ) strcpy(s,"He" );
  if ( i == 3  ) strcpy(s,"Li" );
  if ( i == 4  ) strcpy(s,"Be" );
  if ( i == 5  ) strcpy(s,"B" );
  if ( i == 6  ) strcpy(s,"C" );
  if ( i == 7  ) strcpy(s,"N" );
  if ( i == 8  ) strcpy(s,"O" );
  if ( i == 9  ) strcpy(s,"F" );
  if ( i == 10 ) strcpy(s,"Ne" );
  if ( i == 11 ) strcpy(s,"Na" );
  if ( i == 12 ) strcpy(s,"Mg" );
  if ( i == 13 ) strcpy(s,"Al" );
  if ( i == 14 ) strcpy(s,"Si" );
  if ( i == 15 ) strcpy(s,"P" );
  if ( i == 16 ) strcpy(s,"S" );
  if ( i == 17 ) strcpy(s,"Cl" );
  if ( i == 18 ) strcpy(s,"Ar" );
  if ( i == 19 ) strcpy(s,"K" );
  if ( i == 20 ) strcpy(s,"Ca" );
  if ( i == 21 ) strcpy(s,"Sc" );
  if ( i == 22 ) strcpy(s,"Ti" );
  if ( i == 23 ) strcpy(s,"V" );
  if ( i == 24 ) strcpy(s,"Cr" );
  if ( i == 25 ) strcpy(s,"Mn" );
  if ( i == 26 ) strcpy(s,"Fe" );
  if ( i == 27 ) strcpy(s,"Co" );
  if ( i == 28 ) strcpy(s,"Ni" );
  if ( i == 29 ) strcpy(s,"Cu" );
  if ( i == 30 ) strcpy(s,"Zn" );
  if ( i == 31 ) strcpy(s,"Ga" );
  if ( i == 32 ) strcpy(s,"Ge" );
  if ( i == 33 ) strcpy(s,"As" );
  if ( i == 34 ) strcpy(s,"Se" );
  if ( i == 35 ) strcpy(s,"Br" );
  if ( i == 36 ) strcpy(s,"Kr" );
  if ( i == 37 ) strcpy(s,"Rb" );
  if ( i == 38 ) strcpy(s,"Sr" );
  if ( i == 39 ) strcpy(s,"Y" );
  if ( i == 40 ) strcpy(s,"Zr" );
  if ( i == 41 ) strcpy(s,"Nb" );
  if ( i == 42 ) strcpy(s,"Mo" );
  if ( i == 43 ) strcpy(s,"Tc" );
  if ( i == 44 ) strcpy(s,"Ru" );
  if ( i == 45 ) strcpy(s,"Rh" );
  if ( i == 46 ) strcpy(s,"Pd" );
  if ( i == 47 ) strcpy(s,"Ag" );
  if ( i == 48 ) strcpy(s,"Cd" );
  if ( i == 49 ) strcpy(s,"In" );
  if ( i == 50 ) strcpy(s,"Sn" );
  if ( i == 51 ) strcpy(s,"Sb" );
  if ( i == 52 ) strcpy(s,"Te" );
  if ( i == 53 ) strcpy(s,"I" );
  if ( i == 54 ) strcpy(s,"Xe" );
  if ( i == 55 ) strcpy(s,"Cs" );
  if ( i == 56 ) strcpy(s,"Ba" );
  if ( i == 57 ) strcpy(s,"La" );
  if ( i == 58 ) strcpy(s,"Ce" );
  if ( i == 59 ) strcpy(s,"Pr" );
  if ( i == 60 ) strcpy(s,"Nd" );
  if ( i == 61 ) strcpy(s,"Pm" );
  if ( i == 62 ) strcpy(s,"Sm" );
  if ( i == 63 ) strcpy(s,"Eu" );
  if ( i == 64 ) strcpy(s,"Gd" );
  if ( i == 65 ) strcpy(s,"Tb" );
  if ( i == 66 ) strcpy(s,"Dy" );
  if ( i == 67 ) strcpy(s,"Ho" );
  if ( i == 68 ) strcpy(s,"Er" );
  if ( i == 69 ) strcpy(s,"Tm" );
  if ( i == 70 ) strcpy(s,"Yb" );
  if ( i == 71 ) strcpy(s,"Lu" );
  if ( i == 72 ) strcpy(s,"Hf" );
  if ( i == 73 ) strcpy(s,"Ta" );
  if ( i == 74 ) strcpy(s,"W" );
  if ( i == 75 ) strcpy(s,"Re" );
  if ( i == 76 ) strcpy(s,"Or" );
  if ( i == 77 ) strcpy(s,"Ir" );
  if ( i == 78 ) strcpy(s,"Pt" );
  if ( i == 79 ) strcpy(s,"Au" );
  if ( i == 80 ) strcpy(s,"Hg" );
  if ( i == 81 ) strcpy(s,"Tl" );
  if ( i == 82 ) strcpy(s,"Pb" );
  if ( i == 83 ) strcpy(s,"Bi" );
  if ( i == 84 ) strcpy(s,"Po" );
  if ( i == 85 ) strcpy(s,"At" );
  if ( i == 86 ) strcpy(s,"Rn" );
  if ( i == 87 ) strcpy(s,"Fr" );
  if ( i == 88 ) strcpy(s,"Ra" );
  if ( i == 89 ) strcpy(s,"Ac" );
  if ( i == 90 ) strcpy(s,"Th" );
  if ( i == 91 ) strcpy(s,"Pa" );
  if ( i == 92 ) strcpy(s,"Ul" );
  if ( i == 93 ) strcpy(s,"Np" );
  if ( i == 94 ) strcpy(s,"Pu" );
  if ( i == 95 ) strcpy(s,"Am" );
  if ( i == 96 ) strcpy(s,"Cm" );
  if ( i == 97 ) strcpy(s,"Bk" );
  if ( i == 98 ) strcpy(s,"Cf" );
  if ( i == 99 ) strcpy(s,"Es" );
  if ( i == 100) strcpy(s,"Fm" );
} 
//==============================================================================
double NP_VOL(Cell *C)
{
  double mass[100]={  1.0078,   3.0160,   6.0151,   9.0122,  10.0129,  12.0000,  14.0031,  15.9949,  18.9984,  19.9924,  22.9898,  23.9850,  26.9815,  27.9769,  30.9738,  31.9721,  34.9689,  35.9675,  38.9637,  39.9626,  44.9559,  45.9526,  49.9472,  49.9460,  54.9380,  53.9396,  58.9332,  57.9353,  62.9296,  63.9291,  68.9256,  69.9242,  74.9216,  73.9225,  78.9183,  77.9204,  84.9118,  83.9134,  88.9058,  89.9047,  92.9064,  91.9068,  96.9064,  95.9076, 102.9055, 101.9056, 106.9051, 105.9065, 112.9041, 111.9048, 120.9038, 119.9041, 126.9045, 123.9059, 132.9055, 129.9063, 137.9071, 135.9071, 140.9077, 141.9077, 144.9128, 143.9120, 150.9199, 151.9198, 158.9254, 155.9243, 164.9303, 161.9288, 168.9342, 167.9339, 174.9408, 173.9400, 179.9475, 179.9467, 184.9530, 183.9525, 190.9606, 189.9599, 196.9666, 195.9658, 202.9723, 203.9730, 208.9804, 208.9824, 209.9871, 210.9906, 223.0197, 223.0185, 227.0278, 230.0331, 231.0359, 233.0396, 236.0466, 238.0496, 241.0568, 243.0614, 247.0703, 249.0749, 252.0830, 257.0951}    ;

  int    i;
  double I[3][3],e[3][3],w[3];
  double M,a,b,c;
  double vol;

  for(i=0; i < C->N;i++)
    if( C->ATMZ[i] == 0 )
    {fprintf(stderr,"Error: fail to read species for moment of intertia matrix calculation %3d %3d!\n",i,C->ATMZ[i]);exit(1);}

  for(i=0,M=0.0; i < C->N;i++)
    M+=mass[C->ATMZ[i]];

  for(i=0,a=0.0,b=0.0,c=0.0; i < C->N;i++)
  {
    a+=mass[C->ATMZ[i]]*C->X[i][0];
    b+=mass[C->ATMZ[i]]*C->X[i][1];
    c+=mass[C->ATMZ[i]]*C->X[i][2];
  }

  //=====  diagonal elements  =====
  for(i=0,I[0][0]=(-(b*b/M)-(c*c/M)); i < C->N;i++)
    I[0][0]+=((mass[C->ATMZ[i]]*C->X[i][1]*C->X[i][1]+mass[C->ATMZ[i]]*C->X[i][2]*C->X[i][2]));

  for(i=0,I[1][1]=(-(a*a/M)-(c*c/M)); i < C->N;i++)
    I[1][1]+=((mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][0]+mass[C->ATMZ[i]]*C->X[i][2]*C->X[i][2]));

  for(i=0,I[2][2]=(-(b*b/M)-(a*a/M)); i < C->N;i++)
    I[2][2]+=((mass[C->ATMZ[i]]*C->X[i][1]*C->X[i][1]+mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][0]));

  //=====  off diagonal elements  =====
  for(i=0,I[0][1]=a*b/M; i < C->N;i++)
    I[0][1]-=(mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][1]);

  I[1][0]=I[0][1];

  for(i=0,I[0][2]=a*c/M; i < C->N;i++)
    I[0][2]-=(mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][2]);

  I[2][0]=I[0][2];

  for(i=0,I[1][2]=b*c/M; i < C->N;i++)
    I[1][2]-=(mass[C->ATMZ[i]]*C->X[i][1]*C->X[i][2]);
  I[2][1]=I[1][2];

  EV(I,e,w);

  vol=(4.0*Pi/3.0)*pow(2.5/M,1.5)*pow(w[0]*w[1]*w[2],0.5);
  if(0)
  {
    printf("Particle volume: % 10.4lf   A^3\n",vol);
    printf("Volume per atom: % 10.4lf   A^3/atom\n",vol/C->N);
    printf("Particle radius: % 10.4lf   A\n",pow(pow(2.5/M,1.5)*pow(w[0]*w[1]*w[2],0.5),0.333333));
  }

  return vol;
}
