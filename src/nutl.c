#include "nutl.h"

//==============================================================================
void READ_STRAT(ANN *R)
{
  int      i,j,k,m;
  char     s[700],s2[200];
  FILE    *mlp;
  int      I,W,T,J;
  int      indx,indm,sub_mdl,sub_ele,SKPR,SKPW;
  int      NG[20],NC[20];
  double  *M;
  double **X;
  int      spc[20],skp[20];

  // initialize sub-system models: depends on the number of elements
  if(R->NSPC == 2)
    for(i=0,k=0; i < R->NSPC;i++)
    {
      atom_symb(R->SPCZ[i],s);
      sprintf(R->file_list[k],"%s/%s.dat",R->otpt,s);
      k++;
    }
  if(R->NSPC == 3)
    for(i=0,k=0; i < (R->NSPC-1);i++)
      for(j=i+1; j < R->NSPC;j++)
      {
	atom_symb(R->SPCZ[i],s );
	atom_symb(R->SPCZ[j],s2);
	sprintf(R->file_list[k],"%s/%s%s.dat",R->otpt,s,s2);
	k++;
      }
    
  // read symm_func information from the first binary model
  for(i=0; i < 20;i++) 
    NG[i]=NC[i]=0;
  if( (mlp=fopen(R->file_list[0],"r")) == 0 )
  {
    printf("ERROR opening %s file to read basis parameters\n",s);         
    exit(1);
  }
  while( fgets(s,200,mlp) )
    if( strncmp(s,"  B2A",5) == 0 )
      break;
  for(i=0; i < 4;i++) 
    fgets(s,200,mlp);
  for(i=0; i < R->NSYM;i++)
  {
    fgets(s,200,mlp);
    sscanf(s,"%d %d %d %d %d %d %d %d %d\n",&k,&j,&k,&k,&k,&k,&k,&k,&k);
    NG[j]++;
  }
  fclose(mlp);

  // initialize variables: depends on the number of elements
  I = W = 0;
  if(R->NSPC == 2)
  {
    R->nw  = R->NSPC*(NG[2]+2*NG[4])*R->NU[1];  // total number of adjustable weights
    R->O   = NG[2] + NG[4];                     // starting point for adjustable weights
    I      = NG[2] + NG[4];                     // input_comp per element: elemental
    SKPR   = 0;                                 // where to skip reading
    SKPW   = 0;                                 // where to skip writing
    sub_mdl= 2;                                 // number of subsystems
    sub_ele= 1;                                 // number of elements per subsystem
    spc[0] = 0; spc[1] = 1;                     // list of elements in all models
    skp[0] = 0; skp[1] = 0;
  }
  if(R->NSPC == 3)
  {
    R->nw  = R->NSPC*(NG[4])*R->NU[1];  // total number of adjustable weights
    R->O   = 3*NG[2] + 5*NG[4];         // starting point for adjustable weights |AA|AAA|...
    I      = 2*NG[2]+3*NG[4];           // input_comp per element: binary  
    SKPR   = NG[2]+NG[4];               // where to skip reading
    SKPW   = 2*NG[2]+3*NG[4];           // where to skip writing
    sub_mdl= 3;                         // number of subsystems
    sub_ele= 2;                         // number of elements per subsystem
    spc[0] = 0; spc[1] = 1; spc[2] = 0; spc[3] = 2; spc[4] = 1; spc[5] = 2;
    skp[0] = 0; skp[1] = 1; skp[2] = 1; skp[3] = 0; skp[4] = 0; skp[5] = 1;
  }
  
  // allocate arrays and list of elements for all sub-models
  for(i=1,W=R->NU[1]*(I+1); i < R->NL-1;i++)
    W += (R->NU[i]+1)*R->NU[i+1];        // weight+bias per element
  M = make_d1D(W*sub_mdl*sub_ele);
  X = make_d2D(I*sub_mdl*sub_ele,3);
  
  // ========== read subsystem 'model' files
  indm = 0;
  indx = 0;
  for(i=0; i < sub_mdl;i++)
  {
    mlp = fopen(R->file_list[i],"r");
    while(fgets(s,200,mlp))
      if(strncmp(s,"        minimum",15)==0)
        break;
    for(j=0;j<(sub_ele*I);j++)
    {
      fgets(s,200,mlp);
      sscanf(s,"%lf %lf %lf\n",&X[indx][0],&X[indx][1],&X[indx][2]);
      indx++;
    }    
    while(fgets(s,200,mlp))
      if(strncmp(s,"|                        neural",31)==0)
        break;
    fgets(s,200,mlp);
   
    for(j=0;j<(sub_ele*W);j++)
    {
      fgets(s,200,mlp);
      sscanf(s,"%lf\n",&M[indm++]);
    }
    fclose(mlp);
  }

  // ========== map parameters from the subsystem into new 'model'
  indm = 0;
  indx = 0;
  for(i=0;i<(sub_mdl*sub_ele);i++)
  {
    T = spc[i];
    J = skp[i];
    // biases first
    for(j=0;j<R->NL-1;j++)
      for(k=0;k<R->NU[j+1];k++)
	R->B[T][j][k] = M[indm++];
    // input layer weights; SKIP?
    for(j=0;j<R->NU[1];j++)
    {
      indm += J*SKPR;
      for(k=0 ; k < (I-J*SKPR) ; k++)
	R->W[T][0][j][k+J*SKPW] = M[indm++];
    }
    // rest of the neurons
    for(j=1;j<R->NL-1;j++)
      for(k=0;k<R->NU[j+1];k++)
	for(m=0;m<R->NU[j];m++)
	  R->W[T][j][k][m] = M[indm++];
    // ranges
    indx += J*SKPR;
    for(j=0;j<(I-J*SKPR);j++)
    {
      R->DR[T][j+J*SKPW]   = X[indx][2];
      R->Rmax[T][j+J*SKPW] = X[indx][1];
      R->Rmin[T][j+J*SKPW] = X[indx][0];
      indx++;
    }
  }
  free_d1D(M);
  free_d2D(X,I*sub_mdl*sub_ele);
}
//==============================================================================
//                                            
//==============================================================================
void ANA_STR(ANN *R)
{
  char   s[400],s1[400];
  FILE   *in;
  int    N,MAX,COMP,ns,spc[10],t[10];
  int    i,n,k,tag=0,version=0;

  sprintf(s,"%s/stamp.dat",R->data);
  if( (in=fopen(s,"r")) == 0 )
  {
    fprintf(stderr,"%s not found. Please parse the data in %s first\n",s,R->data);
    exit(0);
  }
  while( fgets(s,200,in) ) 
  {
    if(strncmp(s,"PRSE",4) == 0) sscanf(s,"%s %d",s1,&N);    
    if(strncmp(s,"COMP",4) == 0) sscanf(s,"%s %d",s1,&COMP);    
    if(strncmp(s,"NSPC",4) == 0) sscanf(s,"%s %d",s1,&ns);    
    if(strncmp(s,"TSPC",4) == 0)
    { 
      for(i=0,n=0,k=4; i < ns;i++,k+=n) sscanf(s+k,"%d%n",&t[i],&n);
      for(i=0; i < ns;i++) spc[i]=t[i];
    }
    if(strncmp(s,"NMAX",4) == 0) sscanf(s,"%s %d",s1,&MAX);    
    if(strncmp(s,"VERS",4) == 0) version = 1;
  }
  fclose(in);

  if(version == 0)
  {printf("Error: the data is parsed with an old version of maise (before 2.4.xx)!\n");exit(1);}

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
int MARK_CL(Cell *C, int M, long seed)
{
  int i,j,m,k;

  if( seed == 0 ) seed=time(NULL);
  PlantSeeds(seed);
  if( M == C->N )
  {
    for(i=0; i < C->N;i++)
      C->FRC[i] = 1;

    return C->N*3;
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
  else
    return k*3;
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
void SAVE_ANN(ANN *R, double E_TIME,double C_TIME)
{
  int    i,spc;
  char   s[600],b[200],t[300],e[200];
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
  sprintf(t," %08lX",(long)(sum*1e11)%(long)(1e8));
  fprintf(out,"|  model unique ID      | %s%s",t,e+22+strlen(t));

  sprintf(t," %d",R->NSPC);
  fprintf(out,"|  number of species    | %s%s",t,e+22+strlen(t));
  for(i=0,t[0]=0; i < R->NSPC;i++)
    sprintf(t+strlen(t)," %2d",R->SPCZ[i]);
  fprintf(out,"|  species types        | %s%s",t,e+22+strlen(t));
  for(i=0,t[0]=0,sprintf(t+strlen(t)," "); i < R->NSPC;i++,sprintf(t+strlen(t)," "))
    atom_symb(R->SPCZ[i],t+strlen(t));
  fprintf(out,"|  species names        | %s%s",t,e+22+strlen(t));

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
  fprintf(out,"|  test  energy data    | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->NF);
  fprintf(out,"|  train force  data    | %s%s",t,e+22+strlen(t));
  sprintf(t," %d",R->TNF);
  fprintf(out,"|  test  force  data    | %s%s",t,e+22+strlen(t));
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
  if ( fabs(C_TIME/(double)R->NP-E_TIME)/E_TIME > 0.5 )
    sprintf(t," %1.2lf  %1.2lf sec-",C_TIME,E_TIME);
  else
    sprintf(t," %1.2lf  %1.2lf sec",C_TIME,E_TIME);
  fprintf(out,"|  cpu and wall times   | %s%s",t,e+22+strlen(t));
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
  char   s[400],s1[200];
  int    i,n,spc,k,l,ver;
  FILE   *in;
  double *V;

  sprintf(s,"model");
  if( R->JOBT/10 == 4 )
    sprintf(s,"%s/model",R->otpt);

  ver=check_ver(s);
  if( ver < 2400)
  {
    printf("Error: the model file is in old format; convert the model first!\n");
    exit(1);
  }
    
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
      printf("Error: the model is in old format (contains reference energies)!\n");
      exit(1);
    }
    else if( strncmp(s,"|  number of layers     |",25) == 0 )
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

  if ((R->JOBT/10)!=4)
  {
    for(spc=0;spc<R->NSPC;spc++)
      for(n=0;n<R->NU[1];n++)
        for(l=0;l<R->NU[0];l++)
          R->W[spc][0][n][l] *= R->DR[spc][l];
    for(spc=0;spc<R->NSPC;spc++)
      for(l=0;l<R->NU[0];l++)
	R->Rmin[spc][l] /= R->DR[spc][l];
  }

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
  int    i;
  double I[3][3],e[3][3],w[3];
  double M,a,b,c;
  double vol;

  for(i=0; i < C->N;i++)
    if( C->ATMZ[i] == 0 )
    {fprintf(stderr,"Error: fail to read species for moment of intertia matrix calculation %3d %3d!\n",i,C->ATMZ[i]);exit(1);}

  for(i=0,M=0.0; i < C->N;i++)
    M+=C->mass[C->ATMZ[i]];

  for(i=0,a=0.0,b=0.0,c=0.0; i < C->N;i++)
  {
    a+=C->mass[C->ATMZ[i]]*C->X[i][0];
    b+=C->mass[C->ATMZ[i]]*C->X[i][1];
    c+=C->mass[C->ATMZ[i]]*C->X[i][2];
  }

  //=====  diagonal elements  =====
  for(i=0,I[0][0]=(-(b*b/M)-(c*c/M)); i < C->N;i++)
    I[0][0]+=((C->mass[C->ATMZ[i]]*C->X[i][1]*C->X[i][1]+C->mass[C->ATMZ[i]]*C->X[i][2]*C->X[i][2]));

  for(i=0,I[1][1]=(-(a*a/M)-(c*c/M)); i < C->N;i++)
    I[1][1]+=((C->mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][0]+C->mass[C->ATMZ[i]]*C->X[i][2]*C->X[i][2]));

  for(i=0,I[2][2]=(-(b*b/M)-(a*a/M)); i < C->N;i++)
    I[2][2]+=((C->mass[C->ATMZ[i]]*C->X[i][1]*C->X[i][1]+C->mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][0]));

  //=====  off diagonal elements  =====
  for(i=0,I[0][1]=a*b/M; i < C->N;i++)
    I[0][1]-=(C->mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][1]);

  I[1][0]=I[0][1];

  for(i=0,I[0][2]=a*c/M; i < C->N;i++)
    I[0][2]-=(C->mass[C->ATMZ[i]]*C->X[i][0]*C->X[i][2]);

  I[2][0]=I[0][2];

  for(i=0,I[1][2]=b*c/M; i < C->N;i++)
    I[1][2]-=(C->mass[C->ATMZ[i]]*C->X[i][1]*C->X[i][2]);
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
