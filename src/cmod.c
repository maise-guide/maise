#include "cmod.h"

extern const double Pi;
extern const double eV2GPa;
const double kb = 11604.0;
extern const double atom_mass[96];
extern const double mindist[96];
extern const double minpres[96];
extern const double maxpres[96];

//========================================================================
double CELL_KIN(Cell *C)
{
  int i,q;

  for(i=0,C->K=0.0;i<C->N;i++)
    for(q=0;q<D3;q++)
      C->K += C->m[C->ATMZ[i]]*C->V[i][q]*C->V[i][q];

  C->K *= 0.5;

  return C->K;
}
//=======================================================================
double Lindemann(Cell *C, int J)
{
  int i,j,q;
  double x,r,LI;

  if (J == 0)
    {
      for(i=0,C->LI=0;i<C->N;i++)
	for(j=i+1;j<C->N;j++)
	  C->R1[i][j] = C->R2[i][j] = 0.0;
      return 0.0;
    }
  
  if( J == 1 )
    {
      C->LI++;
      LI = 0.0;
      for(i=0;i<C->N;i++)
	for(j=i+1;j<C->N;j++)
	  {
	    for(q=0,r=0.0;q<D3;q++)
	      {
		x  = C->X[i][q] - C->X[j][q];
		r += x*x;
	      }
	    C->R2[i][j] += r;
	    C->R1[i][j] += sqrt(r);
	    if( (r=C->R2[i][j]*(double)C->LI - C->R1[i][j]*C->R1[i][j]) > 1e-12 )
	      LI += sqrt( r ) / C->R1[i][j];
	  }
      return 2.0*LI/(double)(C->N*(C->N-1));
    }  

  return -1.0;
}
//=========================================================================
double CELL_ENE(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{   
  int i;

  if( R->MODT==1 )
  {
    P->EFS = 0;
    PARS_STR(P,W,C,L,0,"."); 
    C->E = ENE_ANN(R,L);

    for(i=0;i<C->N;i++)
      C->EA[i] = L->EA[i];
  }
  if( R->MODT>1 )
    C->E = ENE_POT(C);

  C->H = C->E;
  if(C->ND==3)
    C->H += C->p*CELL_VOL(C);
  return C->H;
}
//=========================================================================
double CELL_FRC(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int i,q;

  if( R->MODT==1 )
  {
    P->EFS = 3;
    PARS_STR(P,W,C,L,0,".");
    C->E = FRC_ANN(R,L);

    for(i=0;i<C->N;i++)
      C->EA[i] = L->EA[i];
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
	C->F[i][q] = L->f[i][q];
    for(q=0;q<6;q++)
      C->U[q] = L->s[q];
  }
  if( R->MODT>1 )
    C->E = FRC_POT(C);

  //=====  remove noise to preserve symmetry in some cases =====
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      if( fabs(C->F[i][q])<1e-8 )
	C->F[i][q] = 0.0;
  for(q=0;q<6;q++)
    if( fabs(C->U[q])<1e-8 )
      C->U[q] = 0.0;

  C->H = C->E;

  if(C->ND==3)
    C->H += C->p*CELL_VOL(C);
  return C->H;
}
//=========================================================================
void  MDOUT(Cell *C,double PE,double KE,double LV,double IP,double DE,double LI, int PIPE, int DISK, int NS, int sam, double T0,double E0,double li)
    // PUTS the center of mass at the origin!
{
  FILE   *out;
  char   file[400],s[200];
  double steps;
  
  steps = (double) (sam);

  if(PIPE && steps > 0.0)
    printf("ave. %8d:  U= % 12.6lf   K= % 12.6lf   IP= % 12.6lf   LV= %12.6lf  DE=% 12.6lf  LI=% 12.6lf\n",NS*(sam-1),PE/steps,KE/steps,IP/steps,LV/steps,DE/steps,li);

  if(DISK==0)
    return;
  sprintf(file,"%s/TEMP-%04d.dat",C->WDIR,(int)T0);
  if(sam == 1)
    out = fopen(file,"w");
  else
    out = fopen(file,"a");
  fprintf(out,"%8d  % 12.6lf  % 12.6lf  % 12.6lf  % 12.6lf  % 12.6lf  % 12.6lf\n",NS*(sam-1),C->H/(double)C->N,C->K*kb/1.5/(double)C->N,C->pint*eV2GPa,C->LAT[0],(C->H+C->K)/(double)C->N-E0,li);
  fclose(out);
  JAR(C);

  // save the current snapsshot
  strcpy(s,C->TAG);
  sprintf(C->TAG," %s (MD) : temp % 10.3lf  steps % 8d",VERSION,C->K*kb/1.5/(double)C->N,NS*(sam-1));
  sprintf(file,"%s/CONTCAR",C->WDIR);
  SAVE_CELL(C,file,0);
  strcpy(C->TAG,s);
}
//=========================================================================
void Maxwell(Cell *C, double T, int therm)
{
  int    i,q,t,fixed;
  double x,v,f,Q,V[3];
  int move[C->N];

  if( (therm % 10) > 0)
  {
    if(C->POS == 1)
      return;
    else
    {printf("Error: the MD job needs velocities provided in POSCAR file!\n");exit(1);}
  }

  for(i=0,fixed=0;i<C->N;i++)
  {
    for(q=0,t=0;q<3;q++)
      if(C->FF[i][q]==0)
	t++;
    if(t==3)
    {move[i]=0;fixed++;}
    else
      move[i]=1;
  }

  v = sqrt(3.0*T);
  for(q=0;q<3;q++)
    V[q] = 0.0;
  for(i=0;i<C->N;i++)
    if(move[i]==1)
    {
      f = Random()*2.0*Pi;
      Q = Random()    *Pi;
      C->V[i][0] = v * sin(Q)*cos(f)/sqrt(C->m[C->ATMZ[i]]);
      C->V[i][1] = v * sin(Q)*sin(f)/sqrt(C->m[C->ATMZ[i]]);
      C->V[i][2] = v * cos(Q)       /sqrt(C->m[C->ATMZ[i]]);
      for(q=0;q<3;q++)
	V[q] += C->V[i][q];
    }
    else
      C->V[i][0] = C->V[i][1] = C->V[i][2] = 0.0;

  for(i=0;i<C->N;i++)
    if(move[i]==1)
      for(q=0;q<3;q++)
	C->V[i][q] -= V[q]/(double)(C->N-fixed);

  for(i=0,x=0.0;i<C->N;i++)
    for(q=0;q<3;q++)
      x += C->V[i][q]*C->V[i][q]*(double)C->m[C->ATMZ[i]];
  x = v/sqrt(x/(double)C->N);
  
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      C->V[i][q] *= x;

}
//========================================================================= 
//========================================================================= 
void Dynamics(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L, double T0, int PIPE,char *fname) 
{ 
  int    i,k,q,n,NS,NI,spc1,spc2,sam,x,NRDF; 
  double dt,fi,xi,fs,tau,li,T,vi,ktt,E0,a; 
  FILE   *out; 
  char   file[400],s[200]; 
  int    move[C->N],fixed;
  double KE,PE,LI,LV,IP,DE;
  double ***H,**B,**D;
  double wall_time,user_time,syst_time,cpus_time;
  struct timeval t1, t2;
  struct rusage t3;
  int    RUNT;

  // to make sure velocities will be outputted
  C->POS = 1;

  // type of the MD run
  RUNT = R->MDTP / 10;

  // initiate some arrays
  B = make_d2D(C->N,3); // average position of atoms
  D = make_d2D(C->N,3); // net displacement of atoms
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {  
      D[i][q] = 0.0;
      B[i][q] = C->X[i][q];
    }

  // find which atoms are fixed!
  for(i=0,fixed=0;i<C->N;i++)
  {
    for(q=0,k=0;q<3;q++)
      if(C->FF[i][q]==0)
	k++;
    if(k==3)
    {move[i]=0;fixed++;}
    else
      move[i]=1;
  }

  // initiate the timing
  gettimeofday(&t1, NULL);
  getrusage(RUSAGE_SELF, &t3);
  syst_time  = (double) t3.ru_stime.tv_sec + (double) t3.ru_stime.tv_usec / (double) 1e6; // = kernel time
  user_time  = (double) t3.ru_utime.tv_sec + (double) t3.ru_utime.tv_usec / (double) 1e6; // = user time
  cpus_time  = (syst_time + user_time);

  // initiate the RDF
  H = make_d3D(C->NRDF,C->MNT,C->MNT);
  NRDF = 0;
  for(k=0;k<C->NRDF;k++)
    for(spc1=0;spc1<C->NSPC;spc1++)
      for(spc2=0;spc2<C->NSPC;spc2++)
	H[k][spc1][spc2] = 0.0;

  // initiate the required variables
  IP   = 0.0;
  LI   = 0.0;
  KE   = 0.0;
  PE   = 0.0;
  DE   = 0.0;
  a    = 0.0;
  dt   = R->DELT;
  T    = T0/kb;
  tau  = R->CPLT;
  ktt  = R->ICMP/R->CPLP;
  sam  = 0;
  x    = 0;
  NS   = 10;
  if(R->NSTP<10)
    NS = 1;
  NI   = R->NSTP/NS;
  LV   = 0.0;
  for(i=0;i<3;i++) 
    C->LAT[i] = C->R0*VectorLen(C->L[i],3);    
  E0   = CELL_KIN(C)/(double)C->N+CELL_ENE(R,P,W,C,L)/(double)C->N;
  // initialize the Lindemann index and scaling factors
  Lindemann(C,0); 
  vi   = 1.0;
  xi   = 0.0;
  fs   = 0.0;

  for( n=0;n<NI;n++)
  {
    C->H    = CELL_FRC(R,P,W,C,L); // for output and book-keeping of energies
    C->K    = CELL_KIN(C);
    C->pint = CELL_PRS(C);
    li      = Lindemann(C,1); 
    LI     += li;
    KE     += C->K*kb/1.5/(double)C->N;
    PE     += C->H/(double)C->N;
    LV     += C->LAT[0];
    IP     += C->pint*eV2GPa;
    DE     += (C->K+C->H)/(double)C->N-E0;
    sam++;
    MDOUT(C,PE,KE,LV,IP,DE,LI,PIPE,1,NS,sam,T0,E0,li);

    for(k=0;k<NS;k++)
    {         
      x++; // MD steps index
      C->H    = CELL_FRC(R,P,W,C,L);
      C->K    = CELL_KIN(C);
      C->pint = CELL_PRS(C);

      for(i=0;i<C->N;i++) // update the RDFs
	if(NDX(C,i,0)>C->Rmax)
	  break;
      if(i==C->N)
      {
	NRDF++;
	RDF(C,1);
	for(q=0;q<C->NRDF;q++)
	  for(spc1=0;spc1<C->NSPC;spc1++)
	    for(spc2=0;spc2<C->NSPC;spc2++)
	      H[q][spc1][spc2] += C->RDF[q][spc1][spc2];
      }

      // ####################### update thermostat and barostat scalings
      if (RUNT == 2 || RUNT == 3)
      {
	fs = 1.0/(tau*tau)*(2.0*C->K/( (3.0*(double)C->N)*T)-1.0 );
	xi += dt*fs;
      }
      if (RUNT == 3 || RUNT == 4)
      {
	vi = pow(1.0 - (ktt*dt/3.0)*(C->p-C->pint)*eV2GPa,1.0/3.0);
      }

      // ####################### NVE and isobaric integrations
      if (RUNT == 1 || RUNT == 4)
      {
	for(i=0;i<C->N;i++)
	  for(q=0;q<3;q++)
	  {
	    if(move[i]==1)
	    {
	      C->X[i][q] += C->V[i][q]*dt + 0.5*C->F[i][q]*dt*dt/C->m[C->ATMZ[i]];
	      B[i][q]    += (C->V[i][q]*dt + 0.5*C->F[i][q]*dt*dt/C->m[C->ATMZ[i]])/2.0;
	      C->V[i][q] += 0.5*C->F[i][q]*dt/C->m[C->ATMZ[i]];      
	    }
	  }
	C->H = CELL_FRC(R,P,W,C,L);
	C->K = CELL_KIN(C);
	for(i=0;i<C->N;i++)
	  for(q=0;q<D3;q++)
	    if(move[i]==1)
	    {
	      C->V[i][q] = (C->V[i][q]+0.5*C->F[i][q]*dt/C->m[C->ATMZ[i]])/(1.0+a);
	      C->F[i][q] += -a*C->V[i][q];
	    }
	
	JAR(C);
      }
      
      // ######################### NVT, NPT integrations
      if (RUNT == 2 || RUNT == 3)
      {
	for(i=0;i<C->N;i++)
	  for(q=0;q<3;q++)   
	    if(move[i]==1)
	    {
	      fi          = C->F[i][q]/C->m[C->ATMZ[i]]-xi*C->V[i][q];    
	      C->X[i][q] += C->V[i][q]*dt + 0.5*fi*dt*dt;
	      B[i][q]    += (C->V[i][q]*dt + 0.5*fi*dt*dt)/2.0;
	      C->V[i][q] += fi*dt;
	    }
	
	JAR(C);
      }

      // ######################### rescaling for constant NPT and isobaric runs
      if (RUNT == 3 || RUNT == 4)
      {
	for(i=0;i<3;i++)
	{
	  for(q=0;q<3;q++)
	  {
	    C->L[i][q] *= vi;
	    C->X[i][q] *= vi;
	  }
	  C->LAT[i] = C->R0*VectorLen(C->L[i],3);
	}
	Reciprocal(C);
      }
      
      // ######################### MD movie
      if( R->MOVI> 0 && x % R->MOVI == 0 )
      {
	sprintf(file,"mkdir -p %s/movie%04d",C->WDIR,(int)T0);
	system(file);
	strcpy(s,C->TAG);
	sprintf(C->TAG," %s (MD) : temp % 10.3lf  steps % 8d",VERSION,C->K*kb/1.5/(double)C->N,NS*(sam-1));
	sprintf(file,"%s/movie%04d/POSCAR-%d",C->WDIR,(int)T0,x);
	SAVE_CELL(C,file,0);
	strcpy(C->TAG,s);
      }

    } 
  }      

  C->H    = CELL_FRC(R,P,W,C,L);
  C->K    = CELL_KIN(C);
  C->pint = CELL_PRS(C);
  li      = Lindemann(C,1);
  LI     += li;
  KE     += C->K*kb/1.5/(double)C->N;
  PE     += C->H/(double)C->N;
  LV     += C->LAT[0];
  DE     += (C->K+C->H)/(double)C->N-E0;
  IP     += C->pint*eV2GPa;
  sam++;
  MDOUT(C,PE,KE,LV,IP,DE,LI,PIPE,1,NS,sam,T0,E0,li);

  // compute wall and cpu time
  getrusage(RUSAGE_SELF, &t3);
  gettimeofday(&t2, NULL);
  syst_time  = (double) t3.ru_stime.tv_sec + (double) t3.ru_stime.tv_usec / (double) 1e6; // = kernel time
  user_time  = (double) t3.ru_utime.tv_sec + (double) t3.ru_utime.tv_usec / (double) 1e6; // = user time
  wall_time  = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / (double) 1e6;
  cpus_time  = (syst_time + user_time) - cpus_time;

  //output the averages of E, T, and Lindemann index for each temperature
  out=fopen(fname,"a");
  fprintf(out,"% 12.4lf   % 12.6lf   % 12.6lf   % 12.6lf   % 12.6lf   % 12.6lf   % 12.6lf   % 10d   % 14.2lf",T0,PE/(double)sam,KE/(double)sam,IP/(double)sam,LV/(double)sam,DE/(double)sam,li,(sam-1)*NS,cpus_time);
  if(cpus_time/(double)R->NP < 0.5*wall_time)
    fprintf(out,"-\n");
  else
    fprintf(out,"\n");
  fclose(out);  
  printf("\n === % lf   % lf   % lf   % lf   % lf   % lf   % lf   % 10d   % 14.2lf\n\n",T0,PE/(double)sam,KE/(double)sam,IP/(double)sam,LV/(double)sam,DE/(double)sam,li,(sam-1)*NS,cpus_time);
  
  // output the final structure
  strcpy(s,C->TAG);
  sprintf(C->TAG," %s (MD) : temp % 10.3lf  steps % 8d",VERSION,C->K*kb/1.5/(double)C->N,NS*(sam-1));
  sprintf(file,"%s/CONTCAR-%04d",C->WDIR,(int)T0);
  SAVE_CELL(C,file,0);
  strcpy(C->TAG,s);

  // output the average position  of atoms
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      D[i][q]    = C->X[i][q];
      C->X[i][q] = B[i][q];
    }
  JAR(C);
  sprintf(file,"%s/AVERAGE-%04d",C->WDIR,(int)T0);
  SAVE_CELL(C,file,0);
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      C->X[i][q] = D[i][q];
  JAR(C);

  // output the average RDF
  for(q=0;q<C->NRDF;q++)
    for(spc1=0;spc1<C->NSPC;spc1++)
      for(spc2=0;spc2<C->NSPC;spc2++)
	C->RDF[q][spc1][spc2] = H[q][spc1][spc2]/(double)NRDF;
  sprintf(file,"%s/RDF-%04d.dat",C->WDIR,(int)T0);
  Print_RDF_FULL(C,file);

  free_d2D(B,C->N);
  free_d2D(D,C->N);
  free_d3D(H,C->NRDF,C->MNT);
}
//=========================================================================  
//=========================================================================
double CELL_PRS(Cell *C)
{
  int i;
  double v1 = CELL_VOL(C);

  for (i =0,C->pint=0.0;i<3;i++) 
    C->pint += C->U[i]/v1/3.0;

  return C->pint;
}
//=========================================================================  
//=========================================================================
void CELL_MD(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int n,N;
  double Tmin, Tmax, dT, T;
  char s[400],fname[400];
  FILE *out;

  if (R->MDTP/10 == 1)
    sprintf(s,"|                        NVE Molecular dynamics                       |");
  if (R->MDTP/10 == 2)
    sprintf(s,"|                        NVT Molecular dynamics                       |");
  if (R->MDTP/10 == 3)
    sprintf(s,"|                        NPT Molecular dynamics                       |");
  if (R->MDTP/10 == 4)
    sprintf(s,"|                      Isobaric Molecular dynamics                    |");

  printf("%s\n",s);
  printf("=======================================================================\n\n");

  if (R->MDTP < 10)
  {
    printf("Error: the MD job type (MDTP=%d) is not supported.\n",R->MDTP);
    exit(1);
  }

  sprintf(fname,"%s/ave-out.dat",C->WDIR);
  if(access(fname,F_OK)==-1)
  {
    out=fopen(fname,"w");
    FPRNT_HEAD(out);
    fprintf(out,"%s\n",s);
    fprintf(out,"=======================================================================\n\n");
    fprintf(out,"   Target_T       Potential_E    Kinetic_E       Intern_P      Lattice_V       E_change      Lindemann        Steps          CPU_time\n");
    fprintf(out,"     (K)           (eV/atom)        (K)           (GPa)           (A)          (eV/atom)                                        (s) \n\n");
    fclose(out);
  }

  Tmin = R->TMIN;
  Tmax = R->TMAX;
  dT   = R->TSTP;
  
  sprintf(s,"mkdir -p %s",C->WDIR);
  system(s);
  N = (int)fabs((Tmax-Tmin)/dT);
  Maxwell(C,Tmin/kb,R->MDTP);
  P->EFS = 1;
  for(n=0;n<=N;n++)
  {
    T = Tmin + (double)n*dT;
    Dynamics(R,P,W,C,L,T,1,fname);
  }
  sprintf(C->TAG," %s (MD) : temp % 10.3lf  steps % 8d",VERSION,C->K*kb/1.5/(double)C->N,R->NSTP);
}
//======================================================
void CELL_OUT(Cell *C)
{
  int    i;
  double V;
  char   s[200],b[200],t[200];
  FILE   *out;
  time_t rawtime;
  struct tm * timeinfo;

  sprintf(s,"OUTCAR");
  out = fopen(s,"a");
  
  abc(C);
  sprintf(b,"----------------------------------------------------------------------------------------------------------\n");

  if( C->it==0 )
  {
    time(&rawtime);
    timeinfo = localtime ( &rawtime );
    sprintf(t," %s",asctime(timeinfo));
    t[strlen(t)-1] = 0;
    fprintf(out,"%s    Cell optimization    with %s   model %s   on %s                     \n%s",b,C->VER,C->ID,t,b);
    sprintf(t," %d",C->MINT);
    fprintf(out,"   optimizer type         %s\n",t);
    sprintf(t," %d",C->RLXT);
    fprintf(out,"   relaxation type        %s\n",t);
    sprintf(t," %1.2lf GPa",C->p*eV2GPa);
    fprintf(out,"   target pressure        %s\n",t);

    sprintf(t," %d",C->ND);
    fprintf(out,"   dimensionality         %s\n",t);
    for(i=0,t[0]=0,sprintf(t+strlen(t)," ");i < C->NSPC;i++,sprintf(t+strlen(t),"    "))
      atom_symb(C->SPCZ[i],t+strlen(t));
    fprintf(out,"   species                %s\n",t);
    sprintf(  t,"   atoms of each species                     ");
    for(i=0;i < C->NSPC;i++)
      sprintf(t+27+6*i,"%d           ",C->SPCN[i]);
    fprintf(out,"%s\n",t);
    sprintf(t," %d",C->N);
    fprintf(out,"   total number of atoms  %s\n",t);
    fprintf(out,"%s\n\n\n",b);
  }  
  V = CELL_VOL(C);
  fprintf(out,"%s                   LATTICE VECTORS                                                 VOLUME (Angst^3/atom) \n%s",b,b);
  fprintf(out,"      % 12.6lf  % 12.6lf  % 12.6lf  % 45.6lf vol\n",C->L[0][0],C->L[0][1],C->L[0][2],V/(double)C->N);
  fprintf(out,"      % 12.6lf  % 12.6lf  % 12.6lf\n",C->L[1][0],C->L[1][1],C->L[1][2]);
  fprintf(out,"      % 12.6lf  % 12.6lf  % 12.6lf\n",C->L[2][0],C->L[2][1],C->L[2][2]);
  if(C->OUT/10>0)
  {
    fprintf(out,"%s                 POSITION (Angst)                     TOTAL-FORCE (eV/Angst)           ATOM ENERGY (eV)\n%s",b,b);
    for(i=0;i<C->N;i++)
      fprintf(out,"   % 12.6lf % 12.6lf % 12.6lf   % 12.6lf % 12.6lf % 12.6lf % 17.12lf\n",C->X[i][0],C->X[i][1],C->X[i][2],C->F[i][0],C->F[i][1],C->F[i][2],C->EA[i]);
    fprintf(out,"%s",b);
    fprintf(out,"  Total   ");
    for(i=0;i<6;i++)
      fprintf(out,"% 15.8lf ",C->U[i]/CELL_VOL(C));
    fprintf(out,"\n");
    fprintf(out,"  in kB   ");
    for(i=0;i<6;i++)
      fprintf(out,"% 15.8lf ",C->U[i]*eV2GPa*10.0/CELL_VOL(C));
    fprintf(out,"\n%s",b);
  }
  fprintf(out,"  iter %3d   total enthalpy= % 14.8lf   energy=  % 14.8lf   % 14.8lf  % 14.8lf\n%s\n\n\n",C->it,C->H,C->E,C->H/(double)C->N,C->E/(double)C->N,b);

  fclose(out);
  sprintf(s,"OSZICAR");
  out = fopen(s,"a");
  fprintf(out,"                            % 1.8lf\n",C->H);
  fclose(out);
}
//======================================================
int CELL_OK(Cell *C)
{
  int i;
  
  for(i=0;i<C->N;i++)
    if(NDX(C,i,0)<C->Rm[C->ATMZ[i]]+C->Rm[C->ATMZ[C->Ni[i][0]]])
      return 0;
  
  return 1;
}
//======================================================
double TOPK(double *x, int N, int K)
{
  int    n,k;
  double X[K];

  for(k=0;k<K;k++)
    X[k] = 0.0;
  for(n=0;n<N;n++)
    if( x[n]>X[0] )
    {
      X[0] = x[n];
      for(k=1;k<K;k++)
	if( X[k-1]>X[k] )
	  dSwap(&X[k-1],&X[k]);
    }
  return X[0];
}
//======================================================
void CELL_PHON(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int    i,q,j,k,N,*I;
  double dx,*A,*B,*e,*b,*r,x[3],v[3],bm,rm;
  FILE   *out;
  double wall_time,user_time,syst_time,cpus_time;
  double thz=1e12;           // THz to Hz
  double aum=1.66054e-27;    // Atomic Unit Mass to Kg
  double ang=1e-10;          // Angstrom to Meter
  double evj=1.60218e-19;    // eV to J
  double thzcm=33.35641;     // THz to 1/Cm
  struct timeval t1, t2;
  struct rusage t3;

  printf("|                          Phonon calculation                         |\n");
  printf("=======================================================================\n\n");

  N = C->N*D3;

  A  = make_d1D(N*N);
  B  = make_d1D(N*N);
  e  = make_d1D(N*N);
  b  = make_d1D(N);
  r  = make_d1D(N);
  I  = make_i1D(N);

  // initiate the timing
  gettimeofday(&t1, NULL);

  if(C->EVOK==0)
  {
    C->ev   = make_d1D(N);
    C->EV   = make_d2D(N,N);
    C->EVOK = 1;
  }

  CELL_FRC(R,P,W,C,L);
  C->OUT = 10;
  system("rm -f OUTCAR");
  CELL_OUT(C);

  dx = C->DISP;

  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      C->X[i][q] += dx;
      CELL_FRC(R,P,W,C,L);
      for(j=0;j<C->N;j++)
        for(k=0;k<3;k++)
          A[(i*3+q)*N+j*3+k] = C->F[j][k];
      C->X[i][q] -= 2.0*dx;
      CELL_FRC(R,P,W,C,L);
      for(j=0;j<C->N;j++)
        for(k=0;k<3;k++)
          A[(i*3+q)*N+j*3+k] = 0.5*(C->F[j][k]-A[(i*3+q)*N+j*3+k])/dx / sqrt(C->mass[C->ATMZ[i]]*C->mass[C->ATMZ[j]]);
      C->X[i][q] += dx;
    }
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      for(j=0;j<C->N;j++)
        for(k=0;k<3;k++)
	  B[(i*3+q)*N+j*3+k] = 0.0;

  if(0)
    for(i=0;i<C->N;i++)
      for(q=0;q<3;q++)
      {
	for(j=0;j<C->N;j++)
	  for(k=0;k<3;k++)
	    printf("% 16.12lf ",A[(i*3+q)*N+j*3+k]);
	printf("\n");
      }

  EVN(A,B,e,b,N);

  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      C->ev[i*3+q] = b[i*3+q] * evj / (thz*thz*ang*ang*aum);
      C->ev[i*3+q] = b[i*3+q] < 0.0 ? -sqrt(fabs(C->ev[i*3+q])) : sqrt(fabs(C->ev[i*3+q]));
      for(j=0;j<C->N;j++)
        for(k=0;k<3;k++)
          C->EV[i*3+q][j*3+k] = e[(i*3+q)*N+j*3+k];
    }

  if(0)
    for(i=0;i<C->N*3;i++)
      printf("%3d % 24.16lf % 24.16lf\n",i,C->ev[i],C->EV[N-7+C->ND][i]);

  if(1)
  {
    Sort(C->ev,I,N);

    //===== find translational modes =====
    for(i=0;i<C->N*3;i++)
      for(q=0,b[i]=0.0;q<3;q++)
      {
        for(j=0,dx=0.0;j<C->N;j++)
          dx += C->EV[i][j*3+q];
	b[i] += dx*dx;
      }      
    bm = TOPK(b,C->N*3,3);
    //===== find rotational modes for nanoparticles =====
    for(i=0;i<C->N*3;i++)
      for(q=0,r[i]=0.0;q<3;q++)
      {
        for(j=0,dx=0.0;j<C->N;j++)
	{
	  for(k=0;k<3;k++)
	    x[k] = C->EV[i][j*3+k];
	  VectorProd(x,C->X[j],v);
          dx += v[q];
	}
        r[i] += dx*dx;
      }
    rm = TOPK(r,C->N*3,3);
    printf(" Eigenvalues at the Gamma point (cm-1)\n\n");
    for(i=0;i<C->N*3;i++)
      if( b[I[N-i-1]]>(bm-1e-14) || (C->ND==0&&r[I[N-i-1]]>(rm-1e-14)) )
	printf("\x1B[32m%3d % 24.16lf \x1B[0m\n",i,C->ev[I[N-i-1]]*thzcm/(2.0*Pi));
      else
	printf("%3d % 24.16lf\n",i,C->ev[I[N-i-1]]*thzcm/(2.0*Pi));
  }

  out = fopen("OUTCAR","a");
  fprintf(out,"----------------------------------------------------------------------------------------------------------\n");
  fprintf(out,"                Eigenvalues and eigenvalues at the Gamma point\n");
  fprintf(out,"----------------------------------------------------------------------------------------------------------\n");

  if(0)
  for(i=0;i<3*C->N;i++,printf("\n"))
    for(j=0;j<C->N*3;j++)
      printf("% 1.2lf ",DotProd(C->EV[i],C->EV[j],3*C->N));
  printf("\n");

  for(i=0;i<C->N*3;i++)
  {
    fprintf(out,"eigenvalue %3d % 24.16lf (cm-1)\n",i,C->ev[I[N-i-1]]*thzcm/(2.0*Pi));
    for(j=0;j<C->N;j++,fprintf(out,"\n"))
      for(q=0,fprintf(out,"%3d ",j);q<3;q++)
	fprintf(out,"% 24.16lf",C->EV[I[N-i-1]][j*3+q]);
  }
  fprintf(out,"----------------------------------------------------------------------------------------------------------\n");

  // compute wall and cpu time
  getrusage(RUSAGE_SELF, &t3);
  gettimeofday(&t2, NULL);
  syst_time  = (double) t3.ru_stime.tv_sec + (double) t3.ru_stime.tv_usec / (double) 1e6; // = kernel time
  user_time  = (double) t3.ru_utime.tv_sec + (double) t3.ru_utime.tv_usec / (double) 1e6; // = user time
  wall_time  = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / (double) 1e6;
  cpus_time  = syst_time + user_time;

  fprintf(out,"  \n Total CPU time used (sec):  % 14.8lf     wall time (sec): % 14.8lf",cpus_time,wall_time);
  if(cpus_time/(double)R->NP < 0.5*wall_time)
    fprintf(out,"-\n");
  else
    fprintf(out,"\n");

  fclose(out);

  free_d1D(A);
  free_d1D(B);
  free_d1D(e);
  free_d1D(b);
  free_d1D(r);
  free_i1D(I);
}
//======================================================
void CELL_RELX(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  double H;
  FILE  *out;
  struct timeval t1, t2;
  struct rusage t3;
  double wall_time,user_time,syst_time,cpus_time;

  printf("|                            Cell relaxation                          |\n");
  printf("=======================================================================\n\n");

  if(C->RLXT == 0 && R->MITR > 0)
  {
    printf("Error: type of relaxation RLXT not specified in setup file.\n");
    exit(1);
  }
  
  system("rm -f OUTCAR OSZICAR");

  // initiate the timing
  gettimeofday(&t1, NULL);

  H = CELL_FRC(R,P,W,C,L);
  
  if(R->MITR==0)
    printf("%5d %24.16lf\n",0,H);
  if(C->OUT%10>0)
    CELL_OUT(C);
  if(R->MITR>0)
  {
    if(R->MINT>=0 && R->MINT<4)
    {
      CELL_MIN(R,P,W,C,L);
    }
    else
    {
      fprintf(stderr,"ERROR: Enter valid value for MINT (0-3)\n");
      exit(1);
    }
    
    if( (C->OUT%10==1) )
      CELL_OUT(C);
  }
  if(C->OUT%10==0)
    CELL_OUT(C);
  LIST(C,0);

  if(!CELL_OK(C))
    system("echo >> OUTCAR;echo ERROR distances are too short >> OUTCAR");
  
  // compute wall and cpu time
  getrusage(RUSAGE_SELF, &t3);
  gettimeofday(&t2, NULL);
  syst_time  = (double) t3.ru_stime.tv_sec + (double) t3.ru_stime.tv_usec / (double) 1e6; // = kernel time
  user_time  = (double) t3.ru_utime.tv_sec + (double) t3.ru_utime.tv_usec / (double) 1e6; // = user time
  wall_time  = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / (double) 1e6;
  cpus_time  = syst_time + user_time;

  out=fopen("OUTCAR","a");
  fprintf(out,"  \n Total CPU time used (sec):  % 14.8lf     wall time (sec): % 14.8lf",cpus_time,wall_time);
  if(cpus_time/(double)R->NP < 0.5*wall_time)
    fprintf(out,"-\n");
  else
    fprintf(out,"\n");
  fclose(out);
  
  C->it = 0;
  return;
}
//======================================================
//
//======================================================
void CELL_TEST(ANN *R, PRS *P, PRS *W, Cell *C, LNK *L)
{
  int    i,q,OK;
  double t,dx,Em,Ep,**F;


  printf("|                              Cell test                              |\n");
  printf("=======================================================================\n\n");

  dx = 0.00001;
  F  = make_d2D(C->N,3);

  OK = 1;
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
    {
      t = C->X[i][q];
      C->X[i][q] += dx;
      Ep = CELL_ENE(R,P,W,C,L);
      C->X[i][q] -= 2.0*dx;
      Em = CELL_ENE(R,P,W,C,L);
      F[i][q] = -(Ep-Em)/(dx*2.0);
      C->X[i][q] = t;
    }
  CELL_FRC(R,P,W,C,L);

  for(i=0;i<C->N;i++,printf("\n"))
  {
    for(q=0,printf("NUM  ");q<3;q++)
      printf("% 24.12lf ",F[i][q]);
    printf("\n");
    for(q=0,printf("ANA  ");q<3;q++)
      printf("% 24.12lf ",C->F[i][q]);
    printf("\n");
    for(q=0,printf("DIF  ");q<3;q++)
      printf("% 24.12lf ",C->F[i][q]-F[i][q]);
    printf("\n");
    for(q=0;q<3;q++)
      if( fabs(C->F[i][q]-F[i][q])>1e-6 )
      	OK = 0;
  }
  if( OK==0 )
    printf("CHECK FAILED\n");
  else
    printf("CHECK PASSED\n");
  free_d2D(F,C->N);
}
//======================================================
//
//======================================================
void CELL_MAIN(ANN *R, PRS *P, Cell *C)
{
  int   i;
  PRS   W[9];
  LNK   L;

  if(R->NSPC==0)
  {
    fprintf(stderr,"Error: species are not defined in setup file!\n");
    exit(1);
  }
  Build_Cell(C,1);

  L.B = 0;
  Build_LNK(&L,C->N,C->NM,P->D,3);

  if(R->MODT==1)
  {
    Build_ANN(R);
    READ_ANN(R);
    Build_PRS(P,W,0);
  }
  else
    READ_POT(C,".");

  READ_CELL(C,"POSCAR");
  if( C->ND < 0 )
  {
    fprintf(stderr,"Error: please specify NDIM in setup file\n");
    exit(1);
  }

  P->IO = 0;
  for(i=0;i<C->N;i++)
    L.MRK[i]=1;
  P->EFS = 0;

  JAR(C);
  sprintf(C->ID,"%s",R->ID);
  C->MINT = R->MINT;

  if( R->JOBT==20 ) CELL_RELX(R,P,W,C,&L);
  if( R->JOBT==21 ) CELL_MD(R,P,W,C,&L);
  if( R->JOBT==22 ) CELL_PHON(R,P,W,C,&L);
  if( R->JOBT==23 ) CELL_TEST(R,P,W,C,&L);

  JAR(C);

  SAVE_CELL(C,"CONTCAR",0);
}
//======================================================
//  Convert external code unit cell to MAISE unit cell
//======================================================
int CAST_CELL(Cell *C,  int *ATMN, double *LAT, double *X)
{
  int     i,q,q1;
  double  x[3];

  //===== count # of species =====
  for(i=0,C->NSPC=1;i<C->N;i++)
    if( i>0 && ATMN[i]!=ATMN[i-1] )
      C->NSPC++;

  for(i=0; i < C->NSPC;i++)
    C->SPCN[i] = 0;

  for(i=0;i<D3;i++)
    for(q=0;q<D3;q++)
      C->L[i][q] = LAT[3*i+q]*C->R0;

  //===== assign species & atomic positions =====
  for(i=0;i<C->N;i++)
  {
    C->ATMN[i] = ATMN[i];
    C->ATMZ[i] = C->SPCZ[ATMN[i]];
    C->SPCN[C->ATMN[i]]++;

    for(q=0;q<D3;q++)
      C->X[i][q] = X[3*i+q];

    //===== convert to Cartesian if in fractional =====
    if(C->XT==0)
    {
      for(q=0; q < D3;q++)
      {
	for(q1=0,x[q]=0.0; q1 < D3;q1++)
	  x[q] += C->X[i][q1]*C->L[q1][q];
      }
      for(q=0; q < D3;q++)
	C->X[i][q] = x[q];
    }
  }

  //===== Cartesian =====
  C->XT = 1;

  for(i=0;i<96;i++)
  {
    C->mass[i] = atom_mass[i];
    C->m[i]    = atom_mass[i]*103.6427;
    C->Rm[i]   = mindist[i];
  }

  return 1;  
}
//======================================================
// Interface for calculating energy, forces, stresses
//------------------------------------------------------
// pointers to be passed from external code
//------------------------------------------------------
// ANN    *R     neural network model
// PRS    *P     parser
// PRS    *W     array of parsers
// LNK    *L     link
// Cell   *C     cell
//------------------------------------------------------
// job/cell settings to be defined by user
//------------------------------------------------------
// int    CODE   external code type (for future use)
// int    N      number of atoms
// int    NM     max number of nearest neighbors
// int    ND     clusters (0) or crystals (3)
// int    NP     number of cores for parallelization
// int    XT     fractional (0) or Cartesian (1)
// int    *ATMN  species types
// double *LAT   lattice in 1D array
// double *X     positions in 1D array
//------------------------------------------------------
// calculated values
//------------------------------------------------------ 
// double FRC  forces
// double STR  stresses
//======================================================
double CALL_MAISE(ANN *R, PRS *P, PRS *W, LNK *L, Cell *C, 
		 int CODE, int N, int NM, int ND, int NP, int XT, int *ATMN, double *LAT, double *X, 
                  double *FRC, double *STR)
{
  int     i, q;
  double  H;

  C->N    = N;
  C->NM   = NM;
  C->ND   = ND;
  C->NP   = NP;
  C->XT   = XT;
  C->R0   = 1.0;
  R->JOBT = C->JOBT =    20;
  R->NB   = C->NB   =     1;
  R->A    = C->A    =  C->N;
  R->NP   = C->NP;
  C->p   /= eV2GPa;

  //===== allocate and initialize everything when CALL_MAISE is called the first time =====
  if(L->B == 0)
  {
    P->LM = P->GM = 1;

    READ_MODEL(R,P,C);
    Build_ANN(R);
    if(R->MODT==1)
      READ_ANN(R);
    else
      READ_POT(C,".");

    Build_Cell(C,1);

    C->NSPC = C->nspc = P->NSPC = R->NSPC;
    for(i=0;i<C->nspc;i++)
      C->spcz[i] = C->SPCZ[i] = P->SPCZ[i] = R->SPCZ[i];

    Build_LNK(L,C->N,C->NM,P->D,3);
    for(i=0;i<C->N;i++)
      L->MRK[i] = 1;

    if(R->MODT==1) 
      Build_PRS(P,W,0); 
  }

  //===== convert external code cell to MAISE cell =====
  CAST_CELL(C, ATMN, LAT, X);

  //===== calculate enthalpy, forces, and stresses =====
  H = CELL_FRC(R, P, W, C, L);

  //===== copy forces and stresses to designated arrays =====
  for(i=0;i<C->N;i++)
    for(q=0;q<D3;q++)
      FRC[3*i+q] = C->F[i][q]; 

  for(q=0;q<6;q++)
    STR[q] = C->U[q];

  return H;
}
//======================================================
