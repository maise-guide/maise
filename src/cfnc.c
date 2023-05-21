#include "user.h"
#include "cfnc.h"
#include "eutl.h"

//==================================================================
//  find the dimensionality of the structure: bulk (3) or nano (0)
//==================================================================
int FIND_NDIM(Cell *C)
{
  int i,j,q;

  C->ND = 3;
  LIST(C,0);
  for(i=0;i<C->N;i++)
    for(j=0;j<C->Nn[i];j++)
      for(q=0;q<3;q++)
	if( abs(C->S[i][j][q])>0.01 )
	  return 3; 
  return 0;
}
//==================================================================
void KILL_DBL(Cell *C, double tol)
{
  int i,j,k,q;

  for(i=0; i<C->N; i++)
    for(q=0;q<D3;q++)
    {
      if( C->X[i][q] < 0.0 )
	C->X[i][q] += 1.0;
      
      if(C->X[i][q]>=1.0)
	C->X[i][q] -= 1.0;
    }
  //===== if the diffference is too close to 0.0 or 1.0 =====
  for(i=0;i<C->N;i++)
    for(j=i+1;j<C->N;j++)
      if( fabs(fabs(fabs(C->X[i][0]-C->X[j][0])-0.5)-0.5)<tol && 
	  fabs(fabs(fabs(C->X[i][1]-C->X[j][1])-0.5)-0.5)<tol && 
	  fabs(fabs(fabs(C->X[i][2]-C->X[j][2])-0.5)-0.5)<tol ) 
      {
	for(k=j;k<C->N-1;k++)
	{
	  for(q=0;q<D3;q++)
	    C->X[k][q] = C->X[k+1][q];
          C->ATMN[k] = C->ATMN[k+1];
	  C->ATMZ[k] = C->ATMZ[k+1];
	}
	C->N--;
	j--;
	C->SPCN[C->ATMN[j]]--;
      }
}
//==================================================================
void APPL_SG(Cell *C, double tol)
{
  int    i,j,k,q;

  for(j=0,C->N=0;j<C->NSG;j++)
    for(i=0;i<C->NS;i++)
    {
      if(C->N>C->A)
      {
	fprintf(stderr,"Error: the number of atoms exceeds the allocated number %d \n",C->A);
	exit(1);
      }
      for(q=0;q<D3;q++)
	for(k=0, C->X[C->N][q] = C->SG[j][D3*D3+q] ;k<D3;k++)
	  C->X[C->N][q] += C->W[i][k]*C->SG[j][k+D3*q];
      C->ATMN[C->N] = C->ATMN[i];
      for(q=0;q<D3;q++)
      {
	if(C->X[C->N][q] <0.0)
	  C->X[C->N][q] += 1.0;

	if(C->X[C->N][q]>=1.0)
	  C->X[C->N][q] -= 1.0;
      }
      //=====  note the rescaling with the lattice constants  =====
      for(k=0;k<C->N;k++)
	if( fabs(fabs(fabs(C->X[C->N][0]-C->X[k][0])-0.5)-0.5)<tol/C->LAT[0] && 
	    fabs(fabs(fabs(C->X[C->N][1]-C->X[k][1])-0.5)-0.5)<tol/C->LAT[1] && 
	    fabs(fabs(fabs(C->X[C->N][2]-C->X[k][2])-0.5)-0.5)<tol/C->LAT[2] )
	{
          C->N--;
	  break;
	}
      C->N++;
    }
  ORDER(C);
  C->XT = 0;
  Real(C);
  LIST(C,0);
}
//==================================================================
void READ_CIF(Cell *C, char file[], double tol, int NM, char input[])    
{ 
  int i,j,n,q,k,p; 

  FILE *in; 
  char buf[200],s[200],t[200],r[200],s0[200],s1[200],s2[200]; 

  if(!(in = fopen(file,"r")))
  {
    fprintf(stderr,"Please provide str.cif file\n");
    exit(1);
  }
  //===== for -cif function get the number of atoms to build Cell =====
  if( C->N==0 )
  {
    while(fgets(buf,200,in))
    {
      s[0] = 0;
      sscanf(buf,"%s %s",t,s);
      if(strcmp(t,"_atom_site_fract_z")==0)
      {
	fgets(buf,200,in);
	if(strncmp(buf,"_",1)==0)
	  fgets(buf,200,in);
	
	for(C->NS=0;strncmp(buf,"loop_",5)!=0&&strncmp(buf,"#End",4)!=0;C->NS++)
	{
	  sscanf(buf,"%s %s %d %s %s %s %s %s %s",s,r,&i,r,s0,s1,s2,r,r);
	  C->N += i;
	  fgets(buf,200,in);
	}
      }
    }
    fclose(in);
    if( C->N==0 )
    {
      fprintf(stderr,"str.cif does not have proper Wyckoff positions\n");
      exit(0);    
    }
    else
      INIT_CELL(C,input,4,NM,0);
  }
  //====================================================================
   
  in = fopen(file,"r");
  n = 0; 
  while(fgets(buf,200,in))  
  {
    s[0] = 0;
    sscanf(buf,"%s %s",t,s);

    if(strcmp(t,"_chemical_formula_sum")==0)    // identify element types and quantities
    {
      if(s[0]==0)
	fgets(s,200,in);
      else
	strncpy(s,buf+strlen(t)+1,strlen(buf)-strlen(t));
      s[strlen(s)-2]=0;            // get rid of ' at the end

      for(i=0,C->NTE=0;i<strlen(s);i++) 
        if(s[i] >=65 && s[i] <=90) 
          C->NTE++; 
      C->NSPC = C->NTE;

      for(i=0,j=0;i<strlen(s);j++)  
      {
	C->ES[j][2] = 0;
	sscanf(s+i+1,"%s",t);
	i += strlen(t)+1;
	C->ES[j][0]=t[0];
	C->SL[j] = 2;
	if(t[1] >=97 && t[1] <=122) 
	  C->ES[j][1] = t[1]; 
	else
	{
	  C->ES[j][1] = 0;
	  C->SL[j]--;
	}
	sscanf(t+C->SL[j],"%d",&n);
      }      
      if( strcmp(C->ES[0],"A")==0 )
	C->SPCZ[0] = 0;
      else
	for(i=0;i<C->NSPC;i++)
	  C->SPCZ[i] = symb_atom(C->ES[i]);
    }

    if(strcmp(t,"_symmetry_Int_Tables_number")==0)      // read space group
    { 
      sscanf(s,"%d",&C->SGN); 
      READ_SG(C);
    }

    if(strcmp(t,"_cell_length_a")==0)                    // get unit cell parameters
    {
      sscanf(s,"%lf",&C->LAT[0]);
      C->R0 = 1.0;
    }

    if(strcmp(t,"_cell_length_b")==0) 
      sscanf(s,"%lf",&C->LAT[1]); 
    if(strcmp(t,"_cell_length_c")==0) 
      sscanf(s,"%lf",&C->LAT[2]); 

    if(strcmp(t,"_cell_angle_alpha")==0) 
      sscanf(s,"%lf",&C->ANG[0]); 
    if(strcmp(t,"_cell_angle_beta")==0) 
      sscanf(s,"%lf",&C->ANG[1]); 
    if(strcmp(t,"_cell_angle_gamma")==0) 
    {
      sscanf(s,"%lf",&C->ANG[2]); 
      for(q=0;q<3;q++)
	C->ANG[q] *= Pi/180.0;
      ABC_LT(C);
    }

    if(strcmp(t,"_atom_site_fract_z")==0)              // get Wyckoff positions
    {

      fgets(buf,200,in);  
      if(strncmp(buf,"_",1)==0)
	fgets(buf,200,in);

      for(C->NS=0;strncmp(buf,"loop_",5)!=0&&strncmp(buf,"#End",4)!=0;C->NS++)
      {      
	sscanf(buf,"%s %s %s %s %s %s %s %s %s",s,r,r,r,s0,s1,s2,r,r); 
	sscanf(s0,"%lf",&C->W[C->NS][0]);
        sscanf(s1,"%lf",&C->W[C->NS][1]); 
        sscanf(s2,"%lf",&C->W[C->NS][2]); 

        for(q=0;q<3;q++)
          if(C->W[C->NS][q]<0.0)
            C->W[C->NS][q] += 1.0;
        for(q=0;q<3;q++)
          for(k=0;k<12;k++)
              if(fabs(C->W[C->NS][q]-(double)k/12.0)<0.0002)
                C->W[C->NS][q] = (double)k/12.0;

        fgets(buf,200,in); 

        p = 1;
        if(s[1] >=97 && s[1] <=122)
          p++;
	for(i=0;i<C->NTE;i++)
	  if(p==C->SL[i]&&strncmp(s,C->ES[i],C->SL[i])==0)
	    C->ATMN[C->NS] = i;
      }
    }
  }
  fclose(in);
  APPL_SG(C,tol);
  ORDER(C);
  SAVE_CELL(C,"CONV",0);
}  
//==================================================================
//    find multiplicity of sites for the folded u.c.
//==================================================================
int FIND_MTY(Cell *C, double tol)
{
  int i,j,M;

  JAR(C);
  LIST(C,1);
  for(i=0,M=1;i<C->N;i++)
  {
    for(j=0;j<C->Nn[i];j++)
    {
      if(NDX(C,i,j)>tol)
	break;
      else
	if(C->ATMN[i]!=C->ATMN[C->Ni[i][j]])
	  return 1;
    }
    if(i==0)
      M = j+1;
    if(j+1!=M)
      return 1;
  }

  Relative(C);
  KILL_DBL(C,tol);
  Real(C);

  return M;
}
//==================================================================
//    my algorithm to find the Pearson symbol and PRIM given CONV
//==================================================================
void FIND_PRS(Cell *C, Cell *D, double tol)
{
  int q,k,M;

  abc(C);
  LIST(C,0);
  RDF(C,1);
  LIST(D,0);
  RDF(D,1);

  if(0)
  if(1.0-CxC(C,D)>0.001)
  {
    Copy_C(D,C);
    sprintf(C->PRS,"xX");
    ORDER(C);
    SAVE_CELL(C,"CONV",0);
    SAVE_CELL(C,"PRIM",0);
    return;
  }
  Copy_C(C,D);

  if(   0 <  C->SGN && C->SGN <   3 ) sprintf(C->PRS,"a");
  if(   2 <  C->SGN && C->SGN <  16 ) sprintf(C->PRS,"m");
  if(  15 <  C->SGN && C->SGN <  75 ) sprintf(C->PRS,"o");
  if(  74 <  C->SGN && C->SGN < 143 ) sprintf(C->PRS,"t");
  if( 142 <  C->SGN && C->SGN < 168 ) sprintf(C->PRS,"h");
  if( 167 <  C->SGN && C->SGN < 195 ) sprintf(C->PRS,"h");
  if( 194 <  C->SGN && C->SGN < 231 ) sprintf(C->PRS,"c");

  sprintf(C->PRS+1,"P");  
  M = 1;

  if( (  2 < C->SGN && C->SGN <  16)|| ( 15 < C->SGN && C->SGN <  75) )  // mS or oS
  {
    for(k=0;k<3;k++)
    {
      Copy_C(C,D);
      for(q=0;q<3;q++)
      {
	D->L[0][q] = 0.5*( C->L[(k+0)%3][q]+C->L[(k+1)%3][q]);
        D->L[1][q] = 0.5*(-C->L[(k+0)%3][q]+C->L[(k+1)%3][q]);
        D->L[2][q] = 1.0*( C->L[(k+2)%3][q]);
	if(k>0)
	  dSwap(&D->L[0][q],&D->L[1][q]);
      }
      if((M=FIND_MTY(D,tol))>1)
	break;
    }
    if(M>1)
      sprintf(C->PRS+1,"S");
  }  

  if(M==1)
  if( ( 15 < C->SGN && C->SGN <  75)|| (194 < C->SGN && C->SGN < 231) )  // oF or cF
  {
    Copy_C(C,D);
    for(q=0;q<3;q++)
    {
      D->L[0][q] = 0.5*(C->L[1][q]+C->L[2][q]);
      D->L[1][q] = 0.5*(C->L[2][q]+C->L[0][q]);
      D->L[2][q] = 0.5*(C->L[0][q]+C->L[1][q]);
    }
    M = FIND_MTY(D,tol);
    if(M>1)
      sprintf(C->PRS+1,"F");
  }

  if(M==1)
  if( ( 15 < C->SGN && C->SGN < 143)|| (194 < C->SGN && C->SGN < 231) )  // oI, tI, or cI
  {
    Copy_C(C,D);
    for(q=0;q<3;q++)
    {
      D->L[0][q] = 0.5*(-C->L[0][q]+C->L[1][q]+C->L[2][q]);
      D->L[1][q] = 0.5*( C->L[0][q]-C->L[1][q]+C->L[2][q]);
      D->L[2][q] = 0.5*( C->L[0][q]+C->L[1][q]-C->L[2][q]);
    }
    M = FIND_MTY(D,tol);
    if(M>1)
      sprintf(C->PRS+1,"I");
  }  

  if(M==1)
  if( (142 < C->SGN && C->SGN < 168) )  // hR
  {
    Copy_C(C,D);
    D->L[0][0] =  0.5*C->LAT[0];  D->L[0][1] = 0.5/sqrt(3.0)*C->LAT[0];  D->L[0][2] = C->LAT[2]/3.0;
    D->L[1][0] = -0.5*C->LAT[0];  D->L[1][1] = 0.5/sqrt(3.0)*C->LAT[0];  D->L[1][2] = C->LAT[2]/3.0;
    D->L[2][0] =  0.0*C->LAT[0];  D->L[2][1] =-1.0/sqrt(3.0)*C->LAT[0];  D->L[2][2] = C->LAT[2]/3.0;
    M = FIND_MTY(D,tol);
    if(M>1)
      sprintf(C->PRS+1,"R");
  }

  if(M==1)
    Copy_C(C,D);

  ORDER(D);
  SAVE_CELL(D,"PRIM",0);
}
//==================================================================
//    compute dot product CxC for two structures based on RDF
//==================================================================
void FIND_CXC(Cell *C, Cell *D, int argc, char argv[20][200])
{
  if(argc > 2)
    C->NM   = D->NM   = ( int  )atoi(argv[2]);
  if(argc > 3)
    C->Rmin = D->Rmin = (double)atof(argv[3]);
  if(argc > 4)
    C->Rc = C->Rmax = D->Rc = D->Rmax = (double)atof(argv[4]);
  if(argc > 5)
    C->DR   = D->DR   = (double)atof(argv[5]);

  LIST(C,0);
  LIST(D,0);
  RDF(C,1);
  RDF(D,1);
  Print_RDF_FULL(C,"RDF0.dat");
  Print_RDF_FULL(D,"RDF1.dat");

  printf("%lf\n",CxC(C,D));

}
//==================================================================
//     plot non-normalized RDF 
//==================================================================
void PLOT_RDF(Cell *C, int argc, char argv[20][200])
{
  if(argc > 2)
    C->NM   = ( int  )atoi(argv[2]);
  if(argc > 3)
    C->Rmin = (double)atof(argv[3]);
  if(argc > 4)
    C->Rc = C->Rmax = (double)atof(argv[4]);
  if(argc > 5)
    C->DR   = (double)atof(argv[5]);  

  printf("Max number of nearest neighbors     % d\n",C->NM);
  printf("Soft cutoff for finding neighbors   % lf\n",C->Rmin);
  printf("Hard cutoff for finding neighbors   % lf\n",C->Rmax);
  printf("Gaussian spread for smearing bonds  % lf\n",C->DR);   

  C->ND = FIND_NDIM(C);
  LIST(C,1);
  PRNT_LIST(C);
  printf("\nNeighbor   list      written to     list.dat\n");
  printf(  "Bond       list      written to     bond.dat\n");
  RDF(C,1);
  Print_RDF_FULL(C,"RDF.dat");
  printf(  "Normalized RDF       written to     RDF.dat\n");
  RDF(C,0);
  Print_RDF_FULL(C,"rdf.dat");
  printf(  "Original   RDF       written to     rdf.dat\n");
}
//==================================================================
//    compare RDF, SG, and Volumes for two structures
//==================================================================
void COMP_STR(Cell *C, Cell *D, int argc, char argv[20][200])
{
  int    k;
  double t,o;

  if(argc > 2)
    C->NM   = D->NM   = ( int  )atoi(argv[2]);
  if(argc > 3)
    C->Rmin = D->Rmin = (double)atof(argv[3]);
  if(argc > 4)
    C->Rc = C->Rmax = D->Rc = D->Rmax = (double)atof(argv[4]);
  if(argc > 5)
    C->DR   = D->DR   = (double)atof(argv[5]);

  C->ND = FIND_NDIM(C);
  D->ND = FIND_NDIM(D);

  LIST(C,1);
  PRNT_LIST(C);
  system("mv list.dat list0.dat");
  system("mv bond.dat bond0.dat");
  LIST(D,1);
  PRNT_LIST(D);
  system("mv list.dat list1.dat");
  system("mv bond.dat bond1.dat");

  RDF(C,1);
  RDF(D,1);
  Print_RDF_FULL(C,"RDF0.dat");
  Print_RDF_FULL(D,"RDF1.dat");
  
  printf("\n");
  printf("  STR     vol/atom           space group number           RDF scalar product\n");
  printf(" number   A^3/atom     10^-1   10^-2   10^-4   10^-8         0          1  \n\n");

  t = o = 0.1;
  printf("   %d  % 12.6lf   ",0,CELL_VOL(C)/(double)C->N);
  for(k=1;k<8;k+=2)
  {
    printf("  % 4d  ",FIND_WYC(C,C,t,0));
    t*= o;
    o*= o;
  }
  printf("  % 10.6lf % 10.6lf\n",CxC(C,C),CxC(C,D));
  t = o = 0.1;
  printf("   %d  % 12.6lf   ",1,CELL_VOL(D)/(double)D->N);
  for(k=1;k<8;k+=2)
  {
    printf("  % 4d  ",FIND_WYC(D,D,t,0));
    t*= o;
    o*= o;
  }
  printf("  % 10.6lf % 10.6lf\n",CxC(D,C),CxC(D,D));

  printf("\n");

}
//==================================================================
//    compare RDF, SG, and Volumes for two structures
//==================================================================
void FIND_SPG(Cell *C, Cell *D, double tol, int NM, char input[])
{
  int    n,k,N,M;
  double o;

  if( fabs(tol)>0.5 )
  {
    printf("Specified tolerance magnitude of %lf is too large\n",tol);
    return;
  }

  if(tol>0.0)
  {
    FIND_WYC(C,D,tol,1);
    READ_CIF(C,"str.cif",tol,NM,input);
    FIND_PRS(C,D,tol);
    READ_CELL(D,input);
    printf("%-6d %s%-6d %-6s   %6.1E  % 8.4lf\n",C->SGN,C->PRS,C->N,C->SGS,tol,CxC(C,D));
    return;
  }

  tol = fabs(tol);
  
  N = 0;
  M = FIND_WYC(C,D,1e-12,0);
  for(n=0;n>-13;n--)
    for(k=10;k>0;k--)
    {
      o = (double)k*pow(10.0,(double)n)-1e-12;
      if( o < tol && N!= M)
      {
	READ_CELL(C,input);
	FIND_WYC(C,D,o,1);
	READ_CIF(C,"str.cif",o,NM,input);
	FIND_PRS(C,D,o);
	READ_CELL(D,input); 
	if( N!= C->SGN )
	  printf("%-6d %s%-6d %-6s   %6.1E  % 8.4lf\n",C->SGN,C->PRS,C->N,C->SGS,o,CxC(C,D));
	N = C->SGN;
      }
    }
}
//==================================================================
//   builds and initializes Cell for FLAG operations
//==================================================================
void INIT_CELL(Cell *C, char filename[], int N, int NM, int J)
{
  FILE *in;
  C->MODT = 0;
  C->A    = 0;

  if(J==1)
  {
    if((in=fopen(filename,"r"))==0)
    {
      perror(filename);
      exit(1);
    }

    READ_CELL(C,filename);
  }
    
  C->N   *= N;
  C->A    = C->N;
  C->NM   = NM;
  C->ND   = 3;
  C->Rc   = 6.0;
  C->Rmax = C->Rc;
  C->Rmin = 0.99*C->Rmax;
  C->DR   = 0.008;
  C->NP   = 1;
  C->NB   = 1;

  Build_Cell(C,0);
  if(J==1)
    READ_CELL(C,filename);
}
//==================================================================
//    make supercell
//==================================================================
void MAKE_SUP(Cell *C, Cell *D, int argc, char ARGV[20][200], int NM, char input[])
{
  int    q,n,m,k,N[10],M[10];

  if(argc!=5 && argc!=11)
  {
    printf("Please specify either 3 or 9 integers to generate a supercell\n");
    exit(0);
  }

  for(n=0;n<argc-2;n++)
    N[n] = M[n] = (int)atoi(ARGV[n+2]);

  if(argc==11)
  {
    for(n=0;n<3;n++)
    {
      for(q=0,m=1000,k=-1000;q<3;q++)
      {
	if(M[n+q*3]>k) k = M[n+q*3];
	if(M[n+q*3]<m) m = M[n+q*3];
      }
      N[n] = k-m;
      if(N[n]<=0)
      {
	printf("The specified expansion has no component along axis %d\n",n);
	exit(0);
      }
    }
  }

  INIT_CELL(C,input,N[0]*N[1]*N[2],NM,1);
  INIT_CELL(D,input,1,NM,1);
  Clone(C,D,N[0],N[1],N[2]);
  C->POS = 0;

  if(argc==11)
  {    
    for(n=0;n<3;n++)
      for(q=0;q<3;q++)
	for(m=0,C->L[n][q]=0.0;m<3;m++)
	  C->L[n][q] += D->L[m][q]*(double)M[3*n+m];
    JAR(C);
    Relative(C);
    KILL_DBL(C,1e-10);
    Lat_Align(C);
    Real(C);
  }

  SAVE_CELL(C,"CONTCAR",0);
  exit(0);
}
//======================================================================
// Order atoms by type given an array
//======================================================================
void SWAP_Z(Cell *C, Cell *D, int *Z, int N)
{
  int    i,j,q,n,z[10];
  double O[10];

  //===== order by actual Z =====
  if(N==0)
  {
    for(n=0;n<C->NSPC;n++)
      O[n] = (double)C->SPCZ[n];
    Sort(O,z,C->NSPC);
    for(n=0;n<C->NSPC;n++)
      Z[z[n]] = n;
  }
  //===== order as specified in Z =====
  for(n=0;n<C->NSPC;n++)
    z[n] = 0;
  for(n=0;n<C->NSPC;n++)
    z[Z[n]] = 1;
  for(n=0;n<C->NSPC;n++)
    if( z[n]==0 )
    {
      printf("Please provide %d different numbers between 0 and %d\n",C->NSPC,C->NSPC-1);
      exit(0);
    }
  
  for(n=0;n<C->NSPC;n++)
  {
    C->SPCZ[Z[n]] = D->SPCZ[n];
    C->SPCN[Z[n]] = D->SPCN[n];
  }
  for(i=0;i<C->N;i++)
    C->ATMN[i] = Z[D->ATMN[i]];
  ORDER(C);
  C->POS = 0;
  SAVE_CELL(C,"CONTCAR",0);
}
//==================================================================
//    randomize unit cell
//==================================================================
void RAND_CELL(Cell *C, int argc, char ARGV[20][200])
{
  long   S;
  double dX,dL;

  dX = 0.1;
  dL = 0.1;
  S  = 1;

  if(argc>2)
    dX = (double)atof(ARGV[2]);
  if(argc>3)
    dL = (double)atof(ARGV[3]);
  if(argc>4)
    S  = ( long )atol(ARGV[4]);

  C->POS = 0;
  PlantSeeds(S);
  SHKE_CL(C,dL,dX);

  SAVE_CELL(C,"CONTCAR",0);
  printf("Generated a randomized CONTCAR with dX=%1.4lf dL=%1.4lf seed=%d\n",dX,dL,S);
  exit(0);
}
//==================================================================
//    shift unit cell along a phonon eigenmode
//==================================================================
void EIGM_CELL(Cell *C, int argc, char ARGV[20][200])
{
  int    i,j,q;
  FILE   *in;
  double dX;

  dX = 0.1;

  if(argc>2)
    dX = (double)atof(ARGV[2]);

  C->POS = 0;

  if( !(in=fopen("EV","r")))
  {
    printf("File EV is not found\n");
    exit(0);
  }
  for(i=0;i<C->N;i++)
  {
    fscanf(in,"%lf %lf %lf %lf %lf %lf\n",&C->F[i][0],&C->F[i][1],&C->F[i][2],&C->V[i][0],&C->V[i][1],&C->V[i][2]);
    for(q=0;q<3;q++)
      if( fabs(C->X[i][q]-C->F[i][q])>1e-5)
      {
        printf("EV file is inconsistent with POSCAR");
        exit(0);
      }
    for(q=0;q<3;q++)
      if( fabs(C->V[i][q])<0.05 )
        C->V[i][q] = 0.0;
    if( fabs(C->V[i][0])>1e-5 || fabs(C->V[i][1])>1e-5 ||fabs(C->V[i][2])>1e-5 )
      printf("%3d % lf % lf % lf\n",i,C->V[i][1],C->V[i][1],C->V[i][2]);

  }
  printf("\n");
  fclose(in);

  Relative(C);
  for(i=0;i<C->N;i++)
    for(q=0;q<3;q++)
      C->X[i][q] += C->V[i][q]*dX;

  Real(C);
  JAR(C);

  SAVE_CELL(C,"CONTCAR",0);
  exit(0);
}
//==================================================================
//    run a job defined by FLAGs
//==================================================================
void CELL_EXAM(Cell *C, Cell *D, int argc, char **argv)
{
  double tol,L;
  int    i,NM,N[3],m,M,Z[10];
  char   input[200],input0[200],input1[200],ARGV[20][200];

  NM = 2000;

  //================   list available options  ===============
  if(argc<2)
  {
    printf("For specified JOBT = 0 you should provide a FLAG:\n\n");
    printf("-rdf    compute   and plot the RDF for POSCAR                             \n");
    printf("-cxc    compute   dot product for POSCAR0 and POSCAR1 using RDF           \n");
    printf("-cmp    compare   RDF, space group, and volume of POSCAR0 and POSCAR1     \n");
    printf("-spg    convert   POSCAR into str.cif, CONV, PRIM                         \n");
    printf("-cif    convert   str.cif into CONV and PRIM                              \n");
    printf("-rot    rotate    a nanoparticle along eigenvectors of moments of inertia \n");
    printf("-mov    move      atoms along one direction by a constant shift           \n");
    printf("-ord    order     atoms by species                                        \n");
    printf("-dim    find      whether POSCAR is periodic (3) or non-periodic (0)      \n");
    printf("-box    reset     the box size for nanoparticles                          \n");
    printf("-sup    make      a supercell specified by na x nb x nc                   \n");
    printf("-vol    compute   volume per atom for crystal or nano structures          \n");
    printf("-rnd    randomize unit cell in POSCAR                                     \n");
    printf("-eig    shift     unit cell in POSCAR along eigenmode in EV               \n");
    exit(0);
  }

  printf("|                          Unit cell analysis                         |\n");
  printf("=======================================================================\n\n");

  strcpy(input,"POSCAR");
  strcpy(input0,"POSCAR0");
  strcpy(input1,"POSCAR1");

  if(strncmp(ARGV[1],"-sym",4)==0)
  {
    strcpy(input0,"POSCAR");
    strcpy(input1,"PRIM");    
  }

  for(m=0;m<argc;m++)
    strcpy(ARGV[m],argv[m]);

  //================  use file names for POSCARs   ==========
  if(argc>2)
  {
    M = -1;
    if(strncmp(argv[1],"1",1)==0&&argc>3)
    {
      M = 1;
      strcpy(input,argv[2]);
    }
    if(strncmp(argv[1],"2",1)==0&&argc>4)
    {
      M = 2;
      strcpy(input0,argv[2]);
      strcpy(input1,argv[3]);
    }
    argc -= M+1;
    for(m=1;m<argc;m++)
      strcpy(ARGV[m],argv[m+M+1]);
  }

  //================    plot RDF for 1 structure   ==========
  if(strncmp(ARGV[1],"-rdf",4)==0)
  {
    if(argc>2)
      NM = (int)atoi(ARGV[2]);
    INIT_CELL(C,input,1,NM,1);
    C->POS = 0;
    PLOT_RDF(C,argc,ARGV);
    exit(0);
  }
  //=============  determine Wyckoff positions  =============
  if(strncmp(ARGV[1],"-spg",4)==0)
  {
    tol = 0.01;
    if(argc>2)
      tol = (double)atof(ARGV[2]);
    if(argc>3)
      NM  = (int)atoi(ARGV[3]);
    INIT_CELL(C,input,4,NM,1);
    INIT_CELL(D,input,4,NM,1);
    C->POS = D->POS = 0;
    FIND_SPG(C,D,tol,NM,input);
    exit(0);
  }
  //================  convert cif into CONV  ================
  if(strncmp(ARGV[1],"-cif",4)==0)
  {
    tol = 0.01;
    if(argc>2)
      tol = (double)atof(ARGV[2]);
    if(argc>3)
      NM  = (int)atoi(ARGV[3]);
    C->N = 0;
    C->POS = 0;
    READ_CIF(C,"str.cif",tol,NM,input);    
    exit(0);
  }
  //================  find RDF-based dot product   ==========
  if(strncmp(ARGV[1],"-cxc",4)==0)
  {
    if(argc>2)
      NM = (int)atoi(ARGV[2]);
    INIT_CELL(C,input0,1,NM,1);
    INIT_CELL(D,input1,1,NM,1);
    C->POS = D->POS = 0;
    FIND_CXC(C,D,argc,ARGV);
    exit(0);
  }
  //================  compare RDF of 2 structures  ==========
  if(strncmp(ARGV[1],"-cmp",4)==0)
  {
    if(argc>2)
      NM = (int)atoi(ARGV[2]);
    INIT_CELL(C,input0,4,NM,1);
    INIT_CELL(D,input1,4,NM,1);
    C->POS = D->POS = 0;
    COMP_STR(C,D,argc,ARGV);
    exit(0);
  }
  //================  rotate a nanoparticle  ================
  if(strncmp(ARGV[1],"-rot",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    C->ND = FIND_NDIM(C);
    CENTER(C,0.0);
    if( C->ND!= 0)
    {
      fprintf(stderr,"POSCAR is periodic\n");
      exit(1);
    }
    NANO_ROT(C,1);
    exit(0);
  }
  //================ move atoms by a shift   ================
  if(strncmp(ARGV[1],"-mov",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    if(argc>2)
      m  = (int)atoi(ARGV[2]);
    if(argc>3)
      L = (double)atof(ARGV[3]);
    Relative(C);
    for(i=0;i<C->N;i++)
      C->X[i][m] += L;
    Real(C);
    JAR(C);
    C->POS = 0;
    SAVE_CELL(C,"CONTCAR",0);
    exit(0);
  }
  //================ order atoms by species  ================
  if(strncmp(ARGV[1],"-ord",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    INIT_CELL(D,input,1,NM,1);
    if( argc != 2 && ((argc-2) != C->NSPC) )
    {
      printf("Please provide %d numbers between 0 and %d to reorder species in %s\n",C->NSPC,C->NSPC-1,input);
      exit(0);
    }
    for(i=0;i<C->NSPC;i++)
      Z[i] = (int)atoi(ARGV[i+2]);
    SWAP_Z(C,D,Z,argc-2);
    exit(0);
  }
  //================  determine dimensionality  =============
  if(strncmp(ARGV[1],"-dim",4)==0)
  {
    if(argc>2)
      C->Rc = (double)atof(ARGV[2]);
    INIT_CELL(C,input,1,NM,1);
    printf("% 3d\n",FIND_NDIM(C));
    exit(0);
  }
  //================  resize the nanparticle box  ===========
  if(strncmp(ARGV[1],"-box",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    C->ND = FIND_NDIM(C);
    if( C->ND!= 0)
      printf("WARNING POSCAR is periodic\n");

    if( argc<3 )
    {
      fprintf(stderr,"Please provide new box size\n");
      exit(0);
    }
    L = (double)atof(ARGV[2]);
    for(i=0;i<3;i++)
      C->L[i][i] = L;
    SAVE_CELL(C,"CONTCAR",0);
    exit(0);
  }
  //===================== make a supercell  ==================
  if(strncmp(ARGV[1],"-sup",4)==0)
  {
    MAKE_SUP(C,D,argc,ARGV,NM,input);
    exit(0);
  }
  //================  compute volume per atom  ================
  if(strncmp(ARGV[1],"-vol",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    C->ND = FIND_NDIM(C);
    printf("% lf\n",CELL_VOL(C)/(double)C->N);
    exit(0);
  }
  //==================  randomize unit cell  ==================
  if(strncmp(ARGV[1],"-rnd",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    RAND_CELL(C,argc,ARGV);
    exit(0);
  }
  //==================  randomize unit cell  ==================
  if(strncmp(ARGV[1],"-eig",4)==0)
  {
    INIT_CELL(C,input,1,NM,1);
    EIGM_CELL(C,argc,ARGV);
    exit(0);
  }
  //================  run user-defined functions   ============
  if(strncmp(ARGV[1],"-usr",4)==0)  
  {
    USER_CELL(C,D,argc,argv);
    exit(0);
  }

  //================   manual for flag operations =============
  if(strncmp(ARGV[1],"-man",4)==0)
  {
    printf("-rdf    compute and plot the RDF for POSCAR                             \n");
    printf("        options: [ max near. neigh., soft cutoff, hard cutoff, spread ] \n");
    printf("-cxc    compute dot product for POSCAR0 and POSCAR1 using RDF           \n");
    printf("        options: [ max near. neigh., soft cutoff, hard cutoff, spread ] \n");
    printf("-cmp    compare RDF, space group, and volume of POSCAR0 and POSCAR1     \n");
    printf("        options: [ max near. neigh., soft cutoff, hard cutoff, spread ] \n");
    printf("-spg    convert POSCAR into str.cif, CONV, PRIM                         \n");
    printf("        options: [ tolerance, max near. neigh. ]                        \n");
    printf("-cif    convert str.cif into CONV and PRIM                              \n");
    printf("        options: [ tolerance, max near. neigh. ]                        \n");
    printf("-rot    rotate  a nanoparticle along eigenvectors of moments of inertia \n");
    printf("-mov    move atoms along one direction by a constant shift              \n");
    printf("-ord    order atoms by species                                          \n");
    printf("-dim    find    whether POSCAR is periodic (3) or non-periodic (0)      \n");
    printf("-box    reset   the box size for nanoparticles                          \n");
    printf("-sup    make    a supercell specified by na x nb x nc                   \n");
    printf("-vol    compute volume per atom for crystal or nano structures          \n");
    printf("-rdf    randomize unit cell in POSCAR                                   \n");
    printf("        options: [ at disortions, uc distortions, seed ]                \n");
    printf("-eig    shift unit cell in POSCAR along eigenmode in EV                 \n");
    printf("        options: [ distortion strength ]                                \n");
    printf("-usr    run user-defined functions                                      \n");    
    exit(0);    
  }

  //================   list available options  ===============
  if(strncmp(ARGV[1],"-man",4)!=0)
    printf("\nThe FLAG is   not recognized. Allowed FLAGS are:\n\n");
  printf("-rdf    compute   and plot the RDF for POSCAR                             \n");
  printf("-cxc    compute   dot product for POSCAR0 and POSCAR1 using RDF           \n");
  printf("-cmp    compare   RDF, space group, and volume of POSCAR0 and POSCAR1     \n");
  printf("-spg    convert   POSCAR into str.cif, CONV, PRIM                         \n");
  printf("-cif    convert   str.cif into CONV and PRIM                              \n");
  printf("-rot    rotate    a nanoparticle along eigenvectors of moments of inertia \n");
  printf("-mov    move      atoms along one direction by a constant shift           \n");
  printf("-ord    order     atoms by species                                        \n");
  printf("-dim    find      whether POSCAR is periodic (3) or non-periodic (0)      \n");
  printf("-box    reset     the box size for nanoparticles                          \n");
  printf("-sup    make      a supercell specified by na x nb x nc                   \n");
  printf("-vol    compute   volume per atom for crystal or nano structures          \n");
  printf("-rnd    randomize unit cell in POSCAR                                     \n");
  printf("-eig    shift     unit cell in POSCAR along eigenmode in EV               \n");
  printf("-usr    run       user-defined functions                                  \n");
  exit(0);
}
//==================================================================
