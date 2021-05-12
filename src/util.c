#include "util.h"

const double Pi      = 3.1415926535897932384626433832795;
const double sqrt3   = 1.7320508075688772935274463415059;
const double sqrt2   = 1.4142135623730950;
const double au      = 0.529177249;
const double pi      = 3.14159265359;
const double eV2GPa  = 160.2176565;
const int    D3      = 3;
//atomic mass in atomic units
double atom_mass[96] = {0.0,1.008,4.002602,6.94,9.0121831,10.81,12.011,14.007,15.999,18.99840316,20.1797,22.98976928,24.305,26.9815384,28.085,30.973762,32.06,35.45,39.948,39.0983,40.078,44.955908,47.867,50.9415,51.9961,54.938043,55.845,58.933194,58.6934,63.546,65.38,69.723,72.63,74.921595,78.971,79.904,83.798,85.4678,87.62,88.90584,91.224,92.90637,95.95,98,101.07,102.90549,106.42,107.8682,112.414,114.818,118.71,121.76,127.6,126.90447,131.293,132.9054519,137.327,138.90547,140.116,140.90766,144.242,145,150.36,151.964,157.25,158.925354,162.5,164.930328,167.259,168.934218,173.045,174.9668,178.49,180.94788,183.84,186.207,190.23,192.217,195.084,196.96657,200.592,204.38,207.2,208.9804,209,210,222,223,226,227,232.0377,231.03588,238.02891,237,244,243};
// minimum dist. for relaxations
double mindist[96] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.266,1.114,1.0,1.0,1.0,1.0,1.0,1.0,1.610,1.366,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.895,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.979,1.030,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.984,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
// minimum equil. dist. at 000 GPa
double minpres[96] = {0.000,0.400,0.281,1.526,1.105,0.836,0.774,0.907,0.820,0.572,0.582,1.867,1.591,1.428,1.170,1.073,1.053,1.023,1.063,2.353,1.951,1.608,1.438,1.298,1.234,1.234,1.229,1.239,1.241,1.279,1.329,1.223,1.203,1.193,1.199,1.203,1.163,2.524,2.131,1.769,1.597,1.439,1.372,1.363,1.335,1.361,1.399,1.472,1.510,1.510,1.394,1.479,1.394,1.394,1.858,2.743,2.178,1.867,2.045,2.035,2.015,1.995,1.985,1.985,1.965,1.945,1.925,1.925,1.895,1.905,1.875,1.875,1.567,1.439,1.381,1.381,1.349,1.370,1.406,1.475,1.745,1.454,1.464,1.535,1.404,1.504,1.504,2.607,2.216,2.156,2.065,2.005,1.965,1.206,1.875,1.805};
// minimum equil. dist. at 110 GPa
double maxpres[96] = {110.0,0.400,0.203,0.808,0.906,0.733,0.718,0.758,0.820,0.413,0.420,0.901,1.118,1.131,0.941,0.775,0.761,0.739,0.768,0.807,1.302,1.125,1.177,1.124,1.108,1.111,1.101,1.099,1.091,1.139,1.151,0.884,0.869,0.862,1.199,0.869,0.840,0.478,1.041,1.174,0.755,1.240,1.234,1.238,1.222,1.223,1.275,1.298,1.222,1.222,1.000,0.731,0.877,1.007,1.858,1.219,1.192,1.110,1.478,1.471,1.456,1.442,1.435,1.435,1.420,1.406,1.391,1.391,1.369,1.377,1.355,1.355,1.272,1.253,1.247,1.278,1.258,1.258,1.267,1.262,1.236,1.051,1.058,1.343,1.014,1.087,1.087,1.884,1.601,1.558,1.493,1.449,1.420,1.101,1.355,1.304};


//==============================================================================
void PRNT_HEAD()
{
  printf("=======================================================================\n");
  printf("|               Module for Ab Initio Structure Evolution              |\n");
  printf("|                        version     %12s                     |\n",VERSION);
  printf("|                 https://github.com/maise-guide/maise                |\n");
  printf("=======================================================================\n");
}
//==============================================================================
void FPRNT_HEAD(FILE *out)
{
  fprintf(out,"=======================================================================\n");
  fprintf(out,"|               Module for Ab Initio Structure Evolution              |\n");
  fprintf(out,"|                        version     %12s                     |\n",VERSION);
  fprintf(out,"|                 https://github.com/maise-guide/maise                |\n");
  fprintf(out,"=======================================================================\n");
}
//==============================================================================
void EVN(double *B, double *C,double *e, double *b, int N)
{
  double *data;
  gsl_complex A[N*N];
  data = make_d1D(N);
  int i,j;

  for(i=0;i<N*N;i++)
    GSL_SET_COMPLEX(&A[i],B[i],C[i]);

  gsl_matrix_complex *m = gsl_matrix_complex_alloc (N, N);

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      gsl_matrix_complex_set(m,i,j,A[(i*N+j)]);
  
  gsl_vector *eval = gsl_vector_alloc (N);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (N, N);
  gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc (N);
  
  gsl_eigen_hermv (m, eval,evec, w);
  gsl_eigen_hermv_free (w);
  gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  
  for(i=0;i<N;i++)
  {
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
    b[i] = (double) gsl_vector_get(eval,i);
    for(j=0;j<N;j++)
    {
      gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
      e[i*N+j] = (double) GSL_REAL(z);
    }
  }

  gsl_vector_free(eval);
  gsl_matrix_complex_free(evec);
  free_d1D(data);
  return;
}
//==============================================================================
void EV(double A[3][3], double e[3][3], double b[3])
{
  int p,n,k;
  double data[9];

  for(p=k=0; p < 3;p++)
    for(n=0; n < 3;n++)
      data[k++] = A[p][n];

  gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (3);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (3, 3);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (3);

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

  {
    int i, j;

    for(i=0; i < 3;i++)
    {
      gsl_complex eval_i             = gsl_vector_complex_get (eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);

      b[i] = (double)GSL_REAL(eval_i);
      if( fabs(GSL_IMAG(eval_i)) > 1e-14 )
	printf(" OOPS %1.16lf\n",GSL_IMAG(eval_i));

      for(j=0; j < 3;++j)
      {
	gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
	e[i][j] = (double)GSL_REAL(z);
      }
    }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  return;
}
//==============================================================================
void EXIT(char *s)
{
  printf("%s\n",s);
  exit(1);
}
//==============================================================================
//
//==============================================================================
double RANG()
{
  double x1,x2,w;

  w = 2.0;
  while(w >= 1.0)
  {
    x1 = 2.0*Random() - 1.0;
    x2 = 2.0*Random() - 1.0;
    w = x1*x1 + x2*x2;
  }

  w = sqrt( -2.0*log(w)/w );

  return w*x1;
}
//==============================================================================
void VectorProd(double *a, double *b, double *c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
//==============================================================================
double VectorLen(double *c, int D)
{
  int q;
  double a=0.0;

  for(q=0; q < D;q++)
    a += c[q]*c[q];

  return sqrt(a);
}
//==============================================================================
//
//==============================================================================
void VectorNorm(double *c) 
{ 
  double a = 1.0/sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

  c[0] *= a;
  c[1] *= a;
  c[2] *= a;
}
//==============================================================================
double CrossProd(double *a, double *b, double *c)
{
  return a[0]*b[1]*c[2] + a[1]*b[2]*c[0] + a[2]*b[0]*c[1] -
    a[0]*b[2]*c[1] - a[1]*b[0]*c[2] - a[2]*b[1]*c[0];
}
//==============================================================================
double COS(double *x1, double *x2, int D)
{
  int q;
  double r,r1,r2;

  r=r1=r2=0.0;

  for(q=0; q < D;q++)
  {
    r  += x1[q]*x2[q];
    r1 += x1[q]*x1[q];
    r2 += x2[q]*x2[q];
  }    

  return r/sqrt(r1*r2);
}
//==============================================================================
double DiffLen(double *x1, double *x2, int D)
{
  int q;
  double r=0.0;

  for(q=0; q < D;q++)
    r += (x2[q]-x1[q])*(x2[q]-x1[q]);

  return sqrt(r);
}
//==============================================================================
double DotProd(double *x1, double *x2, int D)
{
  int q;
  double r=0.0;

  for(q=0; q < D;q++)
    r += x1[q]*x2[q];

  return r;
}
//==============================================================================
double sign(double a)
{
  if(a>0.0)
    return 1.0;

  return -1.0;
}
//==============================================================================
long TIME(char file[200])
{
  long t;
  FILE *in;

  if( (in = fopen(file,"r"))==0 )
    return -1;
  fscanf(in,"%ld",&t);
  fclose(in);

  return t;
}
//==============================================================================
int TIME_DIF(int t1, int t2)
{
  while( (t1-t2) < 0 )
    t1 += 604800;

  return t1 - t2;
}
//==============================================================================
void Sort(double *x, int *I, int N)
{
  int i,j;

  for(i=0; i < N;i++)
    I[i] = i;

  for(i=0; i < N-1;i++)
    for(j=i+1; j < N;j++)
      if(x[j] < x[I[i]])
        iSwap(&I[i],&I[j]);
}
//==============================================================================
int FIXED(int *vec)
{
  int swap=1,q;

  for(q=0; q < 3;q++)
  {
    swap=swap && vec[q];
  }

  return(swap);
}
//==============================================================================
void inSwap(double D[1000], int I[1000][3], int i, int n)
{
  double t;
  int    j,q;

  t    = D[i];
  D[i] = D[n];
  D[n] = t;

  for(q=0; q < 3;q++)
  {
    j       = I[i][q];
    I[i][q] = I[n][q];
    I[n][q] = j;
  }
}
//==============================================================================
void iSwap(int *i, int *j)
{
  int k;

  k = *i;
  *i = *j;
  *j = k;

  return;
}
//==============================================================================
void dSwap(double *x, double *y)
{
  double t;

  t = *x;
  *x = *y;
  *y = t;

  return;
}
//==============================================================================
// To organize a Q-dimentional "for"-loop
//==============================================================================
void Counter(int *ic, int *N, int i)
{
  ic[i]++;
  if(ic[i]==N[i])
  {
    ic[i] = 0;
    Counter(ic,N,i-1);
  }
}
//==============================================================================
double dMin(double a, double b)
{
  if(a<b)
    return a;

  return b;
}
//==============================================================================
void tprintf(struct Node** tmp_ref, char new_line[200])
{
  struct Node* last = *tmp_ref;
  struct Node* new_node; 

  if( (new_node = (struct Node*) calloc(1,sizeof(struct Node))) == NULL )
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  sprintf(new_node->line,"%s",new_line);
  new_node->next = NULL;

  if(*tmp_ref == NULL)
  {
    *tmp_ref = new_node;
    return;
  }

  while(last->next != NULL)
    last = last->next;

  last->next = new_node;

  return;
}
//==============================================================================
char * tgets(char *buf,int N,struct Node **n)
{
  struct Node* tmp = *n;

  if( *(n) == NULL)
    return NULL;

  if( N == 0 )
    return NULL;

  if( strncpy(buf,tmp->line,N) == NULL )
  {
    return NULL;
  }

  *n=tmp->next;

  return buf;
}
//==============================================================================
void tclose(struct Node **n)
{
  struct Node* tmp = *(n);

  while(tmp)
  {
    *(n) = tmp->next;
    free(tmp);
  }
}
//==============================================================================
//====================  double type memory routines  ===========================
//==============================================================================
double *make_d1D(int x)
{
  double *p;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  return p;
}
//==============================================================================
void free_d1D(double *data)
{
  free(data);
}
//==============================================================================
double **make_d2D(int x, int y)
{
  double **p;
  int    i;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }
  for(i=0; i < x; ++i)   
    p[i] = NULL;

  for(i=0; i < x; ++i)
    if( (p[i] = malloc(y * sizeof *p[i])) == NULL )  
    {
      fprintf(stderr,"Error: memory allocation failed!\n");
      exit(1);
    };

  return p;
}
//==============================================================================
void free_d2D(double **data, int x)
{
  int    i;

  for(i=0; i < x; ++i)
    if( data[i] != NULL )
      free(data[i]);

  free(data);
}
//==============================================================================
double ***make_d3D(int x, int y, int z)
{
  double ***p;
  int    i,j;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }
  for(i=0; i < x; ++i)   
    p[i] = NULL;

  for(i=0; i < x; ++i)
    if( (p[i] = malloc(y * sizeof *p[i])) == NULL )  
    {
      fprintf(stderr,"Error: memory allocation failed!\n");
      exit(1);
    };
  for(i=0; i < x; ++i)    
    for(j=0; j < y; ++j)      
      p[i][j] = NULL;

  for(i=0; i < x; ++i)
    for(j=0; j < y; ++j)
      if( (p[i][j] = malloc(z * sizeof *p[i][j])) == NULL ) 
      {
	fprintf(stderr,"Error: memory allocation failed!\n");
	exit(1);
      }

  return p;
}
//==============================================================================
void free_d3D(double ***data, int x, int y)
{
  int    i,j;

  for (i=0; i < x; ++i)
    if( data[i] != NULL )
    {
      for(j=0; j < y; ++j)
	free(data[i][j]);
      free(data[i]);
    }

  free(data);
}
//==============================================================================
double ****make_d4D(int x, int y, int z,int w)
{
  double ****p;
  int    i,j,k;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  for(i=0; i < x; ++i)   
    p[i] = NULL;

  for(i=0; i < x; ++i)
    if( (p[i] = malloc(y * sizeof *p[i])) == NULL )  
    {
      fprintf(stderr,"Error: memory allocation failed!\n");
      exit(1);
    };

  for(i=0; i < x; ++i)    
    for(j=0; j < y; ++j)      
      p[i][j] = NULL;

  for(i=0; i < x; ++i)
    for(j=0; j < y; ++j)
      if( (p[i][j] = malloc(z * sizeof *p[i][j])) == NULL) 
      {
	fprintf(stderr,"Error: memory allocation failed!\n");
	exit(1);
      };

  for(i=0; i < x; ++i)   
    for(j=0; j < y; ++j)     
      for(k=0; k < z; ++k) p[i][j][k]=NULL;

  for(i=0; i < x; ++i)
    for(j=0; j < y; ++j)
      for(k=0; k < z; ++k)
        if( (p[i][j][k] = malloc(w * sizeof *p[i][j][k])) == NULL) 
	{
	  fprintf(stderr,"Error: memory allocation failed!\n");
	  exit(1);
	};

  return p;
}
//==============================================================================
void free_d4D(double ****data, int x, int y,int z)
{
  int    i,j,k;

  for(i=0; i < x; ++i)
    if(data[i] != NULL)
    {
      for(j=0; j < y; ++j)
	if( data[i][j]!=NULL )
	{
	  for(k=0;k<z;++k)
	    free(data[i][j][k]);
	  free(data[i][j]);
	}

      free(data[i]);
    }

  free(data);
}
//==============================================================================
//=====================  int type memory routines  =============================
//==============================================================================
int *make_i1D(int x)
{
  int    *p;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  return p;
}
//==============================================================================
void free_i1D(int *data)
{
  free(data);
}
//==============================================================================
int **make_i2D(int x, int y)
{
  int    **p;
  int    i;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  for(i=0; i < x; ++i)   
    p[i] = NULL;

  for(i=0; i < x; ++i)
    if( (p[i] = malloc(y * sizeof *p[i])) == NULL )  
    {
      fprintf(stderr,"Error: memory allocation failed!\n");
      exit(1);
    };

  return p;
}
//==============================================================================
void free_i2D(int **data, int x)
{
  int i;

  for (i=0; i < x; ++i)
    if ( data[i] != NULL )
      free(data[i]);

  free(data);
}
//==============================================================================
int ***make_i3D(int x, int y, int z)
{
  int    ***p;
  int    i, j;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  for(i=0; i < x; ++i)   
    p[i] = NULL;

  for(i=0; i < x; ++i)
    if( (p[i] = malloc(y * sizeof *p[i])) == NULL )  
    {
      fprintf(stderr,"Error: memory allocation failed!\n");
      exit(1);
    };

  for(i=0; i < x; ++i)    
    for(j=0; j < y; ++j)      
      p[i][j] = NULL;

  for(i=0; i < x; ++i)
    for(j=0; j < y; ++j)
      if( (p[i][j] = malloc(z * sizeof *p[i][j])) == NULL ) 
      {
	fprintf(stderr,"Error: memory allocation failed!\n");
	exit(1);
      };

  return p;
}
//==============================================================================
void free_i3D(int ***data, int x, int y)
{
  int    i,j;

  for(i=0; i < x; ++i)
    if( data[i] != NULL )
    {
      for(j=0; j < y; ++j)
	free(data[i][j]);

      free(data[i]);
    }

  free(data);
}
//==============================================================================
int ****make_i4D(int x, int y, int z,int w)
{
  int    ****p;
  int    i, j, k;

  if( (p = malloc(x * sizeof *p)) == NULL )  
  {
    fprintf(stderr,"Error: memory allocation failed!\n");
    exit(1);
  }

  for(i=0; i < x; ++i)   
    p[i] = NULL;

  for(i=0; i < x; ++i)
    if( (p[i] = malloc(y * sizeof *p[i])) == NULL )  
    {
      fprintf(stderr,"Error: memory allocation failed!\n");
      exit(1);
    }

  for(i=0; i < x; ++i)    
    for(j=0; j < y; ++j)      
      p[i][j] = NULL;

  for(i=0; i < x; ++i)
    for(j=0; j < y; ++j)
      if( (p[i][j] = malloc(z * sizeof *p[i][j])) == NULL ) 
      {
	fprintf(stderr,"Error: memory allocation failed!\n");
	exit(1);
      }

  for(i=0; i < x; ++i)   
    for(j=0; j < y; ++j)     
      for(k=0; k < z; ++k) 
	p[i][j][k]=NULL;

  for(i=0; i < x; ++i)
    for(j=0; j < y; ++j)
      for(k=0; k < z; ++k)
        if( (p[i][j][k] = malloc(w * sizeof *p[i][j][k])) == NULL ) 
	{
	  fprintf(stderr,"Error: memory allocation failed!\n");
	  exit(1);
	}

  return p;
}
//==============================================================================
void free_i4D(int ****data, int x, int y,int z)
{
  int i, j, k;

  for(i=0; i < x; ++i)
    if( data[i] != NULL )
    {
      for (j=0; j < y; ++j)
	if( data[i][j]!=NULL )
	{
	  for(k=0; k < z;++k)
	    free(data[i][j][k]);

	  free(data[i][j]);
	}

      free(data[i]);
    }

  free(data);
}
//==============================================================================
int READ_CELL(Cell *C, char filename[])
{
  int i,j,k,q,p,q1,TN[10];
  double x[3],L[3][3];
  FILE *in;
  char buf[200],s[100],t[10][100];
  char FF[D3][2]; //True false for moving atom in x, y, and z dirs.
  int tag=0;
  int sp=0;

  if( !(in = fopen(filename,"r")) )
  {
    //perror(filename);
    return 0;
  }

  fgets(buf,200,in);
  buf[strlen(buf)-1] = '\0';
  strcpy(C->TAG,buf);
  fscanf(in,"%lf\n",&C->R0);     

  for(i=0; i < D3;i++)
  {
    for(q=0; q < D3;q++)
      fscanf(in,"%lf ",&L[i][q]);

    fscanf(in,"\n");
  }
  //======reading species; if any
  fgets(buf,81,in);
  sscanf(buf,"%s",s);
  if( !(s[0] >= 48 && s[0] <= 57) )
  {
    tag++;
    sp=sscanf(buf,"%s %s %s %s %s %s %s %s %s %s\n",t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7],t[8],t[9]);
    fgets(buf,81,in);
  }
  //-----------------------------
  sscanf(buf,"%s",s);
  C->NSPC = sscanf(buf,"%d %d %d %d %d %d %d %d %d %d",&TN[0],&TN[1],&TN[2],&TN[3],&TN[4],&TN[5],&TN[6],&TN[7],&TN[8],&TN[9]);

  //===== just read the number of atoms to allocate arrays =====
  if( C->A == 0 )
  {
    for(k=0,C->N=0; k < C->NSPC;k++)
      C->N += TN[k];
    fclose(in);
    C->A = C->N;
    return C->N;
  }
  for(i=0; i < D3;i++)
    for(q=0; q < D3;q++)
      C->L[i][q] = L[i][q]*C->R0;

  if( C->NSPC != sp ) tag=0;

  for(i=0; i < C->NSPC;i++)
    C->SPCN[i] = TN[i];

  for(C->N=0,i=0,j=0; i < C->NSPC;i++)
  {
    C->N += C->SPCN[i];
    for(; j < C->N;j++)
      C->ATMN[j] = i;
  }
  if( C->N > C->A )
  {
    fflush(stdout);
    fprintf(stderr,"Please increase the number of atoms NMAX in setup from %d to %d\n",C->A,C->N);
    exit(1);
  }
  //===assigning species
  if( tag )
  {
    for(i=0; i < C->NSPC;i++)
      C->SPCZ[i] = symb_atom(t[i]);

    for(j=i=0; j < C->NSPC;j++)
      for(k=0; k < C->SPCN[j];k++)
	C->ATMZ[i++]=C->SPCZ[j];
  }
  else
  {
    fflush(stdout);
    fprintf(stderr,"Please specify the type of species in %s !\n",filename);
    exit(1);
  }

  //======================
  fgets(buf,100,in); 
  if( (strncmp(buf,"S",1) == 0) || (strncmp(buf,"s",1) == 0) )
    fgets(buf,100,in);

  sscanf(buf,"%s",buf);

  for(i=0; i < C->N;i++)
  {
    fgets(s,100,in);  
    for(q=0; q < D3;q++)
      strcpy(FF[q],"T\0");

    sscanf(s,"%lf %lf %lf %s %s %s",&C->X[i][0],&C->X[i][1],&C->X[i][2],FF[0],FF[1],FF[2]); 
    if( (strncmp(buf,"d",1) == 0) || (strncmp(buf,"D",1) == 0) )
    {
      for(q=0; q < D3;q++)
      {
	for(q1=0,x[q]=0.0; q1 < D3;q1++)
	  x[q] += C->X[i][q1]*C->L[q1][q];
      }
      for(q=0; q < D3;q++)
	C->X[i][q] = x[q];
    }

    for(q=0; q < D3;q++)
      if( strcmp(FF[q],"T\0")==0 )
	C->FF[i][q] = 1;
      else if (strcmp(FF[q],"F\0")==0 )
        C->FF[i][q] = 0;

    fscanf(in,"\n");
  }

  // read in the velocities if they are listed in the POSCAR file
  p = 0;
  i = 0;
  while ( fgets(s, 100, in) )
  {
    p += sscanf(s,"%lf %lf %lf",&C->V[i][0],&C->V[i][1],&C->V[i][2]);
    i++;
  }
  if( p == 3*C->N )
    C->POS = 1;

  fclose(in);

  C->XT = 1;  // finishes in real coordinates

  if( tag > 0 )
    for(i=0; i < C->N;i++)
    {
      for(k=j=0; k < C->nspc;k++)
	if (C->ATMZ[i] == C->spcz[k] )
	{
	  j++;
	  C->ATMN[i] = k;
	}
      if( j == 0 && C->JOBT/10 > 1 )
      {
        fflush(stdout);
	if( C->JOBT/10 == 3 )
	  fprintf(stderr,"ERROR: species %d in %s is not among TSPC in setup file\n",C->ATMZ[i],filename);
	else
          fprintf(stderr,"ERROR: species %d in %s is not among species in model file\n",C->ATMZ[i],filename);
	exit(1);
      }
    }

  if( C->ND == 0 )
    CENTER(C,0.0);

  // load data corresponding to "table" file
  for(i=0;i<96;i++)
  {
    C->mass[i] = atom_mass[i];
    C->m[i]    = atom_mass[i]*103.6427;
    C->Rm[i]   = mindist[i];
  }

  return 1;  
}
//==============================================================================
void SAVE_CELL(Cell *C, char filename[], int SD) 
{ 
  int i,k,q; 
  FILE *out; 
  char spc[5];

  if( C->ND == 0 )
    CENTER(C,0.5);

  for(i=1; i < C->N;i++)
    if( C->ATMN[i] < C->ATMN[i-1] )
    {
      fprintf(stderr,"Error: cell is not ordered by species!\n");
      exit(1);
    }

  out = fopen(filename,"w"); 
  fprintf(out,"%s\n",C->TAG); 
  fprintf(out,"% 17.14lf\n",C->R0); 
  for(i=0; i < D3;i++) 
  { 
    for(q=0; q < D3;q++) 
      fprintf(out,"% 19.14lf ",C->L[i][q]/C->R0); 

    fprintf(out,"\n"); 
  }

  if( C->SPCZ[0] > 0 )
  {
    for(k=0; k < C->NSPC;k++)
      if( C->SPCN[k] > 0 )
      {
	atom_symb(C->SPCZ[k],spc);
	fprintf(out,"%s ",spc);
      }

    fprintf(out,"\n");
  }


  for(k=0; k < C->NSPC;k++)
    if( C->SPCN[k] > 0 )
      fprintf(out,"%d ",C->SPCN[k]);

  if( SD )
    fprintf(out,"\nSelective dynamics");

  fprintf(out,"\ndirect\n"); 

  Relative(C); 

  for(i=0; i<C->N && SD==0 ; i++)
    for(q=0; q<3; q++)
      if( C->FF[i][q] == 0 )
	SD = 1;

  for(i=0; i < C->N;i++) 
  { 
    for(q=0; q < D3;q++) 
      fprintf(out,"% 19.14lf ",C->X[i][q]); 

    if( SD ) 
    {
      for(q=0; q < 3;q++)
	if( C->FF[i][q] == 0 )
	  fprintf(out,"  F");
	else
	  fprintf(out,"  T");
    }
    fprintf(out,"\n"); 
  } 
  Real(C); 
  if( C->POS == 1 )
    for(i=0,fprintf(out,"\n"); i < C->N;i++)
      fprintf(out,"% 19.14lf % 19.14lf % 19.14lf\n",C->V[i][0],C->V[i][1],C->V[i][2]);

  fclose(out);

  if( C->ND == 0 )
    CENTER(C,0.0);
} 
//==============================================================================
void Read_OUTCAR(Cell *C, char file[], int NC)
{
  int i,j,n,ok,okv;
  double E;
  FILE *in;
  char buf[100],s1[100];

  C->N = 0;

  in = fopen(file,"r");
  ok = okv = 1;

  for(i=0,n=0; i < NC;i++,n++)
  {
    while( ok && i < NC)
    {
      if( fgets(buf,100,in) == 0)
	i = NC;
      else
      {
        sscanf(buf,"%s",s1);

        if( 1 && strcmp(s1,"FREE") == 0 )     // VASP46
        {
          sprintf(s1,"0"); 
	  for(j=0; j < 3;j++)
	    fgets(buf,100,in);

	  fgets(buf,65,in);
	  fgets(buf,100,in);
	  sscanf(buf,"%lf\n",&E);
          fgets(buf,100,in);
        }

        if(1)
          if( strcmp(s1,"POSITION") == 0 )
          {
            sprintf(s1,"0"); 
	    fgets(buf,100,in);
	    for(j=0; j < C->N;j++)
	    {
	      fgets(buf,100,in);
	      sscanf(buf,"%lf %lf %lf %lf %lf %lf\n",
		  &C->X[j][0],&C->X[j][1],&C->X[j][2],&C->F[j][0],&C->F[j][1],&C->F[j][2]);
	    }
	    ok = 0;
          }
      }
    }

    ok = 1;

  }    

  fclose(in);
}
//==============================================================================
double Read_OSZI(char *file)
{
  FILE *in;
  double E;
  char buf[200];

  if( !(in = fopen(file,"r")) )
  {
    perror(file);
    return 1.0;
  }

  while( fgets(buf,100,in) )
    sscanf(buf+27,"%lf\n",&E);
  fclose(in);

  return E;
}
//==============================================================================
void PRNT_LIST(Cell *C)
{
  int i,j,m,n,k,K,N,*I;
  double *R;
  FILE *out;

  out = fopen("list.dat","w");
  for(i=0; i < C->N;i++)
  {
    fprintf(out,"%3d  ",i);
    for(j=0; j < C->Nn[i];j++)
    {
      fprintf(out,"%3d % 2.5lf ",C->Ni[i][j],NDX(C,i,j));
      if( ( (j+1)%12 == 0 && j != C->Nn[i]-1) || C->Nn[i] < 2 )
	fprintf(out,"\n     ");
    }
    fprintf(out,"\n");
  }
  fclose(out);

  out = fopen("bond.dat","w");
  for(i=0; i < C->N;i++)
    for(j=0; j < C->Nn[i];j++)
      fprintf(out,"%5d %3d % 2.12lf\n",i,C->Ni[i][j],NDX(C,i,j));
  fclose(out);

  out = fopen("hist.dat","w");
  for(n=K=0;n<C->NSPC;n++)
    for(m=n;m<C->NSPC;m++)
      for(i=0; i < C->N;i++)
        if(C->ATMN[i]==n)
	  for(j=0; j < C->Nn[i];j++)
	    if(C->ATMN[C->Ni[i][j]]==m)
	      K++;
  I = make_i1D(K);
  R = make_d1D(K);
  for(n=0;n<C->NSPC;n++)
    for(m=n;m<C->NSPC;m++)
    {
      for(i=K=0; i < C->N;i++)
	if(C->ATMN[i]==n)
	for(j=0; j < C->Nn[i];j++)
	  if(C->ATMN[C->Ni[i][j]]==m)
	    R[K++] = NDX(C,i,j);
      Sort(R,I,K);
      if(K>0)
	fprintf(out,"% 12.8lf",R[0]);
      for(k=1,N=1;k<K;k++)
      {
	if( R[I[k]]-R[I[k-1]]<1e-6 )
	  N++;
	else
	{
	  if(n==m)
	    fprintf(out,"%6d   %d %d\n% 12.8lf",N/2,n+1,m+1,R[I[k]]);
	  else
	    fprintf(out,"%6d   %d %d\n% 12.8lf",N/1,n+1,m+1,R[I[k]]);
	  N=1;
	}
      }
      if(K>0)
      {
	if(n==m)
	  fprintf(out,"%6d   %d %d\n",N/2,n+1,m+1);
	else
	  fprintf(out,"%6d   %d %d\n",N/1,n+1,m+1);
      }
    }
  free_d1D(R);
  free_i1D(I);
  fclose(out);

}
//==============================================================================
void Print_LOG(char buf[])
{
  FILE *out;
  out = fopen("log.out","a");
  fprintf(out,"%s",buf);
  fclose(out);
}
//==============================================================================
//    Reads in model parameters from 'model' file
//==============================================================================
void READ_MODEL(ANN *R, PRS *P, Cell *C)
{
  int    i,n,k,ver;
  FILE  *in;
  char   s[225];

  if(R->JOBT/10 == 5)
    sprintf(s,"%s/model",R->otpt);
  else
    sprintf(s,"model");
  
  ver=check_ver(s);
  
  if( (in=fopen(s,"r")) == 0 )
  {
    fprintf(stderr,"ERROR opening model file!\n");
    exit(1);
  }
  R->MODT = 0;
  fgets(s,200,in);
  fgets(s,200,in);
  sscanf(s+2,"%s",s);
  if( strncmp(s,"neural",6) == 0 )
    C->MODT = R->MODT = 1;
  if( strncmp(s,"Gupta",5) == 0 )
    C->MODT = R->MODT = 2;
  if( strncmp(s,"Sutton-Chen",11) == 0 )
    C->MODT = R->MODT = 3;
  if( strncmp(s,"Lennard-Jones",13) == 0 )
    C->MODT = R->MODT = 4;
  
  if( R->MODT==0 )
  {
    fprintf(stderr,"ERROR reading %s\n",s);
    exit(1);
  }
  
  R->NSPC = 0;
  
  while( fgets(s,200,in) )
    if( strncmp(s,"|  number of species    |",25) == 0 )
    {
      sscanf(s+26,"%d" ,&R->NSPC);
      break;
    }
  
  while( fgets(s,200,in) )
    if( strncmp(s,"|  species types        |",25) == 0 )
    {
      for(i=0,n=0,k=26; i < R->NSPC;i++,k+=n)
	sscanf(s+k,"%d%n",&R->SPCZ[i],&n);
      break;
    }
  
  if( R->MODT==1 )
  {
    if( ver < 2400 )
    {
      printf("Error: the model file is in old format; convert the model first!\n");
      exit(1);
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
    
    while( fgets(s,200,in) )
      if( strncmp(s,"  B2A",5) == 0 )
	break;
    for(i=0; i < 4;i++)
      fgets(s,200,in);
    i = -1;
    while( fgets(s,200,in) )
      if( strncmp(s," Rc  ",5) == 0 )
	break;
      else
	i++;
    P->NSYM = R->NSYM = i;
    P->D    = R->NU[0];
    R->D    = R->NU[0];
  }
  fclose(in);


  if( R->MODT == 1 )
    P->DSCR = 2;

  if( R->MODT >  1 )
    P->DSCR = 3;

  return;
}
//==============================================================================
//    reads in EVOS parameters from "setup" file
//==============================================================================
void READ_MAIN(Tribe *T, ANN *R, PRS *P, Cell *C, int J, int ARGC)
{
  int i,n,k,ver;
  FILE *in;
  char buf[200],s[225];
  double t;

  // load data corresponding to "table" file
  for(i=0;i<96;i++)
  {
    C->mass[i]    = atom_mass[i];
    C->m[i]       = atom_mass[i]*103.6427;
    C->Rm[i]      = mindist[i];
    T->maxpres[i] = maxpres[i];
    T->minpres[i] = minpres[i];
  }

  if( ARGC > 1 )
  {
    C->JOBT  = 0;
    T->JOBT  = 0;
    T->N     = 0;
    return;
  }

  if( !(in = fopen("setup","r")) )
  {
    sprintf(buf,"\nPlease either provide setup or use allowed FLAGS (type maise -man for help)\n");
    fprintf(stderr,"\nPlease either provide setup or use FLAGS (type maise -man for help)\n");
    Print_LOG(buf);
    exit(1);
  }

  for(i=0; i < 10;i++)
    T->SPCN[i]=0;

  R->NL    = 0;            // EA quits because NL is not defined and becomes random upon secondary reading of setup
  R->NSPC  = 0;
  R->seed2 = 0;
  T->seed  = 0;
  T->ND    = -1;
  C->ND    = -1;
  R->NP    = 1;
  C->NP    = 1;
  R->NB    = 1;
  C->NB    = 1;
  C->OUT   = 10;           // output EFS at the end of run only
  C->POS   = 0;
  C->A     = 100;
  C->N     = 100;
  R->A     = 100;          //default value of NMAX
  R->MOVI  = 0;
  R->CPLT  = 25.0;
  R->CPLP  = 50.0;  
  R->ICMP  = 0.005;
  R->LREG  = 0.0;
  R->WE    = 1.0;
  R->WF    = 0.01;         // 1 eV/A => 0.010meV/atom
  R->WS    = 0.001;        // 1 kB => ...
  T->te    = 0.0;
  R->ECUT  = 0.9;
  R->EMAX  = 5.0;
  R->WENE  = 1000.0;
  R->FMAX  = 50.0;
  R->FMIN  = 1e-5;
  R->VMIN  = 0.0;
  R->VMAX  = 45.0;
  T->p     = 0.0;          // pressure in GPa
  C->p     = 0.0;          // pressure in eV
  C->Rmax  = 6.0;
  C->Rmin  = 5.8;
  C->DR    = 0.008;
  C->NM    = 250;          // max number of nearest neighbors
  T->No    = 20;           // max number of tries to find the right slice
  T->Ns    = 20;           // max number of random shifts
  T->Nm    = 1000000;      // max number of matings per generation
  T->Nc    = 1000;         // max number of matings per couple
  T->Nu    = 100;          // max number of mutations per couple
  T->NB    = 1;            // number of bests always kept in a generation
  T->Rhc   = 0.80;         // scaling factor for hard-core radius
  T->HE    = 0.10;         // discard p highest enthalpy structures
  T->CUT   = 0.95;         // discard similar structures if CxC is above this
  P->LM    = 1;            // max number of pl for PS descriptor in INI/pl.dat
  P->GM    = 1;            // max number of gn for PS descritpro in INI/w.dat
  R->PENE  = 0;            // parse based on energy; for enthalpy should be 1
  C->RLXT  = 0;
  R->MODT  = 1;
  C->MODT  = 1;
  C->ND    = -1;
  T->ND    = -1;
  T->QT    = 0;            // default queue type is torque
  C->DISP  = 1e-3;         // displacement in Ang for frozen phonon calculation
  strcpy(C->WDIR,".");
  strcpy(R->depo,".");
  strcpy(R->data,".");
  strcpy(R->eval,"~/EVL");
  strcpy(R->otpt,".");

  while( fgets(buf,200,in) )
  {
    //======================  JOB TYPE parameter  ===========================
    if( strncmp(buf,"JOBT",4) == 0 ) { sscanf(buf+4,"%d" , &R->JOBT); T->JOBT = C->JOBT = R->JOBT; }
    if( strncmp(buf,"NPAR",4) == 0 ) { sscanf(buf+4,"%d"  ,&R->NP  );       }
    if( strncmp(buf,"NPAR",4) == 0 ) { sscanf(buf+4,"%d"  ,&C->NP  );       }
    if( strncmp(buf,"NBLK",4) == 0 ) { sscanf(buf+4,"%d"  ,&R->NB  );       }
    //======================  MD parameters  ================================
    if( strncmp(buf,"TMIN",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->TMIN);       }
    if( strncmp(buf,"TMAX",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->TMAX);       }
    if( strncmp(buf,"TSTP",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->TSTP);       }
    if( strncmp(buf,"NSTP",4) == 0 ) { sscanf(buf+4,"%d"  ,&R->NSTP);       }
    if( strncmp(buf,"DELT",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->DELT);       }
    if( strncmp(buf,"CPLT",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->CPLT);       }
    if( strncmp(buf,"CPLP",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->CPLP);       }
    if( strncmp(buf,"ICMP",4) == 0 ) { sscanf(buf+4,"%lf" ,&R->ICMP);       }
    if( strncmp(buf,"MDTP",4) == 0 ) { sscanf(buf+4,"%d"  ,&R->MDTP);       }
    if( strncmp(buf,"MOVI",4) == 0 ) { sscanf(buf+4,"%d"  ,&R->MOVI);       }
    //====================  main structure parameters  ======================
    if( strncmp(buf,"NBLK",4) == 0 ) { sscanf(buf+4,"%d"  ,&C->NB  );       }
    if( strncmp(buf,"COUT",4) == 0 ) { sscanf(buf+4,"%d"  ,&C->OUT );       }
    if( strncmp(buf,"NMAX",4) == 0 ) { sscanf(buf+4,"%d"  ,&C->A   );       }
    if( strncmp(buf,"NMAX",4) == 0 ) { sscanf(buf+4,"%d"  ,&C->N   );       }
    if( strncmp(buf,"NMAX",4) == 0 ) { sscanf(buf+4,"%d"  ,&R->A   );       }
    if( strncmp(buf,"MMAX",4) == 0 ) { sscanf(buf+4,"%d"  ,&C->NM  );       }
    //====================  main EVOS parameters  ===========================
    if( strncmp(buf,"CODE",4) ==0 ) { sscanf(buf+4,"%d" ,&T->CODE  );       }
    if( strncmp(buf,"QUET",4) ==0 ) { sscanf(buf+4,"%d" ,&T->QT    );       }
    if( strncmp(buf,"NDIM",4) ==0 ) { sscanf(buf+4,"%d" ,&T->ND    );       }
    if( strncmp(buf,"NDIM",4) ==0 ) { sscanf(buf+4,"%d" ,&C->ND    );       }
    if( strncmp(buf,"LBOX",4) ==0 ) { sscanf(buf+4,"%lf",&T->B     );       }
    if( strncmp(buf,"NPOP",4) ==0 && J == 1 ) { sscanf(buf+4,"%d" ,&T->N ); }
    if( strncmp(buf,"SITR",4) ==0 && J == 1 ) { sscanf(buf+4,"%d" ,&T->n ); }
    if( strncmp(buf,"NITR",4) == 0          ) { sscanf(buf+4,"%d" ,&T->NI); }
    if( strncmp(buf,"TINI",4) == 0 && J == 1) { sscanf(buf+4,"%d" ,&T->JS); }
    if( strncmp(buf,"RAND",4) == 0) { sscanf(buf+4,"%ld" ,&R->seed2);       }
    if( strncmp(buf,"SEED",4) == 0) { sscanf(buf+4,"%ld" ,&T->seed);        }
    //==================== operation EVOS parameters  =======================
    if( strncmp(buf,"TETR",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 0] = (int)(t*(double)T->N); } sprintf(T->NES[ 0],"TETR");  
    if( strncmp(buf,"PLNT",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 1] = (int)(t*(double)T->N); } sprintf(T->NES[ 1],"PLNT");  
    if( strncmp(buf,"PACK",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 2] = (int)(t*(double)T->N); } sprintf(T->NES[ 2],"PACK");  
    if( strncmp(buf,"BLOB",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 3] = (int)(t*(double)T->N); } sprintf(T->NES[ 3],"BLOB");
    if( strncmp(buf,"MATE",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 4] = (int)(t*(double)T->N); } sprintf(T->NES[ 4],"MATE");  
    if( strncmp(buf,"SWAP",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 5] = (int)(t*(double)T->N); } sprintf(T->NES[ 5],"SWAP");  
    if( strncmp(buf,"RUBE",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 6] = (int)(t*(double)T->N); } sprintf(T->NES[ 6],"RUBE");  
    if( strncmp(buf,"REFL",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 7] = (int)(t*(double)T->N); } sprintf(T->NES[ 7],"REFL");  
    if( strncmp(buf,"INVS",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 8] = (int)(t*(double)T->N); } sprintf(T->NES[ 8],"INVS");  
    if( strncmp(buf,"CHOP",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[ 9] = (int)(t*(double)T->N); } sprintf(T->NES[ 9],"CHOP");  
    if( strncmp(buf,"MUTE",4) == 0 ) { sscanf(buf+4,"%lf",&t ); T->TES[10] = (int)(t*(double)T->N); } sprintf(T->NES[10],"MUTE");  
    //==================== crossover/mutation EVOS parameters  ==============
    if( strncmp(buf,"MCRS",4) == 0 ) { sscanf(buf+4,"%lf",&T->pm  );        }
    if( strncmp(buf,"SCRS",4) == 0 ) { sscanf(buf+4,"%lf",&T->ps  );        }
    if( strncmp(buf,"LCRS",4) == 0 ) { sscanf(buf+4,"%lf",&T->ml  );        }
    if( strncmp(buf,"ACRS",4) == 0 ) { sscanf(buf+4,"%lf",&T->ma  );        }
    if( strncmp(buf,"SDST",4) == 0 ) { sscanf(buf+4,"%lf",&T->cs  );        }
    if( strncmp(buf,"LDST",4) == 0 ) { sscanf(buf+4,"%lf",&T->cl  );        }
    if( strncmp(buf,"ADST",4) == 0 ) { sscanf(buf+4,"%lf",&T->ca  );        }
    if( strncmp(buf,"ELPS",4) == 0 ) { sscanf(buf+4,"%lf",&T->te  );        }
    //==================== species parameters  ==============================
    if( strncmp(buf,"NSPC",4) == 0 ) { sscanf(buf+4,"%d" ,&T->NSPC);        }
    if( strncmp(buf,"NSPC",4) == 0 ) { sscanf(buf+4,"%d" ,&R->NSPC);        }
    if( strncmp(buf,"NSPC",4) == 0 ) { sscanf(buf+4,"%d" ,&P->NSPC);        }
    if( strncmp(buf,"NSPC",4) == 0 ) { sscanf(buf+4,"%d" ,&C->NSPC);        }
    if( strncmp(buf,"NSPC",4) == 0 ) { sscanf(buf+4,"%d" ,&C->nspc);        }
    if( strncmp(buf,"TSPC",4) == 0 ) 
    {
      for(i=0,n=0,k=4; i < T->NSPC;i++,k+=n) sscanf(buf+k,"%d%n",&T->SPCZ[i],&n);
      for(i=0; i < T->NSPC;i++) 
	R->SPCZ[i] = P->SPCZ[i] = C->SPCZ[i] = C->spcz[i] = T->SPCZ[i];

      strcpy(R->compound,"");
      for(i=0; i < R->NSPC;i++)
      {
	atom_symb(R->SPCZ[i],s);
	strcat(R->compound,s);
      }
    }
    if( strncmp(buf,"ASPC",4) == 0 ) 
    {
      for(i=0,n=0,k=4; i < T->NSPC;i++,k+=n) 
	sscanf(buf+k,"%d%n",&T->SPCN[i],&n);

      for(i=0; i < T->NSPC;i++) 
	C->spcn[i]=C->SPCN[i]=T->SPCN[i];
    }
    //======================= I/O parameters  ===============================
    if( strncmp(buf,"DEPO",4) == 0 ) { sscanf(buf+4,"%s",R->depo   );       }
    if( strncmp(buf,"DATA",4) == 0 ) { sscanf(buf+4,"%s",R->data   );       }
    if( strncmp(buf,"EVAL",4) == 0 ) { sscanf(buf+4,"%s",R->eval   );       }
    if( strncmp(buf,"OTPT",4) == 0 ) { sscanf(buf+4,"%s",R->otpt   );       }
    if( strncmp(buf,"WDIR",4) == 0 ) { sscanf(buf+4,"%s",C->WDIR   );       }
    //======================= model parameters  =============================
    if( strncmp(buf,"MODT",4) == 0 ) { sscanf(buf+4,"%d" ,&R->MODT );       }
    if( strncmp(buf,"MODT",4) == 0 ) { sscanf(buf+4,"%d" ,&C->MODT );       }
    //======================= Neural Network model parameters  ==============
    if( strncmp(buf,"NNNN",4) == 0 ) { sscanf(buf+4,"%d" ,&R->NL   );       }
    if( strncmp(buf,"NNNU",4) == 0 ) for(i=1,n=0,k=4; i <= R->NL;i++,k+=n) sscanf(buf+k,"%d%n",&R->NU[i],&n);
    if( strncmp(buf,"NCMP",4) == 0 ) { sscanf(buf+4,"%d" ,&R->D    );       }
    if( strncmp(buf,"NCMP",4) == 0 ) { sscanf(buf+4,"%d" ,&P->D    );       }
    if( strncmp(buf,"NSYM",4) == 0 ) { sscanf(buf+4,"%d" ,&P->NSYM );       }
    if( strncmp(buf,"NSYM",4) == 0 ) { sscanf(buf+4,"%d" ,&R->NSYM );       }
    if( strncmp(buf,"NNGT",4) == 0 ) for(i=1,n=0,k=4; i <= R->NL;i++,k+=n) sscanf(buf+k,"%d%n",&R->GT[i],&n);
    //==================== Neural Network training parameters  ==============
    if( strncmp(buf,"NTRN",4) == 0 ) { sscanf(buf+4,"%d", &R->N    );       }
    if( strncmp(buf,"NTST",4) == 0 ) { sscanf(buf+4,"%d", &R->TN   );       }
    if( strncmp(buf,"TEFS",4) == 0 ) { sscanf(buf+4,"%d" ,&R->EFS  );       }
    if( strncmp(buf,"PEFS",4) == 0 ) { sscanf(buf+4,"%d" ,&P->EFS  );       }
    if( strncmp(buf,"FMRK",4) == 0 ) { sscanf(buf+4,"%lf",&P->FMRK );       }
    if( strncmp(buf,"NNRC",4) == 0 ) { sscanf(buf+4,"%lf",&P->RC   );       }  // for PS only
    if( strncmp(buf,"LREG",4) == 0 ) { sscanf(buf+4,"%lf",&R->LREG );       }
    if( strncmp(buf,"WFRC",4) == 0 ) { sscanf(buf+4,"%lf",&R->WF   );       }
    if( strncmp(buf,"WENE",4) == 0 ) { sscanf(buf+4,"%lf",&R->WENE );       }
    //==================== data parameters  =================================
    if( strncmp(buf,"PENE",4) == 0 ) { sscanf(buf+4,"%d", &R->PENE );       }
    if( strncmp(buf,"EMAX",4) == 0 ) { sscanf(buf+4,"%lf",&R->EMAX );       }
    if( strncmp(buf,"FMAX",4) == 0 ) { sscanf(buf+4,"%lf",&R->FMAX );       }
    if( strncmp(buf,"FMIN",4) == 0 ) { sscanf(buf+4,"%lf",&R->FMIN );       }
    if( strncmp(buf,"VMIN",4) == 0 ) { sscanf(buf+4,"%lf",&R->VMIN );       }
    if( strncmp(buf,"VMAX",4) == 0 ) { sscanf(buf+4,"%lf",&R->VMAX );       }
    //==================== optimization nparameters  ========================
    if( strncmp(buf,"TIME",4) == 0 ) { sscanf(buf+4,"%d" ,&T->time);        }
    if( strncmp(buf,"PGPA",4) == 0 ) { sscanf(buf+4,"%lf",&T->p   );        }
    if( strncmp(buf,"ETOL",4) == 0 ) { sscanf(buf+4,"%lf",&R->ETOL );       }
    if( strncmp(buf,"DISP",4) == 0 ) { sscanf(buf+4,"%lf",&C->DISP );       }
    if( strncmp(buf,"MINT",4) == 0 ) { sscanf(buf+4,"%d" ,&R->MINT );       }
    if( strncmp(buf,"MITR",4) == 0 ) { sscanf(buf+4,"%d" ,&R->MITR );       }
    if( strncmp(buf,"RLXT",4) == 0 ) { sscanf(buf+4,"%d", &C->RLXT );       }
    //===================== rdf parameters ==================================
    if( strncmp(buf,"RMAX",4) == 0 ) { sscanf(buf+4,"%lf",&C->Rmax );       } // Rmax for storing RDF 
    if( strncmp(buf,"RMIN",4) == 0 ) { sscanf(buf+4,"%lf",&C->Rmin );       } // Smooth cutoff between Rmin and Rmax in RDF
    if( strncmp(buf,"DELR",4) == 0 ) { sscanf(buf+4,"%lf",&C->DR   );       } // Gaussian spread 
    //======================= misc tribe parameters  =========================
    if( strncmp(buf,"KMSH",4) == 0 ) { sscanf(buf+4,"%lf",&T->KM  );        }
    if( strncmp(buf,"DENE",4) == 0 ) { sscanf(buf+4,"%lf",&T->DE   );       }
    if( strncmp(buf,"SLIC",4) == 0 ) { sscanf(buf+4,"%d", &T->No   );       }
    if( strncmp(buf,"SHFT",4) == 0 ) { sscanf(buf+4,"%d", &T->Ns   );       }
    if( strncmp(buf,"MAGN",4) == 0 ) { sscanf(buf+4,"%d", &T->Nm   );       }
    if( strncmp(buf,"MACP",4) == 0 ) { sscanf(buf+4,"%d", &T->Nc   );       }
    if( strncmp(buf,"MUCP",4) == 0 ) { sscanf(buf+4,"%d", &T->Nu   );       }
    if( strncmp(buf,"BEST",4) == 0 ) { sscanf(buf+4,"%d", &T->NB   );       }
    if( strncmp(buf,"RHRD",4) == 0 ) { sscanf(buf+4,"%lf",&T->Rhc  );       }
    if( strncmp(buf,"ECUT",4) == 0 ) { sscanf(buf+4,"%lf",&T->HE   );R->ECUT = T->HE;}
    if( strncmp(buf,"SCUT",4) == 0 ) { sscanf(buf+4,"%lf",&T->CUT  );       }
    //==================== misc PARS parameters  ============================
    if( strncmp(buf,"PSLM",4) == 0 ) { sscanf(buf+4,"%d", &P->LM   );       }
    if( strncmp(buf,"PSGM",4) == 0 ) { sscanf(buf+4,"%d", &P->GM   );       }
    if( strncmp(buf,"LNPS",4) == 0 ) { sscanf(buf+4,"%d", &P->LN   );       }
    if( strncmp(buf,"GNPS",4) == 0 ) { sscanf(buf+4,"%d", &P->GN   );       }
    if( strncmp(buf,"NNIO",4) == 0 ) { sscanf(buf+4,"%d", &P->IO   );       }
    //=======================================================================
  }
  
  if( T->seed == 0 ) 
  {
    T->seed = time(NULL);
  }    
  R->seed=T->seed;
  PlantSeeds(T->seed);

  C->p = T->p/eV2GPa;

  if( R->MODT == 1 )
    P->DSCR = 2;

  if( R->MODT >  1 )
    P->DSCR = 3;

  if( P->DSCR == 1 )
    R->D = P->D = P->LN*P->GN;

  if( R->NL > 2 )
  {
    fprintf(stderr,"Number of hidden layers %3d should not exceed 2. Please change NNNN\n",R->NL);
    exit(1);
  }
  
  R->NU[0]        = R->D;
  R->NU[R->NL+1]  = 1;
  R->GT[0]        = 1;
  R->GT[R->NL+1]  = 0;
  R->NL          += 2;

  fclose(in);

  for(n=1,T->SES[0]=T->N,T->FES[0]=T->N+T->TES[0]; n < 11;n++)
  {
    T->SES[n] = T->SES[n-1] + T->TES[n-1];
    T->FES[n] = T->FES[n-1] + T->TES[n  ];
  }
  //===== For INVS the only odd number of atoms allowed is for 1 species in the core =====
  if( T->JOBT == 1 && T->ND == 0 && T->TES[8] > 0 )
  {
    for(n=i=0; i < T->NSPC;i++)
      n += T->SPCN[i]%2;

    if( n > 1 || (T->JS == 0 && (n-T->SPCN[0]%2) > 0) )
    {
      printf("The INVS operation does not work for such number of atoms, resetting INVS to REFL\n");
      fprintf(stderr,"The INVS operation does not work for such number of atoms, resetting INVS to REFL\n");
      T->TES[7] += T->TES[8];
      T->TES[8]  = 0;
    }
  }      
  if( T->JOBT == 1 )
    for(n=0; n < 11;n++)
      if( T->TES[n] > 0 )
        printf("EVLV %3d %s %3d %3d %3d\n",n,T->NES[n],T->TES[n],T->SES[n],T->FES[n]);

  if( T->ND == 0 && T->JOBT == 1 && T->FES[10] != 2*T->N )
  {
    printf("Please make sure that the sum of all operations for nanoparticles is 1\n");
    fprintf(stderr,"Please make sure that the sum of all operations for nanoparticles is 1\n");
    exit(1);
  }

  sprintf(T->VER,VERSION);
  sprintf(R->VER,VERSION);
  sprintf(C->VER,VERSION);

  if( P->DSCR == 3 )
    R->D = P->D = C->NM*C->NSPC;
  
  //===== If the job is a cell operation, this will overwrite the setup and basis =====
  //===== parameters with those in the "model" file in the current directory. =====
  
  if( R->JOBT/10==2 || R->JOBT/10==5 )
  {

    if(R->JOBT/10 == 5)
      sprintf(s,"%s/model",R->otpt);
    else
      sprintf(s,"model");

    ver=check_ver(s);

    if( (in=fopen(s,"r")) == 0 )
    {
      fprintf(stderr,"ERROR opening model file!\n");
      exit(1);
    }
    R->MODT = 0;
    fgets(s,200,in);
    fgets(s,200,in);
    sscanf(s+2,"%s",s);
    if( strncmp(s,"neural",6) == 0 )
      C->MODT = R->MODT = 1;
    if( strncmp(s,"Gupta",5) == 0 )
      C->MODT = R->MODT = 2;
    if( strncmp(s,"Sutton-Chen",11) == 0 )
      C->MODT = R->MODT = 3;
    if( strncmp(s,"Lennard-Jones",13) == 0 )
      C->MODT = R->MODT = 4;

    if( R->MODT==0 )
    {
      fprintf(stderr,"ERROR reading %s\n",s);
      exit(1);
    }

    R->NSPC = 0;

    while( fgets(s,200,in) )
      if( strncmp(s,"|  number of species    |",25) == 0 )
      {
	sscanf(s+26,"%d" ,&R->NSPC);
	break;
      }
    
    while( fgets(s,200,in) )
      if( strncmp(s,"|  species types        |",25) == 0 )
      {
	for(i=0,n=0,k=26; i < R->NSPC;i++,k+=n)
	  sscanf(s+k,"%d%n",&R->SPCZ[i],&n);
	break;
      }
    
    if( R->MODT==1 )
    {
      if( ver < 2400 )
      {
	printf("Error: the model file is in old format; convert the model first!\n");
	exit(1);
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
      
      while( fgets(s,200,in) )
	if( strncmp(s,"  B2A",5) == 0 )
	  break;
      for(i=0; i < 4;i++)
	fgets(s,200,in);
      i = -1;
      while( fgets(s,200,in) )
	if( strncmp(s," Rc  ",5) == 0 )
	  break;
	else
	  i++;
      P->NSYM = R->NSYM = i;
      P->D    = R->NU[0];
      R->D    = R->NU[0];
    }
    fclose(in);

    if( R->JOBT/10==2 )
      {
	C->A = 0;
	READ_CELL(C,"POSCAR");
	R->A = C->A;
	P->NSPC = R->NSPC;
	C->nspc = R->NSPC;
	for(i=0;i<C->nspc;i++)
	  C->spcz[i] = P->SPCZ[i] = R->SPCZ[i];    
      }

    if (R->JOBT/10==5 )
      {
	C->nspc = P->NSPC = R->NSPC;
	for(i=0;i<R->NSPC;i++)
	  C->spcz[i] = C->SPCZ[i] = P->SPCZ[i] = R->SPCZ[i];
      }
  }
}
//==============================================================================
// This is an ANSI C library for multi-stream random number generation.  
// The use of this library is recommended as a replacement for the ANSI C 
// rand() and srand() functions, particularly in simulation applications 
// where the statistical 'goodness' of the random number generator is 
// important.  The library supplies 256 streams of random numbers; use 
// SelectStream(s) to switch between streams indexed s = 0,1,...,255.
//
// The streams must be initialized.  The recommended way to do this is by
// using the function PlantSeeds(x) with the value of x used to initialize 
// the default stream and all other streams initialized automatically with
// values dependent on the value of x.  The following convention is used 
// to initialize the default stream:
//    if x > 0 then x is the state
//    if x < 0 then the state is obtained from the system clock
//    if x = 0 then the state is to be supplied interactively.
//
// The generator used in this library is a so-called 'Lehmer random number
// generator' which returns a pseudo-random number uniformly distributed
// 0.0 and 1.0.  The period is (m - 1) where m = 2,147,483,647 and the
// smallest and largest possible values are (1 / m) and 1 - (1 / m)
// respectively.  For more details see:
// 
//       "Random Number Generators: Good Ones Are Hard To Find"
//                   Steve Park and Keith Miller
//              Communications of the ACM, October 1988
//
// Name            : rngs.c  (Random Number Generation - Multiple Streams)
// Authors         : Steve Park & Dave Geyer
// Language        : ANSI C
// Latest Revision : 09-22-98
//==============================================================================

#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
      
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */

//==============================================================================
double Random(void)
//==============================================================================
// Random returns a pseudo-random real number uniformly distributed 
// between 0.0 and 1.0. 
//==============================================================================
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0) 
    seed[stream] = t;
  else 
    seed[stream] = t + MODULUS;
  return ((double) seed[stream] / MODULUS);
}
//==============================================================================
void PlantSeeds(long x)
//==============================================================================
// Use this function to set the state of all the random number generator 
// streams by "planting" a sequence of states (seeds), one per stream, 
// with all states dictated by the state of the default stream. 
// The sequence of planted states is separated one from the next by 
// 8,367,782 calls to Random().
//==============================================================================
{
  const long Q = MODULUS / A256;
  const long R = MODULUS % A256;
        int  j;
        int  s;

  initialized = 1;
  s = stream;                            /* remember the current stream */
  SelectStream(0);                       /* change to stream 0          */
  PutSeed(x);                            /* set seed[0]                 */
  stream = s;                            /* reset the current stream    */
  for (j = 1; j < STREAMS; j++) {
    x = A256 * (seed[j - 1] % Q) - R * (seed[j - 1] / Q);
    if (x > 0)
      seed[j] = x;
    else
      seed[j] = x + MODULUS;
   }
}
//==============================================================================
void PutSeed(long x)
//==============================================================================
// Use this function to set the state of the current random number 
// generator stream according to the following conventions:
//    if x > 0 then x is the state (unless too large)
//    if x < 0 then the state is obtained from the system clock
//    if x = 0 then the state is to be supplied interactively
//==============================================================================
{
  char ok = 0;

  if (x > 0)
    x = x % MODULUS;                       /* correct if x is too large  */
  if (x < 0)                                 
    x = ((unsigned long) time((time_t *) NULL)) % MODULUS;              
  if (x == 0)                                
    while (!ok) {
      fprintf(stderr,"\nEnter a positive integer seed (9 digits or less) >> ");
      scanf("%ld", &x);
      ok = (0 < x) && (x < MODULUS);
      if (!ok)
        fprintf(stderr,"\nInput out of range ... try again\n");
    }
  seed[stream] = x;
}
//==============================================================================
void GetSeed(long *x)
//==============================================================================
// Use this function to get the state of the current random number 
// generator stream.                                                   
//==============================================================================
{
  *x = seed[stream];
}
//==============================================================================
void SelectStream(int index)
//==============================================================================
//
// Use this function to set the current random number generator
// stream -- that stream from which the next random number will come.
//==============================================================================
{
  stream = ((unsigned int) index) % STREAMS;
  if ((initialized == 0) && (stream != 0))   /* protect against        */
    PlantSeeds(DEFAULT);                     /* un-initialized streams */
}
//==============================================================================
void TestRandom(void)
//==============================================================================
// Use this (optional) function to test for a correct implementation.
//==============================================================================
{
  long   i;
  long   x;
  char   ok = 0;  

  SelectStream(0);                  /* select the default stream */
  PutSeed(1);                       /* and set the state to 1    */
  for(i = 0; i < 10000; i++)
    Random();
  GetSeed(&x);                      /* get the new state value   */
  ok = (x == CHECK);                /* and check for correctness */

  SelectStream(1);                  /* select stream 1                 */ 
  PlantSeeds(1);                    /* set the state of all streams    */
  GetSeed(&x);                      /* get the state of stream 1       */
  ok = ok && (x == A256);           /* x should be the jump multiplier */    
  if (ok)
    printf("\n The implementation of rngs.c is correct.\n\n");
  else
    printf("\n\a ERROR -- the implementation of rngs.c is not correct.\n\n");
}
//==============================================================================
int check_ver(char *fname)
{
  FILE *in;
  char s[200],t[200];
  int n,k,d,i;


  n = k = d = 0;

  if( !(in=fopen(fname,"r")))
    return -2;

  fgets(s,200,in);
  fgets(s,200,in);
  sscanf(s+2,"%s",s);
  if( strncmp(s,"neural",6) != 0 )
  {
    fclose(in);
    return -1;
  }

  while( fgets(s,200,in) )
    if( strncmp(s,"|  maise version",16) == 0 )
      break;

  sscanf(s,"|  maise version | %s |",t);
  i = 6;
  if(i < strlen(t))
  {
    sprintf(s,"%.*s", (int)strcspn(t+i, "."), t+i);
    sscanf(s,"%d",&n);
    i += strlen(s)+1;
    if(i<strlen(t))
    {
      sprintf(s,"%.*s", (int)strcspn(t+i, "."), t+i);
      sscanf(s,"%d",&k);
      i += strlen(s)+1;
      if(i<strlen(t))
      {
	  sprintf(s,"%.*s", (int)strcspn(t+i, "."), t+i);
	  sscanf(s,"%d",&d);
      }
    } 
  }
  
  fclose(in);
  
  return n*1000+k*100+d;
}
