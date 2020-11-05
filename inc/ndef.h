#ifndef NDEF
#define NDEF

typedef struct
{
  char   ID[200];    // model unique ID
  int    A;          // MAX # of atoms in all structure
  int    STR;        // # of training and testing structures
  int    N;          // # of training data structures
  int    TN;         // # of test data structures
  int    M;          // # of ANN (gaussian) centers
  int    D;          // # of dimensions in configuration space 
  int    O;          // starting index for unconstrained inputs in stratified training
  int    NN;         // # of first nearest neighbors
  int    Q;          // # of output variables
  int    NE;
  int    NF;         // # of forces in the training dataset
  int    TNF;        // # of forces in the testing  dataset
  double *E;         // energy
  double *Fc;        // force
  double *R1;        // input vector
  double **C;        // ANN centers
  double **F;        // matrix (to be solved: F * W = T)
  double **U;        // not necessary, could use F instead
  double **V;        // V (NOT transposed!) matrix in SVD
  double *G;         // weight of input variables
  double **Rmin;     // lower and upper limits 
  double **Rmax;     // of input variables
  double **DR;       // 1.0/(Rmax-Rmin)
  int    JOBT;
  int    MODT;
  int    DSCR;
  double RE;         // residual error for E
  double RF;         // residual error for F
  double RS;         // residual error for S
  double RT;         // residual error for EFS (total)
  double EE;         // test     error for E
  double EF;         // test     error for F
  double ES;         // test     error for S
  double ET;         // test     error for EFS (total)
  int    *train;
  int    *test;

  int    DNm;     // max # of nearest neighbors in training set
  double Eavg;
  double Edev;
  double *TEeval;
  double WE;      // weights for EFS
  double WF;
  double WS;
  double WENE;    // energy range for calculating error weights
  double time;
  int    EFS;     // E (0) EF (1) ES (2) EFS (3)
  //===== for data analysis =====
  double ECUT;
  double EMAX;
  double FMAX;
  double FMIN;
  double VMIN;
  double VMAX;
  char   depo[200];
  char   data[200];
  char   eval[200];
  char   otpt[200];
  //====== MLP parameters ======  
  int    NW;
  int    NL;
  int    NU[4];    // maxiumum # of hiddel layers is 2
  int    GT[4];    // maxiumum # of hiddel layers is 2
  double ***B;
  double ****W;    // weights
  double ****e;    // errors  in each layer
  double ****d;    // derivatives of errors
  double ****c;    // curvatures  of errors
  double ***Bp;
  double ****Wp;    
  double **WW;     // parameters for classical potentials
  double **WWp;    // parameter derivatives
  double **We;     // temporary array for storing rho_i
  int    MINT;     // minimizer type
  int    MITR;     // max number of iterations
  double ETOL;     // error tolerance for convergence
  double LREG;     // regularization parameter lamba
  int    seed;
  int    ITER;
  char   file_name[200];  //for extra test in evaluation part, should get from file "test.dat"
  char   test_add[200]; //for extra test in evaluation part, should get from file "test.dat"

  int    NSPC;     // number of species     (setup)
  int    SPCZ[10]; // Z of each species     (setup)

  int    NSYM;     // number of symmetry functions
  char   compound[200];
  long   seed2;
  int    N0;     // max # for W, Wp index 0
  int    N1;     // max # for W, Wp index 1
  int    NP;     // # of threads for NN parallelization over atoms
  int    NB;        // # of blocks  for NN parallelization over atoms
  //======= to mix-training
  char   file_list[11][200];   //10=MAX num of elements in stretified training; also in disc.c reading MLPS hardcoded
  int    **mask; //to define mask for cmp: 0=fixed and 1=adjustable
  int    nw;
  int    MIX;
  //====== MD parameters                                                                                                             
  int    NSTP;
  double TMIN;
  double TMAX;
  double TSTP;
  double DELT;
  double CPLT; //coupling constant: for thermostat
  double CPLP;
  double ICMP; //isothermal compressibility: for barostat-> effective: ICMP/COPL
  int    MDTP; // type of MD run
  int    MOVI;
  char   VER[200]; // version
  double Rc;
  int    PENE;    // 0 for parse based on energy; 1 for enthalpy
}ANN;

typedef struct
{
  int    N;
  int    *ATMN;  //vector to keep species of atoms codes: 0,1,2,...
  int    NM;
  int    NF;
  int    *DNn;
  int    **DNi;
  int    **DNs;
  double E;
  double e;     // to keep E in MLP training
  double **F;
  double **f;   // to keep F in MLP training 
  double *S;
  double *EA;   // energy of each atom
  double *s;    // to keep S in MLP training
  double **Cn;
  double ****Fn;
  double ***Sn;
  double p;     // hydrostatic pressure
  int    *MRK;
  int    *Fi;
  int    B;     // built or not yet
  char   path[200]; // path to POSCAR.0
  double DE;    // energy above the lowest-energy structure in the subset
  double W;     // weight for each structure in the error function
} LNK;

typedef struct
{
  int    B;
  int    D;
  double **Pl;    // power spectrum coefficients             
  int    LN;      // lmax                                   
  int    LM;      // LMAX tabulated terms               
  double RC;      // Rcut for Cn
  double **GW;    // power spectrum coefficients         
  int    GM;      // GMax        
  int    GN;      // G            
  int    *GT;     // max # of different PB paramters
  int    **GF;    // PB descriptor functions    
  double **GP;    // PB descriptor parameters 
  double Rc;      // PB cut-off distance (if only one)    
  int    IO;      // 0: parse  1: parse/write 2: read
  int    EFS;     // E (0) EF (1) ES (2) EFS (3)
  int    DSCR;    // descriptor type
  double FMRK;    // fraction of atoms to mark for FRC training
  /////////// for species
  int    NG2;     // actual number of G2 functions extracted from GN.dat file
  int    NG4;     // actual number of G4 functions extracted from GN.dat file
  int    NSYM;    // total number of symmetry functions in GN.dat file: read from "setup"
  int    NSPC;    // number of species     (setup)
  int    SPCZ[10];// Z of each species     (setup)
} PRS;

typedef struct
{
  double  ***Bp;
  double ****Wp;
  double  ***e;
  double  ***d;
  double  ***c;
  double  ***xn;
  double  ***xm;
} PAR;
#endif
