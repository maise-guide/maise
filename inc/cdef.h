#ifndef CDEF
#define CDEF

extern const int D3; // dimensionality of my universe
extern const int linesize;

//=====  Cell parameters  =====

typedef struct
{
  //Lindemann paramateres 
  double Ninv;    // 1.0/(double)C->N 
  double **R1;    // Rij1 for Lindermann index    
  double **R2;    // Rij2 for Lindermann index  
  int    LI;      // # of Lindermann index calls   
  int    *VH;     // velocity histogram     
  int    NH;      // # of bins           
  double VM;      // V max 

  int    A;       // max # of atoms as specified in setup
  int    ND;      // number of periodic dimensions
  int    N;       // # of atomes in the cell
  int    J;       // relaxation type
  int    M;       // number of relaxation DOFs (including lat vecs)

  double R0;      // equilibrium bondlengh
  double dE;      // minimization tolerance

  double **X;     // coordinates
  double **V;     // velocities
  double **F;     // forces
  double *U;      // stresses
  double **W;     // Wyckoff positions
  double ANG[3];
  double LAT[3];
  double m[96];    // atomic mass (adjusted mass)
  double mass[96]; // atomic mass (atomic units)
  int    NRDF;    // number of bins in the RDF
  double ***RDF;  // RDF to calculate CdotC
  double p;       // hydrostatic pressure
  double pint;    // internal pressure

  double P;       // currently enthalpy = E + pV
  double K;       // kinetic energy
  double H;       // to be used in future

  double **L;     // lattice vectors
  double **R;     // reciprical lattice vectors
  int    NM;      // max # of nearest neighbors
  int    *Nn;     // # of nearest neighbors (==max #)
  int    *nn;     // # of NN for rc
  int    MNT;     // max atom type number = 5 (CHECK if it is consistent with NT!)

  int    nspc;    // number of species     (setup)
  int    spcz[10];// Z of each species     (setup)
  int    spcn[10];// atoms of each species (setup)
  int    NSPC;    // number of species     (POSCAR)
  int    SPCZ[10];// Z of each species     (POSCAR)
  int    SPCN[10];// atoms of each species (POSCAR)
  int   *ATMZ;    // Z of individual atoms     
  int   *ATMN;    // ordered species of atoms 0,1,2,... (POSCAR w/r to setup)

  int    **FF;    // tells which atoms are fixed 
  double *BC;     // Bader charge on each atom
  double NE;      // number of electrons
  double E;       // energy
  double *EA;     // energy for each atom
  double *min;
  double *max;
  int    **Ni;    // nearest neighbor number within Rc
  double ***S;    // shifts to handle periodic boundary conditions

  double **e;
  int    *MM;     // initial magnetic moment
  int    NCH[3];
  double dCH[3];
  double *CHG;
  double **NDX;
  double ***DX;

  double Rmax;    // Rmax for CxC
  double Rmin;    // Rmin for CxC
  double DR;      // Gaussian spread for RDF
  double Rm[96];  // min allowed radius distance in INI/atoms
  
  int    SGN;       // space group number
  char   SGS[10];   // space group symbol
  int    NSG;       // # of space group operations
  double **SG;      // space group operations
  int    *SL;       // symbol lengths
  int    NTE;       // # of elements types
  char   ES[9][3];  // element symbols
  char   PRS[3];    // Pearson symbol

  int    NS;        // number of irriduceable sites

  char   TAG[200];
  char   TYPE[10][20];

  int    D;         // dimensionality of descriptor output and NN input
  double *G;
  double ***ndx;    // for storing distances between j and k of atom i
  double ***ndxj;   // for storing distances between i and kth neighbor of C->Ni[i][j]
  double ***cos;    // for storing cosines between ij and ik
  double ***cosj;    // for storing cosines between ij and ik where k is neighbor of C->Ni[i][j]
  double **fc;      // for storing PB cut-off functions
  double Rc;        // cut-off for NN list
  double rc;        // cut-off for a shorter NN list
  int    **Nj;      // Nj[i][j] is neighbor of j corresponding to i

  int    *FRC;      // mask for atom forces to be used in NN training and testing  
  int    RLXT;      // relaxation type after ISIF in VASP
  int    XT;        // type of coordinates: real of fractional
  int    JOBT;      // job type
  int    MODT;      // model type
  int    it;        // iteration number in cell relaxation
  double LJa;       // Lennard-Jones parameter, to be generalized
  double LJe;       // Lennard-Jones parameter, to be generalized 
  double DISP;      // displacement in Ang for frozen phonon calculation

  int    NW;
  double WW[16][50];   // parameters in classical potentials with up to 4 species
  char   WDIR[200];
  int    OUT;       // flags for EFS output
  int    POS;       // flags for VASP output
  int    NP;        // # of threads for NN parallelization over atoms
  int    NB;        // # of blocks  for NN parallelization over atoms
  char   VER[200];  // maise version
  char   ID[200];   // model unique ID
  int    MINT;      // cell optimizer type

  double **EV;      // phonon eigenvectors
  double *ev;       // phonon eigenvalues
  int    EVOK;      // allocation flag
}Cell;

#endif
