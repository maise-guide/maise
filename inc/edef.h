#ifndef EDEF
#define EDEF

// Tribe parameters
typedef struct
{
  int    JOBT;    // job type
  int    CODE;    // code to use in ES
  int    ND;      // system type 3=crystal 2=film 0=particle
  double B;       // box size for particles
  int    N;       // population size
  int    M;       // number of children
  int    NI;      // number of generations
  int    P;       // number of parents (N-M)
  int    NB;      // number of bests always kept in a generation

  int    TES[11];    // type   of evolutionary search operation
  char   NES[11][5]; // name   of evolutionary search operation
  int    SES[11];    // srart  of evolutionary search operation
  int    FES[11];    // finish of evolutionary search operation

  int    QT;      // queuing system type: 0=torque 1=slurm
  int    JS;      // start type at n=0: 0=random 1=specified 2=specified cell 3=spec. cell and fixed coords
  int    n;       // generation number
  int    time;    // max time per relaxation
  long   seed;    // starting seed

  Cell   *C;      // pointer to tribe's cells
  Cell   *Z;      // pointer to tribe's cells to deal with clusters
  Cell   CC;      // cell for generic purposes

  double R0;      // equilibrium bondlengh
  double dE;      // minimization tolerance
  double DE;      // energy/atom window w.r.t. lowest energy to collect all distinct structures in POOL at the end
  double VOL;     // total volume of hard spheres of Rm

  int    NSPC;      // number of types of species (MAX=5)
  int    SPCZ[10];   // species type
  int    SPCN[10];   // number of each type
  char   SS[10][3];// species symbol
  int    MM[5];    // starting magnetic moment: if negative, random sign
  int    MAG;     // magnetic (1) or non-magnetic (0)

  double **X;     // coordinates
  double **V;     // velocities
  double **F;     // forces
  double *m;      // mass

  int    No;      // max number of tries to find the right slice
  int    Ns;      // max number of random shifts
  int    Nm;      // max number of matings per generation
  int    Nc;      // max number of matings per couple
  int    Nu;      // max number of mutations per couple
  double cs;      // cloning: swapping rate
  double cl;      // cloning: mutation strength for lattice vectors
  double ca;      // cloning: mutation strength for atomic potisions
  double pm;      // mutation rate
  double ps;      // mating:  swapping rate
  double ml;      // mating:  mutation strength for lattice vectors
  double ma;      // mating:  mutation strength for atomic potisions
  double te;      // tetris:  nanoparticle ellipticity
  int    pos;     // plant:   number of POSCAR files in INI to seed

  double p;       // pressure
  double KM;      // k-mesh
  double *E;      // current energy 
  double *P1;     // energy of parent1
  double *P2;     // energy of parent2
  double **T;     // all energies
  double *f;      // fitness
  double Rm[5];   // min radius at given P
  double Rhc;     // scale for hard-core radius
  double HE;      // discard highest enthalpy structures
  double CUT;     // similarity cutoff (if > CUT (~0.9) declare similar)
  int    *S;      // position of the selected members in the population
  int    **I;     // number in the full population
  int    **G;     // parents number
  int    **Ni;    // 
  char   VER[200];// version
  char   ISO[200];// path to isotropy
  double maxpres[96];  // min. equi. dist. at 000 GPa
  double minpres[96];  // min. equi. dist. at 110 GPa
  int    MODE;    // operation mode for manual EVOS (JOBT=14)
}Tribe;
#endif
