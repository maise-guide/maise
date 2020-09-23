#include "mlib.h"

int main(int argc, char **argv)
{
  int    i, q;

  //===== pointers to MAISE structures ======
  ANN    R;
  PRS    P; 
  PRS    W[9];
  LNK    L;
  Cell   C;

  //======== cell variables to pass =========
  int    CODE;                                 // external code type
  int    N;                                    // number of atoms
  int    NM;                                   // max number of nearest neighbors
  int    ND;                                   // clusters (0) or crystals (3)
  int    NP;                                   // number of cores for parallelization
  int    XT;                                   // fractional (0) or Cartesian (1)
  int    ATMN[4];                              // int y int yspecies types
  double LAT[3][3];                            // lattice vectors
  double POS[4][3];                            // positions

  //== returned enthalpy, forces, stresses ==
  double H;
  double FRC[4][3];
  double STR[6];

  double LATf[3*3];                            // lattice vectors flattened 1d array
  double POSf[4*3];                            //       positions flattened 1d array
  double FRCf[4*3];                            //          forces flattened 1d array

  CODE =  0;                                   // currently 0 for all code types
  N    =  4; 
  NM   =  500; 
  ND   =  0; 
  NP   =  4; 
  XT   =  1;
  C.p  =  0.0;                                 // pressure in GPa for ND = 3
  
  LAT[ 0][ 0] = 20.0; LAT[ 0][ 1] =  0.0; LAT[ 0][ 2] =  0.0;
  LAT[ 1][ 0] =  0.0; LAT[ 1][ 1] = 20.0; LAT[ 1][ 2] =  0.0;
  LAT[ 2][ 0] =  0.0; LAT[ 2][ 1] =  0.0; LAT[ 2][ 2] = 20.0;

  POS[ 0][ 0] =  2.0; POS[ 0][ 1] =  0.0; POS[ 0][ 2] =  0.0;
  POS[ 1][ 0] = -2.0; POS[ 1][ 1] =  0.0; POS[ 1][ 2] =  0.0;
  POS[ 2][ 0] =  0.0; POS[ 2][ 1] =  2.0; POS[ 2][ 2] =  0.0;
  POS[ 3][ 0] =  0.0; POS[ 3][ 1] = -2.0; POS[ 3][ 2] =  0.0;

  ATMN[0] = 0; ATMN[1] = 0; ATMN[2] = 0; ATMN[3] = 0;

  L.B = 0;

  for(i=0;i<3;i++)
    for(q=0;q<3;q++)
      LATf[3*i+q] = LAT[i][q];

  for(i=0;i<N;i++)
    for(q=0;q<3;q++)
      POSf[3*i+q] = POS[i][q];

  H = CALL_MAISE(&R, &P, W, &L, &C, CODE, N, NM, ND, NP, XT, ATMN, LATf, POSf, FRCf, STR);

  for(i=0;i<N;i++)
    for(q=0;q<3;q++)
      FRC[i][q] = FRCf[3*i+q];

  printf("Force components\n");
  for(i=0;i<N;i++)
    printf("%3d  %18.14lf  %18.14lf  %18.14lf\n", i, FRC[i][0], FRC[i][1], FRC[i][2]);

  printf("Enthalpy  %18.14lf\n\n", H);
  printf("The forces and enthalpy should be consistent with the values in ./test/OUTCAR.\n\n");

  return 0;
}
