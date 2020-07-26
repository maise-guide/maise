#include "user.h"

//==================================================================
// Examples of cell analysis operations executed with 'maise -usr'
//
//             Explanation of cell variables for Cu3Ag1
//
// C->NSPC        2               number of species
// C->N           4               total number of atoms
// 
//     spc        0   1
// C->SPCZ[spc]  29  47           Z of species spc
// C->SPCN[spc]   3   1           number of atoms by species
//
//      i         0   1   2   3
// C->ATMN[i]     0   0   0   1   species types for each atom
// C->ATMZ[i]    29  29  29  47   species Z     for each atom
// C->Nn[i]      78  78  78  78   number of neighbors within C->Rc (6 Ang)
//
//==================================================================
void USER_CELL(Cell *C, Cell *D, int argc, char **argv)
{
  int    i,NM;
  double A;

  C->A = C->N = 1000;                // set max number of atoms
  NM   = 300;                        // set max number of neighbors
  INIT_CELL(C,"",1,NM,0);            // allocate arrays for C->N and C->NM
  if(READ_CELL(C,"POSCAR")==0)       // read VASP-format structure
  {
    perror("POSCAR");
    exit(0);
  }

  LIST(C,1);                         // find nearest neighbors
  PRNT_LIST(C);                      // print the nearest neighbor list into list.dat
  SAVE_CELL(C,"CONTCAR",0);          // save VASP-format structure

  printf("  i   Ni   Rij_min    Ang_ijk\n");
  for(i=0;i<C->N;i++)
  {
    printf("%3d  %3d ",i,C->Nn[i]);  // print the number of neighbors within C->Rc for atom i
    printf("%9.5lf  ",C->NDX[i][0]); // print the nearest neighbor distance for atom i
    A = acos(Cos(C,i,0,1))*180.0/Pi; // find the angle between bonds to two nearest neighbors
    printf("% 9.5lf\n",A);           // print the angle
  }
}
//==================================================================
