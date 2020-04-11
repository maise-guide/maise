//========================================================================//
//                    GNU GENERAL PUBLIC LICENSE                          //
//                      Version 3, 29 June 2007                           //
//                                                                        //
//   Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>  //
//     Everyone is permitted to copy and distribute verbatim copies       //
//      of this license document, but changing it is not allowed.         //
//========================================================================//

#include "mlib.h"

int main(int argc, char** argv)
{
  Tribe T;
  ANN   R;
  PRS   P;
  Cell  C,D;

  //========== print the header of maise ===========
  PRNT_HEAD();

  //=====  read setup, atoms, and table files  =====
  READ_MAIN(&T,&R,&P,&C,1,argc);

  //========  run structure analysis jobs  =========
  if(T.JOBT/10==0)
    CELL_EXAM(&C,&D,argc,argv);

  //======  launch evolutionary search module  =====
  if(T.JOBT/10==1)
    EVOS_MAIN(&T,&R,&P,&C);

  //========  launch cell simulation module  =======
  if(T.JOBT/10==2)
    CELL_MAIN(&R,&P,&C);

  //========  launch neural network module  ========
  if(T.JOBT/10>2) 
    NNET_MAIN(&R,&P,&C);

  return 0;
}
