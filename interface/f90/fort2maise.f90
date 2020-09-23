program fort2maise
  
  use maise
  use iso_c_binding

  implicit none

  integer(c_int)                     :: i,j

  !===== pointers to MAISE structures ======
  type(ANN)                          :: R
  type(PRS)                          :: P
  type(PRS)                          :: W(9)
  type(LNK)                          :: L
  type(CELL)                         :: C

  !======== cell variables to pass =========
  integer(c_int)                     :: CODE  ! external code type
  integer(c_int)                     :: N     ! number of atoms
  integer(c_int)                     :: NM    ! max number of nearest neighbors
  integer(c_int)                     :: ND    ! clusters (0) or crystals (3)
  integer(c_int)                     :: NP    ! number of cores for parallelization
  integer(c_int)                     :: XT    ! fractional (0) or Cartesian (1)
  integer(c_int), dimension(4)       :: atmn  ! species types
  real(c_double), dimension(3, 3)    :: LAT   ! lattice vectors
  real(c_double), dimension(3, 4)    :: POS   ! positions

  !== returned enthalpy, forces, stresses ==
  real(c_double)                     :: H
  real(c_double), dimension(3, 4)    :: FRC
  real(c_double), dimension(6)       :: STR

  CODE   =  0;                                  ! currently 0 for all code types
  N      =  4; 
  NM     =  500; 
  ND     =  0; 
  NP     =  4; 
  XT     =  1;
  C % p  =  0.0;                                ! pressure in GPa for ND = 3
  
  LAT( 1, 1) = 20.d0; LAT( 2, 1) =  0.d0; LAT( 3, 1) =  0.d0;
  LAT( 1, 2) =  0.d0; LAT( 2, 2) = 20.d0; LAT( 3, 2) =  0.d0;
  LAT( 1, 3) =  0.d0; LAT( 2, 3) =  0.d0; LAT( 3, 3) = 20.d0;

  POS( 1, 1) =  2.d0; POS( 2, 1) =  0.0d0; POS( 3, 1) = 0.d0;
  POS( 1, 2) = -2.d0; POS( 2, 2) =  0.0d0; POS( 3, 2) = 0.d0;
  POS( 1, 3) =  0.d0; POS( 2, 3) =  2.0d0; POS( 3, 3) = 0.d0;
  POS( 1, 4) =  0.d0; POS( 2, 4) = -2.0d0; POS( 3, 4) = 0.d0;
  
  atmn(1) = 0; atmn(2) = 0; atmn(3) = 0; atmn(4) = 0;
  
  L % B = 0; 

  H = CALL_MAISE(R, P, W, L, C, CODE, N, NM, ND, NP, XT, atmn, LAT, POS, FRC, STR)
  
  print "(a)","Force components"
  do j = 1, 4
     print "(i3,xxf18.14,xxf18.14,xxf18.14)", j, FRC(1,j),FRC(2,j),FRC(3,j)
  end do

  print "(a,f18.14,/)","Enthalpy  ",H
  print "(a,/)","The forces and enthalpy should be consistent with the values in ./test/OUTCAR."
  
end program fort2maise
