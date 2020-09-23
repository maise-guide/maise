module maisedef

  use iso_c_binding

  implicit none

  type, bind(C) :: CELL

    real(c_double)                   ::  Ninv
    type(c_ptr)                      ::  R1
    type(c_ptr)                      ::  R2
    integer(c_int)                   ::  LI
    type(c_ptr)                      ::  VH
    integer(c_int)                   ::  NH
    real(c_double)                   ::  VM

    integer(c_int)                   ::  A
    integer(c_int)                   ::  ND
    integer(c_int)                   ::  N
    integer(c_int)                   ::  J
    integer(c_int)                   ::  M1

    real(c_double)                   ::  R0
    real(c_double)                   ::  dE

    type(c_ptr)                      ::  X
    type(c_ptr)                      ::  V
    type(c_ptr)                      ::  F
    type(c_ptr)                      ::  U
    type(c_ptr)                      ::  W
    real(c_double)                   ::  ANG(3)
    real(c_double)                   ::  LAT(3)
    real(c_double)                   ::  m(96)
    real(c_double)                   ::  mass(96)
    integer(c_int)                   ::  NRDF
    type(c_ptr)                      ::  RDF
    real(c_double)                   ::  p1
    real(c_double)                   ::  pint

    real(c_double)                   ::  P
    real(c_double)                   ::  K
    real(c_double)                   ::  H

    type(c_ptr)                      ::  L
    type(c_ptr)                      ::  R
    integer(c_int)                   ::  NM
    type(c_ptr)                      ::  Nn1
    type(c_ptr)                      ::  nn
    integer(c_int)                   ::  MNT

    integer(c_int)                   ::  nspc1
    integer(c_int)                   ::  spcz1(10)
    integer(c_int)                   ::  spcn1(10)
    integer(c_int)                   ::  NSPC
    integer(c_int)                   ::  SPCZ(10)
    integer(c_int)                   ::  SPCN(10)
    type(c_ptr)                      ::  ATMZ
    type(c_ptr)                      ::  ATMN

    type(c_ptr)                      ::  FF
    type(c_ptr)                      ::  BC
    real(c_double)                   ::  NE
    real(c_double)                   ::  E1
    type(c_ptr)                      ::  EA
    type(c_ptr)                      ::  min
    type(c_ptr)                      ::  max
    type(c_ptr)                      ::  Ni
    type(c_ptr)                      ::  S

    type(c_ptr)                      ::  e
    type(c_ptr)                      ::  MM
    integer(c_int)                   ::  NCH(3)
    real(c_double)                   ::  dCH(3)
    type(c_ptr)                      ::  CHG
    type(c_ptr)                      ::  NDX1
    type(c_ptr)                      ::  DX

    real(c_double)                   ::  Rmax
    real(c_double)                   ::  Rmin
    real(c_double)                   ::  DR
    real(c_double)                   ::  Rm(96)
      
    integer(c_int)                   ::  SGN
    character(kind=c_char)           ::  SGS(10)
    integer(c_int)                   ::  NSG
    type(c_ptr)                      ::  SG
    type(c_ptr)                      ::  SL
    integer(c_int)                   ::  NTE
    character(kind=c_char)           ::  ES(27)
    character(kind=c_char)           ::  PRS(3)

    integer(c_int)                   ::  NS

    character(kind=c_char)           ::  TAG
    character(kind=c_char)           ::  TYPE

    integer(c_int)                   ::  D
    type(c_ptr)                      ::  G
    type(c_ptr)                      ::  ndx
    type(c_ptr)                      ::  ndxj
    type(c_ptr)                      ::  cos
    type(c_ptr)                      ::  cosj
    type(c_ptr)                      ::  fc
    real(c_double)                   ::  Rc1
    real(c_double)                   ::  rc
    type(c_ptr)                      ::  Nj

    type(c_ptr)                      ::  FRC
    integer(c_int)                   ::  RLXT
    integer(c_int)                   ::  XT
    integer(c_int)                   ::  JOBT
    integer(c_int)                   ::  MODT
    integer(c_int)                   ::  it
    real(c_double)                   ::  LJa
    real(c_double)                   ::  LJe
    real(c_double)                   ::  DISP

    integer(c_int)                   ::  NW
    real(c_double), dimension(16,50) ::  WW
    character(kind=c_char)           ::  WDIR
    integer(c_int)                   ::  OUT
    integer(c_int)                   ::  POS
    integer(c_int)                   ::  NP
    integer(c_int)                   ::  NB
    character(kind=c_char)           ::  VER
    character(kind=c_char)           ::  ID
    integer(c_int)                   ::  MINT

    type(c_ptr)                      ::  EV1
    type(c_ptr)                      ::  ev
    integer(c_int)                   ::  EVOK

  end type CELL

  type, bind(C)  :: ANN

    character(kind=c_char)        ::  ID(200)
    integer(c_int)                ::  A
    integer(c_int)                ::  STR
    integer(c_int)                ::  N
    integer(c_int)                ::  TN
    integer(c_int)                ::  M
    integer(c_int)                ::  D1
    integer(c_int)                ::  O
    integer(c_int)                ::  NN
    integer(c_int)                ::  Q
    integer(c_int)                ::  NE
    integer(c_int)                ::  NF
    integer(c_int)                ::  TNF

    type(c_ptr)                   ::  E1
    type(c_ptr)                   ::  Fc
    type(c_ptr)                   ::  R1
    type(c_ptr)                   ::  C1
    type(c_ptr)                   ::  F
    type(c_ptr)                   ::  U
    type(c_ptr)                   ::  V
    type(c_ptr)                   ::  G
    type(c_ptr)                   ::  Rmin
    type(c_ptr)                   ::  Rmax
    type(c_ptr)                   ::  DR

    integer(c_int)                ::  JOBT
    integer(c_int)                ::  MODT
    integer(c_int)                ::  DSCR

    real(c_double)                ::  RE
    real(c_double)                ::  RF
    real(c_double)                ::  RS
    real(c_double)                ::  RT
    real(c_double)                ::  EE
    real(c_double)                ::  EF
    real(c_double)                ::  ES
    real(c_double)                ::  ET

    type(c_ptr)                   ::  train
    type(c_ptr)                   ::  test

    integer(c_int)                ::  DNm

    real(c_double)                ::  Eavg
    real(c_double)                ::  Edev

    type(c_ptr)                   ::  TEeval

    real(c_double)                ::  WE1
    real(c_double)                ::  WF
    real(c_double)                ::  WS

    real(c_double)                ::  time

    integer(c_int)                ::  EFS

    real(c_double)                ::  ECUT
    real(c_double)                ::  EMAX
    real(c_double)                ::  FMAX
    real(c_double)                ::  FMIN
    real(c_double)                ::  VMIN
    real(c_double)                ::  VMAX

    character(kind=c_char)        ::  depo(200)
    character(kind=c_char)        ::  data(200)
    character(kind=c_char)        ::  eval(200)
    character(kind=c_char)        ::  otpt(200)

    integer(c_int)                ::  NW
    integer(c_int)                ::  NL

    integer(c_int)                ::  NU(4)
    integer(c_int)                ::  GT(4)

    type(c_ptr)                   ::  B
    type(c_ptr)                   ::  W
    type(c_ptr)                   ::  e
    type(c_ptr)                   ::  d
    type(c_ptr)                   ::  c
    type(c_ptr)                   ::  Bp
    type(c_ptr)                   ::  Wp
    type(c_ptr)                   ::  WW
    type(c_ptr)                   ::  WWp
    type(c_ptr)                   ::  We

    integer(c_int)                ::  MINT
    integer(c_int)                ::  MITR

    real(c_double)                ::  ETOL
    real(c_double)                ::  LREG

    integer(c_int)                ::  seed
    integer(c_int)                ::  ITER

    character(kind=c_char)        ::  file_name(200)
    character(kind=c_char)        ::  test_add(200)

    integer(c_int)                ::  NSPC
    integer(c_int)                ::  SPCZ(10)

    integer(c_int)                ::  NSYM
    character(kind=c_char)        ::  compound(200)
    real(c_float)                 ::  seed2
    integer(c_int)                ::  N0
    integer(c_int)                ::  N1
    integer(c_int)                ::  NP
    integer(c_int)                ::  NB

    character(kind=c_char)        ::  file_list(11,200)

    type(c_ptr)                   ::  mask
    integer(c_int)                ::  nw1
    integer(c_int)                ::  MIX 
    integer(c_int)                ::  NSTP

    real(c_double)                ::  TMIN
    real(c_double)                ::  TMAX
    real(c_double)                ::  TSTP
    real(c_double)                ::  DELT
    real(c_double)                ::  CPLT
    real(c_double)                ::  CPLP
    real(c_double)                ::  ICMP

    integer(c_int)                ::  MDTP
    integer(c_int)                ::  MOVI
    character(kind=c_char)        ::  VER(200)
    real(c_double)                ::  Rc
    integer(c_int)                ::  PENE

  end type ANN

  type, bind(C) :: LNK

    integer(c_int)                ::  N
    type(c_ptr)                   ::  ATMN
    integer(c_int)                ::  NM
    integer(c_int)                ::  NF
    type(c_ptr)                   ::  DNn
    type(c_ptr)                   ::  DNi
    type(c_ptr)                   ::  DNs
    real(c_double)                ::  E 
    real(c_double)                ::  e1
    type(c_ptr)                   ::  F
    type(c_ptr)                   ::  f1
    type(c_ptr)                   ::  S
    type(c_ptr)                   ::  EA
    type(c_ptr)                   ::  s1
    type(c_ptr)                   ::  Cn
    type(c_ptr)                   ::  Fn
    type(c_ptr)                   ::  Sn
    real(c_double)                ::  p 
    type(c_ptr)                   ::  MRK
    type(c_ptr)                   ::  Fi
    integer(c_int)                ::  B
    character(kind=c_char)        ::  path(200)

  end type LNK

  type, bind(C) :: PRS

    integer(c_int)                :: B
    integer(c_int)                :: D
    type(c_ptr)                   :: Pl
    integer(c_int)                :: LN
    integer(c_int)                :: LM
    real(c_double)                :: RC
    type(c_ptr)                   :: GW
    integer(c_int)                :: GM
    integer(c_int)                :: GN
    type(c_ptr)                   :: GT
    type(c_ptr)                   :: GF
    type(c_ptr)                   :: GP
    real(c_double)                :: Rc1
    integer(c_int)                :: IO
    integer(c_int)                :: EFS
    integer(c_int)                :: DSCR
    real(c_double)                :: FMRK

    integer(c_int)                :: NG2
    integer(c_int)                :: NG4
    integer(c_int)                :: NSYM
    integer(c_int)                :: NSPC
    integer(c_int)                :: SPCZ(10)

  end type PRS

end module maisedef
