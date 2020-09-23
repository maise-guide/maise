module maise

 use maisedef
 use iso_c_binding

 implicit none

 interface
   function CALL_MAISE(R, P, W, L, C, CODE, N, NM, ND, NP, XT, ATMN, LAT, POS, FRC, STR) bind(C, name='CALL_MAISE')
     import c_int, c_double, c_float, c_char, c_ptr, ANN, PRS, LNK, CELL
     integer(c_int), intent(in), value :: CODE, N, NM, ND, NP, XT
     integer(c_int), intent(in)        :: ATMN(*) 
     real(c_double), intent(in)        :: LAT(3,3)
     type(CELL)                        :: C
     type(ANN)                         :: R
     type(PRS)                         :: P
     type(PRS)                         :: W(9)
     type(LNK)                         :: L
     real(c_double), intent(inout)     :: POS(3,*)
     real(c_double), intent(out)       :: FRC(3,*)
     real(c_double), intent(out)       :: STR(6)
     real(c_double)                    :: CALL_MAISE
   end function CALL_MAISE
 end interface

end module maise
