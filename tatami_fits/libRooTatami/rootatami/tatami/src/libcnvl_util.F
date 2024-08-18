c               
c    libcnvl_util.f   
c       ---- Fortran utilities
c    $Id: libcnvl_util.F 8220 2002-09-05 01:27:07Z katayama $ 
c                
c    $Log$
c    Revision 1.1  2002/09/05 01:27:05  katayama
c    New package tatami from Nakadaira/Sumisawa san
c
c    Revision 1.2  2002/05/09 07:08:35  nakadair
c    Add add_outlier(...).
c
c    Revision 1.1  2002/05/07 21:39:21  nakadair
c    Add convolution with resolution function and fortran interface.
c 
c                
      
c
c     FUNCTION SETRESP
c     
c      fill /COMMON/ RESPER w/ default values
c 
c     argument EXPMC := 1 for real Data
c                    := 2 for MC data

      INTEGER FUNCTION SETRESP(EXPMC)
      Implicit none
      INTEGER EXPMC
c     
#include "tatami/libcnvl_for.inc"
c     
      SETRESP = cpdefpar(EXPMC)
#if 0
      IF ( EXPMC .eq. 1 ) THEN
         SREC(1)     = 0.809
         SREC(2)     = 0.154
         SASC(1)     = 0.753
         SASC(2)     = 0.064
         SMNREC      = 0.647
         STLREC      = 3.00
         FTLREC      = 0.083
         SMNASC      = 0.647
         STLASC      = 3.00
         FTLASC      = 0.083
         FDNPSG(1)   = 0.787
         FDNPSG(2)   = 0.763
         FPNPSG(1)   = 0.790
         FPNPSG(2)   = 0.757
         TNPPSG(1,1) = 0.108
         TNPPSG(2,1) = 1.321
         TNPPSG(1,2) = -0.0019
         TNPPSG(2,2) = 1.113
         TNPNSG(1,1) = -0.281
         TNPNSG(2,1) = 1.583
         TNPNSG(1,2) = -0.375
         TNPNSG(2,2) = 1.548
         FDNPML(1)   = 0.676
         FDNPML(2)   = 0.650
         FPNPML(1)   = 0.955
         FPNPML(2)   = 0.963
         TNPPML(1,1) = -0.010
         TNPPML(2,1) = 0.927
         TNPPML(1,2) = 0.037
         TNPPML(2,2) = 0.674
         TNPNML(1,1) = -0.194
         TNPNML(2,1) = 1.990
         TNPNML(1,2) = -0.269
         TNPNML(2,2) = 2.070
         TAUK1P(1)   = 0.15639
         TAUK1P(2)   = 0.017949
         TAUK2P(1)   = 0.048578
         TAUK2P(2)   = 0.17347E-17
         TAUK1N(1)   = 0.15639
         TAUK1N(2)   = 0.017949
         TAUK2N(1)   = 0.048578
         TAUK2N(2)   = 0.17347E-17
         SIGMAK(1)   = 0.000
         SIGMAK(2)   = 6.556701e-02
         FK2         = 0.2675694
         SIGOL       = 36.2
         FOLSG       = 0.0306
         FOLML       = 5.83e-4
         SETRESP     = 0
      ELSE IF (EXPMC .eq. 2 ) THEN
         SREC(1)     = 0.990
         SREC(2)     = 0.054
         SASC(1)     = 0.972
         SASC(2)     = 0.028
         SMNREC      = 1.052
         STLREC      = 4.436
         FTLREC      = 0.038
         SMNASC      = 1.052
         STLASC      = 4.436
         FTLASC      = 0.038
         FDNPSG(1)   = 0.787
         FDNPSG(2)   = 0.763
         FPNPSG(1)   = 0.790
         FPNPSG(2)   = 0.757
         TNPPSG(1,1) = 0.108
         TNPPSG(2,1) = 1.321
         TNPPSG(1,2) = -0.0019
         TNPPSG(2,2) = 1.113
         TNPNSG(1,1) = -0.281
         TNPNSG(2,1) = 1.583
         TNPNSG(1,2) = -0.375
         TNPNSG(2,2) = 1.548
         FDNPML(1)   = 0.676
         FDNPML(2)   = 0.650
         FPNPML(1)   = 0.955
         FPNPML(2)   = 0.963
         TNPPML(1,1) = -0.010
         TNPPML(2,1) = 0.927
         TNPPML(1,2) = 0.037
         TNPPML(2,2) = 0.674
         TNPNML(1,1) = -0.194
         TNPNML(2,1) = 1.990
         TNPNML(1,2) = -0.269
         TNPNML(2,2) = 2.070
         TAUK1P(1)   = 0.15639
         TAUK1P(2)   = 0.017949
         TAUK2P(1)   = 0.048578
         TAUK2P(2)   = 0.17347E-17
         TAUK1N(1)   = 0.15639
         TAUK1N(2)   = 0.017949
         TAUK2N(1)   = 0.048578
         TAUK2N(2)   = 0.17347E-17
         SIGMAK(1)   = 0.000
         SIGMAK(2)   = 6.556701e-02
         FK2         = 0.2675694
         SIGOL       = 43.63
         FOLSG       = 0.0022
         FOLML       = 9.9e-5
         SETRESP     = 0
      ELSE
         SETRESP = 1
      ENDIF
#endif
      RETURN
      END

