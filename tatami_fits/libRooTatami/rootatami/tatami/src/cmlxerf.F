c
c     cmlxerf.F
c        --- Wrapper of complex exp,erf for C,C++
c
c     $Id: cmlxerf.F 8220 2002-09-05 01:27:07Z katayama $
c
c     $Log$
c     Revision 1.1  2002/09/05 01:27:03  katayama
c     New package tatami from Nakadaira/Sumisawa san
c
c     Revision 1.1.1.1  2002/02/21 18:56:07  nakadair
c     Initial import(Branch from lmfit)
c
c     Revision 1.1  2001/06/29 14:38:33  nakadair
c     Complex error function
c
c
      FUNCTION RECEXP(X, Y)
      IMPLICIT NONE
      DOUBLE PRECISION RECEXP, X, Y
      COMPLEX*16 Z, R
c     complex error function (CERNlib)
c
      Z = DCMPLX(X, Y)
      R = EXP(Z)
c      
      RECEXP = DBLE(R)
c
      RETURN
      END

      FUNCTION IMCEXP(X, Y)
      IMPLICIT NONE
      DOUBLE PRECISION IMCEXP, X, Y
      COMPLEX*16 Z, R
c     complex error function (CERNlib)
c
      Z = DCMPLX(X, Y)
      R = EXP(Z)
c      
      IMCEXP = DIMAG(R)
c
      RETURN
      END

      FUNCTION REWERF(X, Y)
      IMPLICIT NONE
      DOUBLE PRECISION REWERF, X, Y
      COMPLEX*16 Z, R, WWERF
c     complex error function (CERNlib)
      EXTERNAL WWERF
c
      Z = DCMPLX(X, Y)
      R = WWERF(Z)
c      
      REWERF = DBLE(R)
c
      RETURN
      END

      FUNCTION IMWERF(X, Y)
      IMPLICIT NONE
      DOUBLE PRECISION IMWERF, X, Y
      COMPLEX*16 Z, R, WWERF
c     complex error function (CERNlib)
      EXTERNAL WWERF
c
      Z = DCMPLX(X, Y)
      R = WWERF(Z)
c      
      IMWERF = DIMAG(R)
c
      RETURN
      END
