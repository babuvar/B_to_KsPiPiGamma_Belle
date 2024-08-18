c               
c    test_tatami_for.F   
c       ---- Test program using tatami with fortran
c    $Id: test_tatami_for.F 9948 2006-12-05 06:11:00Z katayama $ 
c                
c    $Log$
c    Revision 1.1  2002/09/05 01:27:02  katayama
c    New package tatami from Nakadaira/Sumisawa san
c
c    Revision 1.1  2002/05/07 21:38:53  nakadair
c    Add test programs.
c 
c                
      PROGRAM TESTCNVL
      Implicit none
      EXTERNAL IDEBUGOUT
      INTEGER IDEBUGOUT
      INTEGER EXPMC
      INTEGER SETRESP
      INTEGER STATUS

      DOUBLE PRECISION llmt,ulmt,granu,ak,ck,dm,tau,x,y1,y2,dy
      integer i,nmax
c     
#include "tatami/libcnvl_for.inc"
c     
      PARAMETER (EXPMC=1)
      PARAMETER (llmt=-20.0)
      PARAMETER (ulmt=20.0)
      PARAMETER (granu=0.1)
      PARAMETER (ak=0.98)
      PARAMETER (ck=0.5)
      PARAMETER (dm=0.4)
      PARAMETER (tau=1.5)
      
      STATUS = SETRESP(EXPMC);

      nmax = int((ulmt-llmt)/granu)
      
      open(99,FILE='epep.dat')
      
      do i=0,nmax-1
         x = llmt + i * granu;
         y1 = MfRk_fullrec(x, tau, dm, ak, ck)
         y2 = MfRk_partial(x, tau, dm, 0.002)
         dy = y2-y1
         write(99,*) x, y1, y2, dy
      enddo

      close(99)

      END
