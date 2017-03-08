C SUBROUTINES WRITTEN BY PETER HOUTEKAMER FOR THE 3 -LEVEL T21
C MODEL VERSION 15-3-1991
      SUBROUTINE TIMINP (JTSTEP,XOLD,XNEW,DXNEW)
C VERSION OF TIMINT WITH ARGUMENTS
C
C *** JTSTEP = 0 : Forward integration, half time step
C *** JTSTEP = 1 : Leapfrog integration, half time step
C *** JTSTEP > 1 : Leapfrog integration + implicit diffusion
C ***              and time filter, full time step
C
c     IMPLICIT REAL*8 (A-H,O-Z)
c     include 'comfra.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      REAL*8 XOLD(NSH2,3),XNEW(NSH2,3),DXNEW(NSH2,3)

      real*8 as(nsh2,3)
cc    real*8 forcing,forcingpsi
cc    common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
cc    logical lforc,lforc2,lforc3,lforc4
cc    common /addforc/ lforc,lforc2,lforc3,lforc4
      INTEGER JTSTEP,I,J,K,L,M,N
      REAL*8 EPS1, EPS2

      EPS2=0.03
      EPS1=1.0D0-2.0D0*EPS2 
C
      IF (JTSTEP.GT.1) THEN
C
        DO 100 L=1,3
        DO 100 K=1,NSH2
        WS(K)=XNEW(K,L)
        XNEW(K,L)=DISS(K,1)*(XOLD(K,L)+DTT*DXNEW(K,L))
        XOLD(K,L)=EPS1*WS(K)+EPS2*(XNEW(K,L)+XOLD(K,L))
  100   CONTINUE
C
      ELSE
C
        IF (JTSTEP.EQ.0) THEN
          DO 200 L=1,3
          DO 200 K=1,NSH2
  200     XOLD(K,L)=XNEW(K,L)
        ENDIF
C
        DO 300 L=1,3
        DO 300 K=1,NSH2
        XNEW(K,L)=XOLD(K,L)+DTT*DXNEW(K,L)
  300   continue
C
c       DTT=2.*DTT
C
      ENDIF
C
      RETURN
      END


      SUBROUTINE LINDMU(DM,LFM,LLM,LM,EIEI)
C *** compute the linear operator for the mu derivate
C *** output
C     DM: packed array with derivative information
C     LFM: index to the first nonzero element in DM given m
C     LLM: index to the last nonzero element in DM given m
C     LM:  gives m if k is given
C     EIEI: the kinetic energy inproduct of the unit vectors
c     IMPLICIT REAL*8 (A-H,O-Z)
c     include 'comfra.f'
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
     
      INTEGER I,J,K,L,M,N 
      INTEGER KEPS(0:NM,0:NM),LFM(0:NM),
     &  LLM(0:NM),LM(NSH)
      REAL*8  ZERO(NSH2),EPS(NSH2),DEPSDM(NSH2),
     &  DM(NSH,0:NM),EIEI(NSH)
      
      DO 8,I=1,NM+1
 8      EIEI(I)=LL(I)*(LL(I)+1)
      DO 9,I=NM+2,NSH
 9      EIEI(I)=2.0*LL(I)*(LL(I)+1)
      M=-1
      K=0
      DO 1,I=NM+1,1,-1
        M=M+1
        LFM(M)=K+1
        DO 6,J=1,I
          K=K+1
          LM(K)=M
          N=LL(K)
          KEPS(M,N)=K
 6      CONTINUE
        LLM(M)=K
 1    CONTINUE
      DO 2,I=1,NSH2
 2      ZERO(I)=0.0D0
      DO 3,N=0,NM
        DO 7,I=1,NSH2
 7        EPS(I)=ZERO(I)
        DO 4,M=0,N
          K=KEPS(M,N)
          EPS(K)=1.0D0
 4      CONTINUE
        CALL DDMU(EPS,DEPSDM) 
C *** the imaginary part is trivial (identical to the real part)
        DO 5,K=1,NSH
 5        DM(K,N)=DEPSDM(K)
 3    CONTINUE
      RETURN
      END

      SUBROUTINE DDTL(DS,DEPSDT,jstep)
C *** COMPUTES THE TANGENT LINEAR EQUATIONS
c     IMPLICIT REAL*8 (A-H,O-Z)
c     include 'comfra.f'
c     include 'comtes.f'
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
c     LITEST VALUES WERE ALL 1, SO NO NEED FOR THIS
c     COMMON /CLINEA/ LITEST(10)
c     LOGICAL LITEST
      integer jstep,nnnstep
      REAL*8 CP(NSH2),WSPACE(NSH2),
     &  AS(NSH2,3),BS(NSH2,3),DS(NSH2,3),
     &  DEPSDT(NSH2,3),DUMMY,DFEPS(NSH2,3),DELPST(NSH2,2),
     &  DELPSI(NSH2,3),LPSI(NSH2,3)

      INTEGER I,J,K,L,M,N

      real*8 xs(nsh2,3),fac
cc    real*8 forcing,forcingpsi
cc    common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
cc    logical lforc,lforc2,lforc3,lforc4
cc    common /addforc/ lforc,lforc2,lforc3,lforc4
      integer nstop,janstep
      common /ntime/ nstop
      
C *** Initializations needed for the computation of c1,c2 and c3
      DO 90 L=1,3
       DO 91 K=1,NSH2
        DFEPS(K,L)=0.0D0
   91  CONTINUE
   90 CONTINUE
      DO 2 K=1,NSH2
        CP(K)=0.0D0
    2 CONTINUE
      CP(2)=1.0D0/SQRT(3.0D0)
      
      CALL COPY(DS,DELPSI)
 
C     MATRIX L TOEPASSEN:
      CALL PSIQ(DELPSI,AS)
      CALL PSIQ(PSI,LPSI)
      DO 200 K=1,NSH2
        DELPST(K,1)=RELT1*(DELPSI(K,1)-DELPSI(K,2))
        DELPST(K,2)=RELT2*(DELPSI(K,2)-DELPSI(K,3))
 200  CONTINUE
      
c     IF (LITEST(1)) THEN
        CALL SUB (AS(1,1),CP,BS(1,1))  
        CALL JACOB(PSI(1,1),BS(1,1),WSPACE(1))
        CALL ADDAB (DFEPS(1,1),WSPACE(1))
c      ENDIF
      
c     IF (LITEST(2)) THEN
        CALL JACOB(DELPSI(1,1),LPSI(1,1),BS(1,1))
        CALL ADDAB (DFEPS(1,1),BS(1,1))
c     ENDIF
      
c     IF (LITEST(3)) THEN
        CALL ADDAB (DFEPS(1,1),DELPST(1,1))
c     ENDIF 
      
c     IF (LITEST(4)) THEN
        CALL SUB (AS(1,2),CP,BS(1,2))  
        CALL JACOB(PSI(1,2),BS(1,2),WSPACE(1))
        CALL ADDAB (DFEPS(1,2),WSPACE(1))
c     ENDIF
      
c     IF (LITEST(5)) THEN
        CALL JACOB(DELPSI(1,2),LPSI(1,2),BS(1,2))
        CALL ADDAB (DFEPS(1,2),BS(1,2))
c     ENDIF

c     IF (LITEST(6)) THEN
        CALL SUBAB (DFEPS(1,2),DELPST(1,1))
        CALL ADDAB (DFEPS(1,2),DELPST(1,2))
c     ENDIF
      
c     IF (LITEST(7)) THEN
        CALL SUB (AS(1,3),CP,BS(1,3))  
        CALL JACOB(PSI(1,3),BS(1,3),WSPACE(1))
        CALL ADDAB (DFEPS(1,3),WSPACE(1))
c     ENDIF

c     IF (LITEST(8).OR.LITEST(10)) THEN
        CALL ADD(LPSI(1,3),OROG,WSPACE(1))
        CALL JACOBP(DELPSI(1,3),WSPACE(1),BS(1,3))
        CALL ADDAB (DFEPS(1,3),BS(1,3))
c     ENDIF
      
c     IF (LITEST(9)) THEN
        CALL SUBAB (DFEPS(1,3),DELPST(1,2))
c     ENDIF

      CALL QPSI(DFEPS,DEPSDT)

c pick appropriate forcing

cc    if (lmultiforce) then
cc    if (0.eq.1) then
cc      do k=1,numforce
	  
cc        if (jstep.eq.(nforceswitch(k))) then  
cc        if (0.eq.1) then 
cc	   write(*,*) 'activate forcing',k,'  at step',nforceswitch(k)
cc	    do i=1,3
cc            do j=1,nsh2
cc              forcingpsi(j,i)=multiforcing(j,i,k)
cc            enddo
cc          enddo
cc	  endif
cc  	enddo
cc    endif

55    FORMAT(A7,5F12.8) 
cc    if (lforc4) then
cc    if (0.eq.1) then
cc      do i=1,3
cc        do j=1,nsh2
cc          DEPSDT(j,i)=DEPSDT(j,i)+forcingpsi(j,i)
cc        enddo
cc      enddo
cc    endif
      
      RETURN
      END

      SUBROUTINE DDMU(EPS,DEPSDM)
C *** eps is a field in spherical harmonics
C *** depsdm is the mu-derivative in spherical harmonics
c     IMPLICIT REAL*8 (A-H,O-Z)
c     include 'comfra.f'
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      REAL*8 EPS(NSH2),DEPSDM(NSH2),DEPSGA(NGP)

C *** compute the derivate on the GAUSS-grid      
      CALL SPTOGG(EPS,DEPSGA,PD)
C *** backtransformation to spherical harmonics
      CALL GGTOSP(DEPSGA,DEPSDM)

      RETURN
      END


      SUBROUTINE TANG(AS,BS,NUUR)
C**   this subroutine integrates the tangent
C**   linear equation with initial condition 
C**   AS(NSH2,3) and integration time NUUR. 
C**   It uses the orbit stored by STOREP.
C**   The result is BS(NSH2,3)
C
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      INTEGER NLMAX
C    *  NPERT,NSKIP,NOUTST,IDIF,NPRINT
      PARAMETER (nlmax=800)
      REAL*8 RMSV
C    *  DT,EPS,TDIS,ADDISL,ADDISH,TREL,TDIF,H0,
c    *  RRDEF1,RRDEF2
c     PARAMETER ( NPERT=0, NSKIP=0, RMSV=1.E-02 )
c     PARAMETER ( NOUTST=36, DT=1./36., EPS=.03,
c    *            TDIS=3.0, ADDISL=0.5, ADDISH=0.5,
c    *            TREL=25., TDIF=2.0, IDIF=4,
c    *            H0=9., RRDEF1=.110, RRDEF2=.070 )
C
C      INCLUDE COMFRA.F
C      INCLUDE COMTES.F
c     include 'comfra.f'
c     include 'comtes.f'

      REAL*8 X(NLAT),W(NLAT),RAD(NLAT),PIGR4,TREAL
      REAL*8 DUM(NSH23),DUMOLD(NSH23),DUMNEW(NSH23)
      REAL*8 AS(NSH2,3),BS(NSH2,3),cs(nsh2,3),e,fac
cc    real*8 forcing,forcingpsi
cc    common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
cc    logical lforc,lforc2,lforc3,lforc4
cc    common /addforc/ lforc,lforc2,lforc3,lforc4
      integer janstep,k,index
      common /stap/ janstep
      integer nstop
      common /ntime/ nstop
      
      real*8 zt(nlat,nlon),dummy,sum
      integer ix,iy,ilevel
      INTEGER LM(NSH)
C
      INTEGER INQ,INP,INH,NTSTEP,NLANCO,JTSTEP,JRUN,
     *  ISEED,I,J,ILANCO,n
      INTEGER NUUR
      NTSTEP=(NUUR*3)/2
      nstop=NTSTEP
      
cc    if(lforc) then
cc    write(*,*) 'ntstep=',janstep
cc    ntstep=janstep
cc    end if
      
      CALL COPY(AS,DUMNEW)
C
      DO 106,I=1,NSH23
      DUMOLD(I)=0.0D0
  106 CONTINUE
      WRITE(*,*) 'TANGENT MODEL'
c      PIGR4=8.0D0*ASIN(1.0D0)
c      DTT=DTT/2.0D0
       print *,"IN TANGENT, DTT=",DTT,NTSTEP
       DO 107,JTSTEP=0,NTSTEP
c        CALL RETRIE(JTSTEP)
         CALL RETRIE(0)
         CALL DDTL(DUMNEW,DUM,jtstep)
         CALL TIMINP(JTSTEP,DUMOLD,DUMNEW,DUM)
	 
	 
c	 if(ltangwr) then
c	   k=jtstep
c           if((mod(k,nsamplestep).eq.0)
c     &       .and.(k.ge.((ntrajstart*3)/2))) then
c             call inprod(dumnew,dumnew,e)
c	     write(45,*) 'TANG',k,e
c	   endif
c	 endif
	 
cc	 if(ltangwr) then
         if(0.eq.1) then
	   k=jtstep
           if(mod(k,9).eq.0) then
	     write(*,*) 'write-out linear evo',k
	     if(k.eq.0) then
	       write(66) as
	     else
	       write(66) dumnew
	     endif
	   endif
	 endif
	 
         if(0.eq.1) then
	   k=jtstep
	   if(k.eq.0) then
	     write(*,*) 'WRITE LIN TRAJ to 51 !!!'
	   endif
           if(mod(k,9).eq.0) then
	     k=(jtstep/3)*2
c             call tometview(51,dumnew,k,138)
              write(66) dumnew
	   endif
	 endif
	 
cc       if(ltrajcost) then
         if(0.eq.1) then	 
c fill array ETRAJ to determine costfunction

	   k=jtstep
cc         if((mod(k,nsamplestep).eq.0)
cc   &       .and.(k.ge.((ntrajstart*3)/2)).and.(k.lt.ntstep)) then
             if(0.eq.1) then
cc	     index=( k-((ntrajstart*3)/2) )/nsamplestep
  	     call copy(dumnew,etraj(1,1,index))
	     
c one gridpoint (vorticity)  --------------------------------
             if (0.eq.1) then
               call lap(etraj(1,3,index),etraj(1,3,index))
               call sptogg(etraj(1,3,index),zt,pp)
	       dummy=zt(24,30)
c	     if (dummy.lt.0) then
c               write(*,*) 'NEGATE'
c	       dummy=-dummy
c	     endif
	       write(*,*) 'etraj(index)',dummy
	       if(jtstep.eq.0) sum=0
	       sum=sum+dummy**2
	       if(jtstep.eq.ntstep) write(*,*) 'costf = ',sum
	       do iy=1,nlat
	         do ix=1,nlon
	           zt(iy,ix)=0.0
	         enddo
	       enddo
	       zt(24,30)=dummy
	     
  	       call ggtosp(zt,etraj(1,3,index))
  	       call lap(etraj(1,3,index),etraj(1,3,index))
	       do ilevel=1,2
	         do n=1,nsh2
  	           etraj(n,ilevel,index)=0.0
	         enddo
	       enddo
	     endif
	  endif
	endif
c	if(jtstep.eq.ntstep) write(*,*) 'costf = ',sum   
 107  CONTINUE
      CALL COPY(DUMNEW,BS)
C
      RETURN
      END


      SUBROUTINE NONL(NUUR)
c----------------------------------------------------------------------
C**   this subroutine computes the T21 orbit with
C**   initial condition QPRIME(NSH2,3). This field can
C**   be obtained by using for example JINPUT.
C**   The forecast period is NUUR.
c---------------------------------------------------------------------
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      INTEGER NLMAX
C    *  NPERT,NSKIP,NOUTST,IDIF,NPRINT
      PARAMETER (nlmax=800)
      REAL*8 RMSV
C    *  DT,EPS,TDIS,ADDISL,ADDISH,TREL,TDIF,H0,
c    *  RRDEF1,RRDEF2
c     PARAMETER ( NPERT=0, NSKIP=0, RMSV=1.E-02 )
c     PARAMETER ( NOUTST=36, DT=1./36., EPS=.03,
c    *            TDIS=3.0, ADDISL=0.5, ADDISH=0.5,
c    *            TREL=25., TDIF=2.0, IDIF=4,
c    *            H0=9., RRDEF1=.110, RRDEF2=.070 )
C
C      INCLUDE COMFRA.F
C      INCLUDE COMTES.F
c     include 'comfra.f'
c     include 'comtes.f'

      REAL*8 X(NLAT),W(NLAT),RAD(NLAT),PIGR4,TREAL
      INTEGER LM(NSH)
C
      INTEGER INQ,INP,INH,NTSTEP,NLANCO,JTSTEP,JRUN,
     *  I,J
      INTEGER NUUR
      logical lforc,lforc2,lforc3,lforc4,lforcnl
      common /addforc/ lforc,lforc2,lforc3,lforc4
      integer janstep
      common /stap/ janstep
      common /addforc2/ lforcnl

      INTEGER K,L,M,N
      REAL*8 QPRIM1(NSH2,NVL)

      NTSTEP=(NUUR*3)/2
cjan
      if(lforc) then
      ntstep=janstep
      end if

      WRITE(*,*) 'NUMBER OF STEPS IN NONL IS: ',NTSTEP 
C INITIALIZE THE MODEL EQUATIONS
      DO 4,I=1,3
        DO 5,J=1,NSH2
          QPRIM1(J,I)=0.0D0
 5      CONTINUE
 4    CONTINUE
      WRITE(*,*) 'START THE INTEGRATION'
c      IF (NTSTEP.GT.NTADJ) THEN
c        WRITE(*,*) 'NUMBER OF STEPS TO LARGE'
c        PAUSE 'PRESS RETURN TO ABORT'
c        STOP
c      ENDIF
      WRITE(*,*) 'THE NON-LINEAR MODEL'
c     PIGR4=8.0D0*ASIN(1.0D0)
c     DTT=DT*PIGR4/2.0D0
c     write(*,*) 'DTT =',DTT
      DO 10 JTSTEP=0,NTSTEP
      if (jtstep.eq.99999) then
        write(*,*) 'switch off forcing'
	lforcnl=.false.
      endif
C *** QPRIME IS KNOWN GENERATE THE REST
        CALL QTOPSI
        CALL STOREP(JTSTEP)
        CALL DDT
ccccc        CALL DDTP
        CALL TIMINP(JTSTEP,QPRIM1,QPRIME,DQPRDT)
   10 CONTINUE
C
      CALL QTOPSI
C
      WRITE(*,*) 'END OF THE NON-LINEAR INTEGRATION'
C
      RETURN
      END


      SUBROUTINE JNONL(AS,NUUR)
C* this routine integrates the nonlinear model
C* with initial condition AS and forecast time NUUR
C* AS is assumed to be in streamfunction coordinates.
C
c     include 'comfra.f'
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      REAL*8 AS(NSH2,3)
      INTEGER NUUR
      INTEGER K,L,M,N
C
      CALL COPY(AS,PSI)
C
      DO 110 L=1,3
      CALL LAP(PSI(1,L),QPRIME(1,L))
 110  CONTINUE
      DO 120 K=1,NSH2
      PSIT(K,1)=PSI(K,1)-PSI(K,2)
      PSIT(K,2)=PSI(K,2)-PSI(K,3)
      QPRIME(K,1)=QPRIME(K,1)-RL1*PSIT(K,1)
      QPRIME(K,2)=QPRIME(K,2)+RL1*PSIT(K,1)-RL2*PSIT(K,2)
      QPRIME(K,3)=QPRIME(K,3)+RL2*PSIT(K,2)    
  120 CONTINUE     
C
      CALL NONL(NUUR)
       
      RETURN
      END
