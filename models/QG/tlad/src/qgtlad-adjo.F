      SUBROUTINE DDTLAD(DS,DSFM,jstep)
C *** Computes the adjoint of the tangent linear equations
C *** stripped down efficient version
C *** assumes a correct version of PSI 
cc    IMPLICIT REAL*8 (A-H,O-Z)
cc    include 'comfra.f'
cc    include 'comtes.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
 
cc    COMMON /CLINEA/ LITEST(10)
cc    LOGICAL LITEST
cc    logical lforc,lforc2,lforc3,lforc4
cc    common /addforc/ lforc,lforc2,lforc3,lforc4
cc    logical lforcwr
cc    common /wrforc/ lforcwr
      integer jstep
      integer nstop
      common /ntime/ nstop
cc    real*8 forcing,forcingpsi
cc    common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)

      REAL*8 CP(NSH2),WSPACE(NSH2),WS1(NSH2),DSMT(NSH2,2),
     &  AS(NSH2,3),CS(NSH2,3),DS(NSH2,3),DSM(NSH2,3),
     &  DSMLAP(NSH2,3),DSFM(NSH2,3),LPSI(NSH2,3),SJACOB(NSH2),
     &  DSMLGA(NLAT,NLON)
      REAL*8 PIGR4,e
      
      integer nstep
      INTEGER I,J,K,L,M,N
      
cc    if (lforc3) then
      
cc      if (jstep.gt.ntadj) then
cc       write(*,*) 'too many forcing steps'
cc       stop
cc      endif
	
cc      call copy(forcingtraj(1,1,jstep-1),forcingpsi)

cc      PIGR4=8.0D0*ASIN(1.0D0)
cc	DT=1./36.
cc      fac=1.0
	
c have a look in the QG paper:
c the identity operator is only active during the application
c time of the forcing (in which case fac=0)
	
cc      do i=1,3
cc        do j=1,nsh2
cc          do kforce=1,numforce
cc            if (kforce.eq.nforce) then
cc              if ( jstep.eq.(nforceswitch(kforce+1)+1)) then
cc		  if(i*j.eq.1) write(*,*) 'nullify',jstep
cc                dadnew(j,i)=0.0
cc	          dadold(j,i)=0.0
cc              endif
cc		if (  ( jstep.le.(nforceswitch(kforce+1))).and.
cc   &              (jstep.gt.(nforceswitch(kforce))) ) then
cc                fac=0.0
cc		endif
cc		if ((kforce.eq.1).and.(jstep.eq.1)) fac=0.0  
cc		forcingpsi(j,i)=(1.0-fac)*forcingpsi(j,i)
cc		
cc            endif
cc          enddo
cc
cc	    DSFM(j,i)=forcingpsi(j,i)
cc	  enddo 
cc      enddo
cc	
cc    else

cc      if(lforcwr) then
cc        write(*,*) 'write',jstep-1
cc        call copy(ds,forcingtraj(1,1,jstep-1))
cc      endif

      DO 90 L=1,3
      DO 91 K=1,NSH2
        DSFM(K,L)=0.0D0
   91 CONTINUE
   90 CONTINUE
      DO 2 K=1,NSH2
        CP(K)=0.0D0
    2 CONTINUE
      CP(2)=1.0D0/SQRT(3.0D0)
      
      CALL PSIQ(PSI,LPSI)
C     MATRIX M TOEPASSEN:
      CALL QPSI(DS,DSM)
      DO 3 L=1,3
 3      CALL LAP(DSM(1,L),DSMLAP(1,L))
      DO 6,K=1,NSH2
        DSMT(K,1)=RELT1*(DSM(K,1)-DSM(K,2))
        DSMT(K,2)=RELT2*(DSM(K,2)-DSM(K,3))
 6    CONTINUE

C *** COMPUTE J(L1(eps),psi1) + J(L2(eps),psi2) + J(L3(eps,psi3))

C       M toepassen
        DO 30,L=1,3
         CALL SUB(PSI(1,L),CP,CS(1,L))
         CALL JACOB(DSMLAP(1,L),CS(1,L),WSPACE(1))
         CALL LAPINV(WSPACE(1),AS(1,L))
 30     CONTINUE
C     L toepassen
        CALL PSIQ(AS,CS)
        DO 4 L=1,3
 4        CALL ADDAB(DSFM(1,L),CS(1,L))

C *** Compute J(L1(psi)+f,eps1) + J(L2(psi)+f,eps2)
        DO 32,L=1,2
         CALL JACOB(DSMLAP(1,L),LPSI(1,L),CS(1,L))
         CALL LAPINV(CS(1,L),AS(1,L))
         CALL SUBAB(DSFM(1,L),AS(1,L))
 32     CONTINUE
      
C *** Compute 1/TuaR (Del2 eps, Del4 eps -Del2 eps, -Del4 eps
        CALL ADDAB(DSFM(1,1),DSMT(1,1))
        CALL SUBAB(DSFM(1,2),DSMT(1,1))
        CALL ADDAB(DSFM(1,2),DSMT(1,2))
        CALL SUBAB(DSFM(1,3),DSMT(1,2))

C *** compute J(l3(psi)+f(1+h/H0),eps3)
        CALL ADD(LPSI(1,3),OROG,WSPACE(1))
        CALL JACOB(DSMLAP(1,3),WSPACE(1),WS1(1))
        CALL LAPINV(WS1(1),AS(1,3))
        CALL SUBAB(DSFM(1,3),AS(1,3))

C *** 1/tau del eps

C    Adjoint of the isotropic dissipation term
        CALL SPTOGG(DSMLAP(1,3),DSMLGA,PP)
        DO 253,I=1,NLAT
          DO 253, J=1,NLON
 253      ININAG(I,J)=RDISS(I,J)*DSMLGA(I,J)
        CALL GGTOSP(ININAG,SJACOB)
        DO 252,K=1,NSH2
 252      WSPACE(K)=DISS(K,2)*SJACOB(K)
        CALL LAPINV(WSPACE,SJACOB)
        CALL ADDAB(DSFM(1,3),SJACOB)
C    Adjoint of the L-derivative
        DO 254,I=1,NLAT
          DO 254 J=1,NLON
 254      ININAG(I,J)=DDISDX(I,J)*DSMLGA(I,J)
        CALL GGTOSP(ININAG,WSPACE)
        CALL DDL(WSPACE,SJACOB)
        CALL LAPINV(SJACOB,WSPACE)
        CALL ADDAB(DSFM(1,3),WSPACE)
C    Adjoint of the M-derivative
        CALL MUDERA(DSM(1,3),SJACOB)
        CALL ADDAB(DSFM(1,3),SJACOB)


cc    endif

      RETURN
      END

      SUBROUTINE DDTAD(JTSTEP)
C THE ADJOINT OF DDT
cc    IMPLICIT REAL*8 (A-H,O-Z)
cc    include 'comfra.f'
cc    include 'comtes.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
      REAL*8 QPRIM(NSH23),DQPRIM(NSH23),as(nsh2,3)
cc    logical lforc,lforc2,lforc3,lforc4
cc    common /addforc/ lforc,lforc2,lforc3,lforc4
      integer nstop
      common /ntime/ nstop
      INTEGER I,J,K,L,M,N
      INTEGER ILONG,JTSTEP
      REAL*8 EPS1,EPS2,PIGR4 
      EPS2=0.04
      EPS1=1.0D0-2.0D0*EPS2
c     PIGR4=8.0D0*ASIN(1.0D0)

      DO 1,I=1,3
        DO 2,J=1,NSH2
          ILONG=J+(I-1)*NSH2
          IF (JTSTEP.GT.2) THEN
           QPRIM(ILONG)=DADNEW(J,I)*DISS(J,1)+
     &      DADOLD(J,I)*DISS(J,1)*EPS2
          ELSE
           QPRIM(ILONG)=DADNEW(J,I)
          ENDIF
 2      CONTINUE
 1    CONTINUE

c       call length(dadold,e)
c       write(*,*) 'dadold =', e
cj        write(*,*) 'ddtad, jtstep =', jtstep,nstop
c      CALL DDTLAD(QPRIM,DQPRIM)
      CALL DDTLAD(QPRIM,DQPRIM,jtstep)
      DO 5,I=1,3
        DO 6,J=1,NSH2
          ILONG=J+(I-1)*NSH2
          DDADNE(J,I)=DQPRIM(ILONG)
 6      CONTINUE
 5    CONTINUE

      
      RETURN
      END

      SUBROUTINE TIMIAD (JTSTEP)
C THE ADJOINT OF TIMINT
C AUTHOR PETER HOUTEKAMER
C DATE 21-1-91
C *** JTSTEP = 0 : Forward integration, half time step
C *** JTSTEP = 1 : Leapfrog integration, half time step
C *** JTSTEP > 1 : Leapfrog integration + implicit diffusion
C ***              and time filter, full time step
C
cc    IMPLICIT REAL*8 (A-H,O-Z)
cc    include 'comfra.f'
cc    include 'comtes.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
c     real*8 forcing,forcingpsi
c     common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
      real*8 amp
      common amp
c     logical lforc,lforc2,lforc3,lforc4
c     common /addforc/ lforc,lforc2,lforc3,lforc4
      real*8 as(nsh2,3),zz(nsh2,3),uu(nsh2,3)
      common /ntime/ nstop
      integer nstop
      INTEGER JTSTEP,I,J,K,L,M,N
      REAL*8 EPS1,EPS2,PIGR4
      EPS2=0.03
      EPS1=1.0D0-2.0D0*EPS2
c     PIGR4=8.0D0*ASIN(1.0D0)

      IF (JTSTEP.GT.2) THEN
C
        print *,"I AM NOT HERE"
        DO 100 L=1,3
        DO 100 K=1,NSH2
        WS(K)=DADOLD(K,L)
        DADOLD(K,L)=EPS2*(1.0+DISS(K,1))*WS(K)+
     &    DISS(K,1)*DADNEW(K,L)
        DADNEW(K,L)=EPS1*WS(K)+DTT*DDADNE(K,L)
 100   CONTINUE
        
      ELSE
C

C to test with ntstep=0
cjan        DTT=0.5*DTT

cc      PIGR4=8.0D0*ASIN(1.0D0)
cc      DT=1./36.
c**        DTT=DT*PIGR4*2.0D0
cc      IF (JTSTEP.EQ.0) DTT=DT*PIGR4*2.0D0/2. 
cc      IF (JTSTEP.EQ.1) DTT=DT*PIGR4*2.0D0/4.
cc      IF (JTSTEP.EQ.2) DTT=DT*PIGR4*2.0D0/2.
C
        IF (JTSTEP.EQ.1) THEN
          PRINT *,"I AM HERE",DTT
c         DTT=DTT*2.
          DO 200 L=1,3
           DO 200 K=1,NSH2
           DADNEW(K,L)=DADOLD(K,L)+DADNEW(K,L)+DTT*DDADNE(K,L)
  200     continue

        ELSE IF (JTSTEP.EQ.2) THEN
c       DTT=2.*DTT
        PRINT *,"I AM PERHAPS HERE",DTT 
C
         DO 300 L=1,3
          DO 300 K=1,NSH2
           DADOLD(K,L)=DADNEW(K,L)+DADOLD(K,L)
           DADNEW(K,L)=DTT*DDADNE(K,L)
  300    CONTINUE

        ENDIF
C
      ENDIF
C
c      write(*,*) 'step in adjoint model =',jtstep,DTT
      RETURN
      END


      SUBROUTINE ADJSQ(AS,BS,NUUR)
C this subroutine integrates the
C adjoint equations with initial condition
C AS and with integration time NUUR
C the result is BS.
C the adjoint here is determined with respect to the 
C squared norm inner product
cc    IMPLICIT NONE
cc    include 'comfra.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
      INTEGER NUUR,L
      REAL*8 AS(NSH2,3),BS(NSH2,3)
      REAL*8 ASC(NSH2,3),BSC(NSH2,3)
 
      DO 1,L=1,3
        CALL LAPINV(AS(1,L),ASC(1,L))
 1    CONTINUE
      CALL ADJ(ASC,BSC,NUUR)
      DO 2,L=1,3
	CALL LAP(BSC(1,L),BS(1,L))
 2    CONTINUE
      RETURN
      END


      SUBROUTINE MUDERA (XS,YS)
C***  computes the 'adjoint' of the 'mu-derivative'
C***  input is XS, output is YS.
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
      REAL*8 XS(NSH2),YS(NSH2),ZS(NSH2)
      INTEGER I,J,K,L,M,N
      REAL*8 DERIMU(NMAT,NMAT)
 
      DO 1 K=1,NSH2
      YS(K)=0.
      ZS(K)=0.
    1 CONTINUE
    
      DO 10 K=2,22
      XS(K)=XS(K)*LL(K)*(LL(K)+1)
   10 CONTINUE
      DO 20 K=23,NSH
      XS(K)=XS(K)*LL(K)*(LL(K)+1)*2
   20 CONTINUE
      DO 30 K=NSH+23,NSH2
      XS(K)=XS(K)*LL(K-NSH)*(LL(K-NSH)+1)*2
   30 CONTINUE
      
      DO 2 K=NSH+1,NMAT 
      XS(K)=XS(K+22)
    2 CONTINUE
    
C**   transpose of DERIMU is used.
    
      DO 3 I=2,NMAT
      DO 4 K=2,NMAT
      ZS(I)=DERIMU(K,I)*XS(K)+ZS(I)
    4 CONTINUE
    3 CONTINUE
    
      DO 5 K=2,22
      YS(K)=ZS(K)/(LL(K)*(LL(K)+1))
    5 CONTINUE
      DO 6 K=23,NSH
      YS(K)=ZS(K)/(LL(K)*(LL(K)+1)*2)
    6 CONTINUE
      DO 7 K=NSH+23,NSH2
      YS(K)=ZS(K-22)/(LL(K-NSH)*(LL(K-NSH)+1)*2)
    7 CONTINUE
      
      RETURN
      END
      

      SUBROUTINE ADJ(AS,BS,NUUR) 
C**   this subroutine integrates the adjoint
C**   equation with initial condition AS(NSH2,3)
C**   and with integration time NUUR. It uses the orbit
C**   stored by STOREP.
C**   The result is BS(NSH2,3). The adjoint is determined with
C**   respect to the kinetic energy inner product.
C BY PETER HOUTEKAMER (PH) AND JAN BARKMEIJER (JB)
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      INTEGER NLMAX,NPERT,NSKIP,NOUTST,NPRINT
      PARAMETER (nlmax=800)
cc    REAL*8 RMSV,DT,EPS,TDIS,ADDISL,ADDISH,TREL,TDIF,H0,
cc   *  RRDEF1,RRDEF2
cc    PARAMETER ( NPERT=0, NSKIP=0, RMSV=1.E-02 )
cc    PARAMETER ( NOUTST=36, DT=1./36., EPS=.03,
cc   *            TDIS=3.0, ADDISL=0.5, ADDISH=0.5,
cc   *            TREL=25., TDIF=2.0, IDIF=4,
cc   *            H0=9., RRDEF1=.110, RRDEF2=.070 )
C
C      INCLUDE COMFRA.F
C      INCLUDE COMTES.F
cc    include 'comfra.f'
cc    include 'comtes.f'

      REAL*8 X(NLAT),W(NLAT),RAD(NLAT),PIGR4,TREAL
      REAL*8 AS(NSH2,3),BS(NSH2,3)
cc    real*8 forcing,forcingpsi
cc    common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
cc    real*8 forcingtraj
cc    common forcingtraj(nsh2,3,0:120)
cc    logical lforc,lforc2,lforc3,lforc4
cc    common /addforc/ lforc,lforc2,lforc3,lforc4
cc    logical lforcwr
cc    common /wrforc/ lforcwr
      integer nstop,janstep,kforce
      common /ntime/ nstop
      common /stap/ janstep
      real*8 uold(nsh2,3),unew(nsh2,3),e,e1,e2,fac,
     &       zold(nsh2,3),znew(nsh2,3),cs(nsh2,3)
      
      INTEGER LM(NSH)
C
      INTEGER INQ,INP,INH,NTSTEP,NLANCO,JTSTEP,JRUN,
     *  ISEED,I,J,ILANCO,m,k,nstep,index
      INTEGER NUUR
      DATA INQ/ 2/, INP/ 12/, INH/ 3/, INF/ 4/
C
      NTSTEP=(NUUR*3)/2
      nstop=ntstep
      
cc    if(lforc) then
cc    write(*,*) 'ntstep',janstep
cc    ntstep=janstep
cc    nstop=ntstep
cc    end if

cc    if(lforcwr) then
cc      do i=1,nsh2
cc        do j=1,3
cc          do m=0,120
cc            forcingtraj(i,j,m)=0.0
cc          enddo
cc        enddo
cc      enddo 
cc    endif
    
      
C
C *** verify the adjoint of the integration scheme
C INITIALIZE THE MODEL EQUATIONS
C
      DO 10 I=1,NSH2
      DO 20 J=1,3
      DADOLD(I,J)=0.0D0
      FORNEW(I,J)=0.0D0
      FOROLD(I,J)=0.0D0
   20 CONTINUE
   10 CONTINUE
      WRITE(*,*) 'ADJOINT MODEL'
c     PIGR4=8.0D0*ASIN(1.0D0)
      DTT=DTT/2.0D0
      print *,"IN ADJOINT, DTT = ",DTT
      CALL COPY(AS,DADNEW)

cc    PIGR4=8.0D0*ASIN(1.0D0)
cc    DTT=DT*PIGR4*2.0D0

C DADNEW is hier de foutenvektor op tijdstip T
      WRITE(*,*) 'START THE BACKWARD INTEGRATION'
C het terugintegreren met de adjoint
C iets anders bedenken voor ntstep <=1

      DO 200,JTSTEP=NTSTEP+1,1,-1	
	  
C *** QPRIME(=PSI-ORBIT) IS STORED
cc      CALL RETRIE(JTSTEP-1)
        CALL RETRIE(0)
        CALL DDTAD(JTSTEP)
        CALL TIMIAD(JTSTEP)
	
cc	if(ltrajcost) then
cc	  k=jtstep-1
cc	  if((mod(k,nsamplestep).eq.0)
cc   &.and.(k.ge.((ntrajstart*3)/2)).and.(k.lt.ntstep)) then
cc	    index=( k-((ntrajstart*3)/2) )/nsamplestep
cc        index=0

cc          fac=1.
cc	    call mult3(etraj(1,1,index),fac,bs)
cc	    call add3(dadnew,bs,dadnew)
cc	  endif
cc	endif
	
 200  CONTINUE

C
      CALL COPY(DADNEW,BS)
      
cc    write(*,*) 'lec=',lec
cc    if (lec) then
c see note in timeNL_forc_knmi.f
c see also OPK

        do i=1,nsh2
          do j=1,3
            cs(i,j)=0.0
          enddo
        enddo
cc	do kforce=1,numforce
cc	  call copy(cs,multiforcing(1,1,kforce))
cc	enddo
      
cc      do jtstep=ntstep+1,1,-1
cc        if (jtstep.eq.1) fac=dt*pigr4*2.0d0/4.
cc        if (jtstep.eq.2) fac=dt*pigr4*2.0d0/2.
cc        if (jtstep.gt.2) fac=dt*pigr4*2.0d0
cc	  
cc	  call mult3(forcingtraj(1,1,jtstep-1),fac,
cc   &               forcingtraj(1,1,jtstep-1))
cc      enddo
cc	 
cc      do jtstep=ntstep+1,1,-1 
cc	  do kforce=1,numforce
cc	    if (  ( jtstep.le.(nforceswitch(kforce+1)))
cc   &    .and.(jtstep.gt.(nforceswitch(kforce))) ) then
cc            call add3(multiforcing(1,1,kforce),
cc   &     forcingtraj(1,1,jtstep-1),multiforcing(1,1,kforce))                    
cc          endif
cc	  enddo
cc	enddo 	
	
cc	call add3(multiforcing(1,1,1),
cc   &     forcingtraj(1,1,0),multiforcing(1,1,1))
          
cc    endif

      RETURN
      END



