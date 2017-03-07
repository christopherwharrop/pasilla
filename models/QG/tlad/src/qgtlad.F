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
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'

      REAL*8 XOLD(NSH2,3),XNEW(NSH2,3),DXNEW(NSH2,3)
      
      real*8 as(nsh2,3)
      real*8 forcing,forcingpsi
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4

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
        DTT=2.*DTT
C
      ENDIF
C
      RETURN
      END

c      SUBROUTINE DDTLAD(DS,DSFM)
      SUBROUTINE DDTLAD(DS,DSFM,jstep)
C *** Computes the adjoint of the tangent linear equations
C *** stripped down efficient version
C *** assumes a correct version of PSI 
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      include 'comtes.f'
      
      COMMON /CLINEA/ LITEST(10)
      LOGICAL LITEST
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4
      logical lforcwr
      common /wrforc/ lforcwr
      integer jstep
      integer nstop
      common /ntime/ nstop
      real*8 forcing,forcingpsi
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)

      REAL*8 CP(NSH2),WSPACE(NSH2),WS1(NSH2),DSMT(NSH2,2),
     &  AS(NSH2,3),CS(NSH2,3),DS(NSH2,3),DSM(NSH2,3),
     &  DSMLAP(NSH2,3),DSFM(NSH2,3),LPSI(NSH2,3),SJACOB(NSH2),
     &  DSMLGA(NGP)
      REAL*8 PIGR4,e
      
      integer nnnstep
      
      if (lforc3) then
      
        if (jstep.gt.ntadj) then
         write(*,*) 'too many forcing steps'
         stop
        endif
	
        call copy(forcingtraj(1,1,jstep-1),forcingpsi)

        PIGR4=8.0D0*ASIN(1.0D0)
	DT=1./36.
        fac=1.0
	
c have a look in the QG paper:
c the identity operator is only active during the application
c time of the forcing (in which case fac=0)
	
        do i=1,3
          do j=1,nsh2
            do kforce=1,numforce
              if (kforce.eq.nforce) then
                if ( jstep.eq.(nforceswitch(kforce+1)+1)) then
		  if(i*j.eq.1) write(*,*) 'nullify',jstep
                  dadnew(j,i)=0.0
	          dadold(j,i)=0.0
                endif
		if (  ( jstep.le.(nforceswitch(kforce+1))).and.
     &              (jstep.gt.(nforceswitch(kforce))) ) then
                  fac=0.0
		endif
		if ((kforce.eq.1).and.(jstep.eq.1)) fac=0.0  
		forcingpsi(j,i)=(1.0-fac)*forcingpsi(j,i)
		
              endif
            enddo

	    DSFM(j,i)=forcingpsi(j,i)
	  enddo 
        enddo
	
      else

        if(lforcwr) then
c        write(*,*) 'write',jstep-1
        call copy(ds,forcingtraj(1,1,jstep-1))
        endif



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
        DO 253,I=1,NGP
 253      GJACOB(I)=RDISS(I)*DSMLGA(I)
        CALL GGTOSP(GJACOB,SJACOB)
        DO 252,K=1,NSH2
 252      WSPACE(K)=DISS(K,2)*SJACOB(K)
        CALL LAPINV(WSPACE,SJACOB)
        CALL ADDAB(DSFM(1,3),SJACOB)
C    Adjoint of the L-derivative
        DO 254,I=1,NGP
 254      GJACOB(I)=DDISDX(I)*DSMLGA(I)
        CALL GGTOSP(GJACOB,WSPACE)
        CALL DDL(WSPACE,SJACOB)
        CALL LAPINV(SJACOB,WSPACE)
        CALL ADDAB(DSFM(1,3),WSPACE)
C    Adjoint of the M-derivative
        CALL MUDERA(DSM(1,3),SJACOB)
        CALL ADDAB(DSFM(1,3),SJACOB)


      endif

      RETURN
      END

      SUBROUTINE DDTAD(JTSTEP)
C THE ADJOINT OF DDT
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      include 'comtes.f'

      REAL*8 QPRIM(NSH23),DQPRIM(NSH23),as(nsh2,3)
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4
      integer nstop
      common /ntime/ nstop

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
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      include 'comtes.f'
      real*8 forcing,forcingpsi
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
      real*8 amp
      common amp
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4
      real*8 as(nsh2,3),zz(nsh2,3),uu(nsh2,3)
      common /ntime/ nstop
      integer nstop

      IF (JTSTEP.GT.2) THEN
C
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

        PIGR4=8.0D0*ASIN(1.0D0)
        DT=1./36.
c**        DTT=DT*PIGR4*2.0D0
        IF (JTSTEP.EQ.1) DTT=DT*PIGR4*2.0D0/4.
        IF (JTSTEP.EQ.2) DTT=DT*PIGR4*2.0D0/2.
C
        IF (JTSTEP.EQ.1) THEN
          DO 200 L=1,3
           DO 200 K=1,NSH2
           DADNEW(K,L)=DADOLD(K,L)+DADNEW(K,L)+DTT*DDADNE(K,L)
  200     continue

        ELSE IF (JTSTEP.EQ.2) THEN
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
c      write(*,*) 'step in adjoint model =',jtstep,dtt
      RETURN
      END

      SUBROUTINE DDMU(EPS,DEPSDM)
C *** eps is a field in spherical harmonics
C *** depsdm is the mu-derivative in spherical harmonics
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'

      REAL*8 EPS(NSH2),DEPSDM(NSH2),DEPSGA(NGP)

C *** compute the derivate on the GAUSS-grid      
      CALL SPTOGG(EPS,DEPSGA,PD)
C *** backtransformation to spherical harmonics
      CALL GGTOSP(DEPSGA,DEPSDM)
      
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
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      
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

      SUBROUTINE MULDMV(DM,EPS,DMEPS,LFM,LLM,LM)
C *** compute the mu-derivative DMEPS of EPS 
C *** this call must be preceded by a call to LINDMU for
C *** the computation of the auxiliary matrices 
C *** DM, LFM, LLM and LM
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      
      INTEGER LFM(0:NM),
     &  LLM(0:NM),LM(NSH)
      REAL*8  ZERO(NSH2),EPS(NSH2),DMEPS(NSH2),
     &  DM(NSH,0:NM)
     
      DO 1,I=1,NSH2
  1     DMEPS(I)=0.0D0
      DO 2,K=1,NSH
        M=LM(K)
        N=LL(K)
        KMIN=LFM(M)
        KMAX=LLM(M)
        DO 3,KF=KMIN,KMAX
 3        DMEPS(KF)=DMEPS(KF)+DM(KF,N)*EPS(K)
 2    CONTINUE
      DO 4,K=NM+2+NSH,NSH2
        M=LM(K-NSH)
        N=LL(K-NSH)
        KMIN=LFM(M)
        KMAX=LLM(M)
        DO 5,KF=KMIN,KMAX
          KFNSH=KF+NSH
          DMEPS(KFNSH)=DMEPS(KFNSH)+DM(KF,N)*EPS(K)
 5      CONTINUE
 4    CONTINUE
      RETURN
      END

      SUBROUTINE MUDMTV(DM,EPS,DMTEPS,LFM,LLM,LM,EIEI)
C *** compute the adjoint of the mu-derivative of EPS 
C *** the result is stored in DMTEPS
C *** this call must be preceeded by a call to LINDMU for
C *** the computation of DM, LFM, LLM, LM and EIEI
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      
      INTEGER LFM(0:NM),
     &  LLM(0:NM),LM(NSH)
      REAL*8  ZERO(NSH2),EPS(NSH2),DMTEPS(NSH2),
     &  DM(NSH,0:NM),EIEI(NSH)
     
      DO 1,I=1,NSH2
  1     DMTEPS(I)=0.0D0
      DO 2,K=1,NSH
        M=LM(K)
        N=LL(K)
        KMIN=LFM(M)
        KMAX=LLM(M)
        DO 3,KF=KMIN,KMAX
C NO ENERGY IN THE FIRST MODE the (0,0)-mode
          IF (EIEI(K).GT.001) THEN
            DMTEPS(K)=DMTEPS(K)+DM(KF,N)*EPS(KF)*
     &      EIEI(KF)/EIEI(K)
          ENDIF
 3      CONTINUE
 2    CONTINUE
      DO 4,K=NM+2+NSH,NSH2
        M=LM(K-NSH)
        N=LL(K-NSH)
        KMIN=LFM(M)
        KMAX=LLM(M)
        DO 5,KF=KMIN,KMAX
          KFNSH=KF+NSH
          IF (EIEI(K-NSH).GT.0.001) THEN
            DMTEPS(K)=DMTEPS(K)+DM(KF,N)*EPS(KFNSH)*
     &      EIEI(KF)/EIEI(K-NSH)
          ENDIF
 5      CONTINUE
 4    CONTINUE
      RETURN
      END
      
c      SUBROUTINE DDTL(DS,DEPSDT)
      SUBROUTINE DDTL(DS,DEPSDT,jstep)
C *** COMPUTES THE TANGENT LINEAR EQUATIONS
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'comfra.f'
      include 'comtes.f'
      COMMON /CLINEA/ LITEST(10)
      LOGICAL LITEST
      integer jstep,nnnstep
      REAL*8 CP(NSH2),WSPACE(NSH2),
     &  AS(NSH2,3),BS(NSH2,3),DS(NSH2,3),
     &  DEPSDT(NSH2,3),DUMMY,DFEPS(NSH2,3),DELPST(NSH2,2),
     &  DELPSI(NSH2,3),LPSI(NSH2,3)

      real*8 xs(nsh2,3),fac
      real*8 forcing,forcingpsi
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4
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
      
      IF (LITEST(1)) THEN
        CALL SUB (AS(1,1),CP,BS(1,1))  
        CALL JACOB(PSI(1,1),BS(1,1),WSPACE(1))
        CALL ADDAB (DFEPS(1,1),WSPACE(1))
       ENDIF
      
      IF (LITEST(2)) THEN
        CALL JACOB(DELPSI(1,1),LPSI(1,1),BS(1,1))
        CALL ADDAB (DFEPS(1,1),BS(1,1))
      ENDIF
      
      IF (LITEST(3)) THEN
        CALL ADDAB (DFEPS(1,1),DELPST(1,1))
      ENDIF 
      
      IF (LITEST(4)) THEN
        CALL SUB (AS(1,2),CP,BS(1,2))  
        CALL JACOB(PSI(1,2),BS(1,2),WSPACE(1))
        CALL ADDAB (DFEPS(1,2),WSPACE(1))
      ENDIF
      
      IF (LITEST(5)) THEN
        CALL JACOB(DELPSI(1,2),LPSI(1,2),BS(1,2))
        CALL ADDAB (DFEPS(1,2),BS(1,2))
      ENDIF

      IF (LITEST(6)) THEN
        CALL SUBAB (DFEPS(1,2),DELPST(1,1))
        CALL ADDAB (DFEPS(1,2),DELPST(1,2))
      ENDIF
      
      IF (LITEST(7)) THEN
        CALL SUB (AS(1,3),CP,BS(1,3))  
        CALL JACOB(PSI(1,3),BS(1,3),WSPACE(1))
        CALL ADDAB (DFEPS(1,3),WSPACE(1))
      ENDIF

      IF (LITEST(8).OR.LITEST(10)) THEN
        CALL ADD(LPSI(1,3),OROG,WSPACE(1))
        CALL JACOBP(DELPSI(1,3),WSPACE(1),BS(1,3))
        CALL ADDAB (DFEPS(1,3),BS(1,3))
      ENDIF
      
      IF (LITEST(9)) THEN
        CALL SUBAB (DFEPS(1,3),DELPST(1,2))
      ENDIF

      CALL QPSI(DFEPS,DEPSDT)

c pick appropriate forcing

      if (lmultiforce) then
        do k=1,numforce
	  
          if (jstep.eq.(nforceswitch(k))) then  
	   write(*,*) 'activate forcing',k,'  at step',nforceswitch(k)
	    do i=1,3
              do j=1,nsh2
                forcingpsi(j,i)=multiforcing(j,i,k)
              enddo
            enddo
	  endif
	enddo
      endif

      
      if (lforc4) then
        do i=1,3
          do j=1,nsh2
	      DEPSDT(j,i)=DEPSDT(j,i)+forcingpsi(j,i)
          enddo
        enddo
      endif
      
      RETURN
      END

      FUNCTION RAN1(IDUM)
C *** Taken from the numerical recipes
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1)PAUSE
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END

      SUBROUTINE ADDAB(AA,BB)
C
C***  computes the sum of AA and BB, 
C***  the result is put in AA.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2)
C
      DO 10 K=1,NSH2
        AA(K)=AA(K)+BB(K)
   10 CONTINUE
C
      RETURN
      END

      SUBROUTINE SUBAB(AA,BB)
C
C***  subtracts BB from AA, 
C***  the result is put in AA.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2)
C
      DO 10 K=1,NSH2
        AA(K)=AA(K)-BB(K)
   10 CONTINUE
C
      RETURN
      END
      
      SUBROUTINE NEG(AA)
C *** Negate AA
      include 'comfra.f'

      REAL*8 AA(NSH2)
      
      DO 10 K=1,NSH2
        AA(K)=-AA(K)
 10   CONTINUE
 
      RETURN
      END
      
      
      SUBROUTINE READCI
C *** Initialize the common CLINEA which is used to select
C *** part of the differential equation only
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CLINEA/ LITEST(10)
      LOGICAL LITEST
      INTEGER IC
      CHARACTER*7 CDUMMY
      
      OPEN(UNIT=21,FILE='CIETJES')
      DO 1 I=1,10
        READ(21,999) CDUMMY,IC
        WRITE(*,*) 'I AND IC: ',I,IC
        IF (IC.EQ.1) THEN
          LITEST(I)=.TRUE.
        ELSE
          LITEST(I)=.FALSE.
        ENDIF
 1    CONTINUE
 999  FORMAT(A7,I1)
      
      CLOSE(UNIT=21)
      
      RETURN
      END
      
      
      SUBROUTINE DDTP
C
C *** Computation of time derivative of the potential vorticity fields
C *** Adapted to be able to test the tangent linear equations
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /CLINEA/ LITEST(10)
      LOGICAL LITEST

      include 'comfra.f'

      real*8 forcing,forcingpsi
      logical lforcnl
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
      common /addforc2/ lforcnl

      DO 400 K=1,3
        DO 400 L=1,NSH2
 400  DQPRDT(L,K)=0.0D0
C *** test the consistency of the logicals
      IF ((LITEST(1).NEQV.LITEST(2)).OR.
     &    (LITEST(4).NEQV.LITEST(5)).OR.
     &    (LITEST(7).NEQV.LITEST(8))) THEN
        WRITE(*,*) 'PECULIAR COMBINATION OF INPUT LOGICALS'
        WRITE(*,*) 'PROGAM TERMINATED IN DDTP '
        PAUSE 'PRESS RETURN TO STOP'
        STOP
      ENDIF
      
C
C *** Advection of potential vorticity at upper level
C
      IF (LITEST(1)) THEN
        WRITE(*,*) 'LITEST(1)=TRUE JACOBIAN FOR LOWER LEVEL'
        CALL JACOB (PSI(1,1),QPRIME(1,1),DQPRDT(1,1))
      ENDIF
C
C *** Advection of potential vorticity at middle level
C
      IF (LITEST(4)) THEN
        WRITE(*,*) 'LITEST(4)=TRUE JACOBIAN FOR MIDDLE LEVEL'
        CALL JACOB (PSI(1,2),QPRIME(1,2),DQPRDT(1,2))
      ENDIF
C
C *** Advection of potential vorticity and dissipation at lower level
C
      DO 100 K=1,NSH2
  100 WS(K)=QPRIME(K,3)+OROG(K)
C
      IF (LGDISS) THEN
C
        IF (LITEST(7).OR.LITEST(10)) THEN
          WRITE(*,*) 'CALLING JACOBP FOR UPPER LEVEL'
          CALL JACOBP (PSI(1,3),WS,DQPRDT(1,3))
        ENDIF
C
      ELSE
C
        WRITE(*,*) 'DDTP: LGDISS UNEXPECTEDLY FALSE'
        WRITE(*,*) 'PROGRAM TERMINATED'
        PAUSE 'PRESS RETURN TO STOP'
        STOP
C        CALL JACOB (PSI(1,3),WS,DQPRDT(1,3))
C
C        DO 150 K=1,NSH2
C  150   DQPRDT(K,3)=DQPRDT(K,3)+DISS(K,2)*PSI(K,3)
C
      ENDIF
C
C *** Relaxation of temperature and forcing
C

      if (lforcnl) then
      DO K=1,NSH2
      DO L=1,3
        FOR(K,L,1)=FOR(K,L,1)+forcing(K,L)
      ENDDO
      ENDDO
      end if

      WRITE(*,*) 'RELAXATION  AND OR FORCING'
      DO 200 K=1,NSH2
        PSIT(K,1)=RELT1*PSIT(K,1)
        PSIT(K,2)=RELT2*PSIT(K,2)
        IF (LITEST(3)) THEN
          DQPRDT(K,1)=DQPRDT(K,1)+PSIT(K,1)          +FOR(K,1,1)
        ELSE
          DQPRDT(K,1)=DQPRDT(K,1) + FOR(K,1,1)
        ENDIF
        IF (LITEST(6)) THEN
          DQPRDT(K,2)=DQPRDT(K,2)-PSIT(K,1)+PSIT(K,2)+FOR(K,2,1)
        ELSE
          DQPRDT(K,2)=DQPRDT(K,2)+FOR(K,2,1)
        ENDIF
        IF (LITEST(9)) THEN
          DQPRDT(K,3)=DQPRDT(K,3)          -PSIT(K,2)+FOR(K,3,1)
        ELSE
          DQPRDT(K,3)=DQPRDT(K,3)+FOR(K,3,1)
        ENDIF
  200 CONTINUE
C
C *** Explicit horizontal diffusion
C
      IF (LEXDIF) THEN
C
        WRITE(*,*) 'DDTP: LEXDIF UNEXPECTEDLY TRUE'
        WRITE(*,*) 'PROGRAM TERMINATED'
        PAUSE 'PRESS RETURN TO STOP'
        STOP
C        DO 300 L=1,3
C        DO 300 K=1,NSH2
C  300   DQPRDT(K,L)=DQPRDT(K,L)+DISS(K,1)*QPRIME(K,L)
C
      ENDIF
C
      RETURN
      END
      
      SUBROUTINE INPRSQ(XS,YS,E)
C Compute the squared norm inner product
      IMPLICIT NONE
      include 'comfra.f'
      REAL*8 XS(NSH2,3),YS(NSH2,3),E
      INTEGER L,K

      E=0.0
      DO 10,L=1,3
       DO 20,K=2,NM+1
	E=E+XS(K,L)*YS(K,L)
 20    CONTINUE
       DO 30,K=NM+2,NSH
	E=E+2*(XS(K,L)*YS(K,L)+XS(K+NSH,L)*YS(K+NSH,L))
 30    CONTINUE
 10   CONTINUE

      RETURN
      END


      SUBROUTINE ADJSQ(AS,BS,NUUR)
C this subroutine integrates the
C adjoint equations with initial condition
C AS and with integration time NUUR
C the result is BS.
C the adjoint here is determined with respect to the 
C squared norm inner product
      IMPLICIT NONE
      include 'comfra.f'
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


      SUBROUTINE ADD(AA,BB,CC)
C
C***  computes the sum of AA and BB, 
C***  the result is put in CC.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2),CC(NSH2)
C
      DO 10 K=1,NSH2
      CC(K)=AA(K)+BB(K)
   10 CONTINUE
C
      RETURN
      END
      
      SUBROUTINE ADD3(AA,BB,CC)
C
C***  computes the sum of AA and BB, 
C***  the result is put in CC.
C
      include 'comfra.f'

      REAL*8 AA(NSH2,3),BB(NSH2,3),CC(NSH2,3)
C
      DO L=1,3
        DO K=1,NSH2
          CC(K,L)=AA(K,L)+BB(K,L)
        ENDDO
      ENDDO
C
      RETURN
      END
      
      
      
      SUBROUTINE SUB(AA,BB,CC)
C
C***  subtracts BB from AA, 
C***  the result is put in CC.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2),CC(NSH2)
C
      DO 10 K=1,NSH2
      CC(K)=AA(K)-BB(K)
   10 CONTINUE
C
      RETURN
      END
      
           
      SUBROUTINE SUB3(AA,BB,CC)
C
C***  subtracts BB from AA, 
C***  the result is put in CC.
C
      include 'comfra.f'

      REAL*8 AA(NSH2,3),BB(NSH2,3),CC(NSH2,3)
C
      DO L=1,3
        DO K=1,NSH2
          CC(K,L)=AA(K,L)-BB(K,L)
        ENDDO
      ENDDO
C
      RETURN
      END   
      
      SUBROUTINE MULT(AA,C,BB)
C
C***  computes the product of AA with real C, 
C***  the result is put in BB.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2),C
C
      
      DO 10 K=1,NSH2
      BB(K)=AA(K)*C
   10 CONTINUE
C
    
      RETURN
      END
      
      SUBROUTINE MULT3(AA,C,BB)
C
C***  computes the product of AA with real C, 
C***  the result is put in BB.
C
      include 'comfra.f'

      REAL*8 AA(NSH2,3),BB(NSH2,3),C
C
      DO L=1,3
        DO K=1,NSH2
          BB(K,L)=AA(K,L)*C
        ENDDO
      ENDDO
C
    
      RETURN
      END
      
      SUBROUTINE COP(AA,BB)
C
C** copies AA(1 level) TO BB(1 level)
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2)
C
      DO 10 K=1,NSH2
      BB(K)=AA(K)
   10 CONTINUE
C
      RETURN
      END

      SUBROUTINE COPY(AA,BB)
C
C***  copies AA into BB.
C
      include 'comfra.f'

      REAL*8 AA(NSH2,3),BB(NSH2,3)
C
      DO 10 L=1,3
      DO 20 K=1,NSH2
      BB(K,L)=AA(K,L)
   20 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      
      
      
      SUBROUTINE INPROD(AA,BB,E)
C
C***  computes the energy innerproduct of AA and BB
C***  the result is put in E.
C
      include 'comfra.f'

      REAL*8 AA(NSH2,3),BB(NSH2,3),E
C
      E=0.
      DO 10 L=1,3
C
      DO 20 K=2,NM+1
      E=E+LL(K)*(LL(K)+1)*AA(K,L)*BB(K,L)
   20 CONTINUE
C
      DO 30 K=NM+2,NSH
      E=E+2*LL(K)*(LL(K)+1)*AA(K,L)*BB(K,L)
     *+2*LL(K)*(LL(K)+1)*AA(K+NSH,L)*BB(K+NSH,L)
   30 CONTINUE
C   
   10 CONTINUE
C
      RETURN
      END
      
      
      SUBROUTINE INPROSN(AA,BB,E)
C
C***   computes the SN inner product
C***  the result is put in E.
C
      include 'comfra.f'

      REAL*8 AA(NSH2,3),BB(NSH2,3),E
C
      E=0.
      DO 10 L=1,3
C
      DO 20 K=2,NM+1
      E=E+AA(K,L)*BB(K,L)
   20 CONTINUE
C
      DO 30 K=NM+2,NSH
      E=E+2*AA(K,L)*BB(K,L)
     *+2*AA(K+NSH,L)*BB(K+NSH,L)
   30 CONTINUE
C
   10 CONTINUE
C
      RETURN
      END

      SUBROUTINE LENGTH(XS,E)
C
C**   computes the squared norm of XS
      include 'comfra.f'

      REAL*8 XS(NSH2,3),E
      E=0.0d0
      DO 10 L=1,3
C
      DO 20 K=2,NM+1
      E=E+XS(K,L)**2
   20 CONTINUE
C
      DO 30 K=NM+2,NSH
      E=E+2*(XS(K,L)**2+XS(K+NSH,L)**2)
   30 CONTINUE
C
   10 CONTINUE
C
      RETURN
      END


      
      SUBROUTINE ENER(XS,E)
C**   COMPUTES THE KINETIC ENERGY AT
C**   ONE LEVEL.
      include  'comfra.f'
      REAL*8 XS(NSH2),E

      E=0.
      DO 20 K=2,NM+1
      E=E+LL(K)*(LL(K)+1)*XS(K)**2
   20 CONTINUE
      DO 30 K=NM+2,NSH
      E=E+2*LL(K)*(LL(K)+1)*XS(K)**2
     *  +2*LL(K)*(LL(K)+1)*XS(K+NSH)**2
   30 CONTINUE
      RETURN
      END

      
      
      
      SUBROUTINE LAP(AA,BB)
C
C***  computes the Laplacian of AA, the result is put in BB.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2)
C
      BB(1)=0.
      BB(NSH+1)=0.
      DO 10 K=2,NSH
      BB(K)=AA(K)/RINHEL(K,1)
      BB(K+NSH)=AA(K+NSH)/RINHEL(K+NSH,1)
   10 CONTINUE
C
      RETURN
      END
      
      
      
      SUBROUTINE LAP3(XS,YS)
C
C***  computes the Laplacian of XS, the result is put in YS.
C
      include 'comfra.f'

      REAL*8 XS(NSH2,3),YS(NSH2,3)
C
      DO 1 L=1,3
      YS(1,L)=0.
      YS(NSH+1,L)=0.
      DO 10 K=2,NSH
      YS(K,L)=XS(K,L)/RINHEL(K,1)
      YS(K+NSH,L)=XS(K+NSH,L)/RINHEL(K+NSH,1)
10    CONTINUE
1     CONTINUE
C
      RETURN
      END

            
      SUBROUTINE LAPINV(AA,BB)
C
C***  computes the inverse Laplacian of AA
C***  the result is put in BB.
C
      include 'comfra.f'

      REAL*8 AA(NSH2),BB(NSH2)
C
      BB(1)=0.
      BB(NSH+1)=0.
      DO 10 K=2,NSH
      BB(K)=RINHEL(K,1)*AA(K)
      BB(K+NSH)=RINHEL(K+NSH,1)*AA(K+NSH)
   10 CONTINUE
C
      RETURN
      END
      
      
      SUBROUTINE LAPINV3(XS,YS)
C
C***  computes the inverse Laplacian of XS
C***  the result is put in YS.
C
      include 'comfra.f'

      REAL*8 XS(NSH2,3),YS(NSH2,3)
C
      DO 1 L=1,3
      YS(1,L)=0.
      YS(NSH+1,L)=0.
      DO 10 K=2,NSH
      YS(K,L)=RINHEL(K,1)*XS(K,L)
      YS(K+NSH,L)=RINHEL(K+NSH,1)*XS(K+NSH,L)
10    CONTINUE
1     CONTINUE
C
      RETURN
      END


      SUBROUTINE MUDERA (XS,YS)
C***  computes the 'adjoint' of the 'mu-derivative'
C***  input is XS, output is YS.
      
      include 'comfra.f'
      
      REAL*8 XS(NSH2),YS(NSH2),ZS(NSH2)
      
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
      

      SUBROUTINE TANG(AS,BS,NUUR)
C**   this subroutine integrates the tangent
C**   linear equation with initial condition 
C**   AS(NSH2,3) and integration time NUUR. 
C**   It uses the orbit stored by STOREP.
C**   The result is BS(NSH2,3)
C
      IMPLICIT NONE

      INTEGER NLMAX,NPERT,NSKIP,NOUTST,IDIF,NPRINT
      PARAMETER (nlmax=800)
      REAL*8 RMSV,DT,EPS,TDIS,ADDISL,ADDISH,TREL,TDIF,H0,
     *  RRDEF1,RRDEF2
      PARAMETER ( NPERT=0, NSKIP=0, RMSV=1.E-02 )
      PARAMETER ( NOUTST=36, DT=1./36., EPS=.03,
     *            TDIS=3.0, ADDISL=0.5, ADDISH=0.5,
     *            TREL=25., TDIF=2.0, IDIF=4,
     *            H0=9., RRDEF1=.110, RRDEF2=.070 )
C
C      INCLUDE COMFRA.F
C      INCLUDE COMTES.F
      include 'comfra.f'
      include 'comtes.f'

      REAL*8 X(NLAT),W(NLAT),RAD(NLAT),PIGR4,TREAL
      REAL*8 DUM(NSH23),DUMOLD(NSH23),DUMNEW(NSH23)
      REAL*8 AS(NSH2,3),BS(NSH2,3),cs(nsh2,3),e,fac
      real*8 forcing,forcingpsi
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4
      integer janstep,k,nday,index
      common /stap/ janstep
      integer nstop
      common /ntime/ nstop
      
      real*8 zt(nlat,nlon),dummy,sum
      integer ix,iy,ilevel
      INTEGER LM(NSH)
C
      INTEGER INQ,INP,INH,INF,NTSTEP,NLANCO,JTSTEP,JRUN,
     *  ISEED,I,J,ILANCO,n
      INTEGER NUUR
      NTSTEP=(NUUR*3)/2
      nstop=NTSTEP
      
      if(lforc) then
      write(*,*) 'ntstep=',janstep
      ntstep=janstep
      end if
      
      CALL COPY(AS,DUMNEW)
C
      DO 106,I=1,NSH23
      DUMOLD(I)=0.0D0
  106 CONTINUE
      WRITE(*,*) 'TANGENT MODEL'
       PIGR4=8.0D0*ASIN(1.0D0)
       DTT=DT*PIGR4/2.0D0
       DO 107,JTSTEP=0,NTSTEP
         CALL RETRIE(JTSTEP)
c         CALL DDTL(DUMNEW,DUM)
         CALL DDTL(DUMNEW,DUM,jtstep)
C***
         CALL TIMINP(JTSTEP,DUMOLD,DUMNEW,DUM)
	 
	 
c	 if(ltangwr) then
c	   k=jtstep
c           if((mod(k,nsamplestep).eq.0)
c     &       .and.(k.ge.((ntrajstart*3)/2))) then
c             call inprod(dumnew,dumnew,e)
c	     write(45,*) 'TANG',k,e
c	   endif
c	 endif
	 
	 if(ltangwr) then
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
	 
         if(ltrajcost) then
	 
c fill array ETRAJ to determine costfunction

	   k=jtstep
           if((mod(k,nsamplestep).eq.0)
     &       .and.(k.ge.((ntrajstart*3)/2)).and.(k.lt.ntstep)) then
	     index=( k-((ntrajstart*3)/2) )/nsamplestep
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

      SUBROUTINE ADJ(AS,BS,NUUR) 
C**   this subroutine integrates the adjoint
C**   equation with initial condition AS(NSH2,3)
C**   and with integration time NUUR. It uses the orbit
C**   stored by STOREP.
C**   The result is BS(NSH2,3). The adjoint is determined with
C**   respect to the kinetic energy inner product.
C BY PETER HOUTEKAMER (PH) AND JAN BARKMEIJER (JB)
      IMPLICIT NONE

      INTEGER NLMAX,NPERT,NSKIP,NOUTST,IDIF,NPRINT
      PARAMETER (nlmax=800)
      REAL*8 RMSV,DT,EPS,TDIS,ADDISL,ADDISH,TREL,TDIF,H0,
     *  RRDEF1,RRDEF2
      PARAMETER ( NPERT=0, NSKIP=0, RMSV=1.E-02 )
      PARAMETER ( NOUTST=36, DT=1./36., EPS=.03,
     *            TDIS=3.0, ADDISL=0.5, ADDISH=0.5,
     *            TREL=25., TDIF=2.0, IDIF=4,
     *            H0=9., RRDEF1=.110, RRDEF2=.070 )
C
C      INCLUDE COMFRA.F
C      INCLUDE COMTES.F
      include 'comfra.f'
      include 'comtes.f'

      REAL*8 X(NLAT),W(NLAT),RAD(NLAT),PIGR4,TREAL
      REAL*8 AS(NSH2,3),BS(NSH2,3)
      real*8 forcing,forcingpsi
      common /forc/ forcing(nsh2,3),forcingpsi(nsh2,3)
c      real*8 forcingtraj
c      common forcingtraj(nsh2,3,0:120)
      logical lforc,lforc2,lforc3,lforc4
      common /addforc/ lforc,lforc2,lforc3,lforc4
      logical lforcwr
      common /wrforc/ lforcwr
      integer nstop,janstep,kforce
      common /ntime/ nstop
      common /stap/ janstep
      real*8 uold(nsh2,3),unew(nsh2,3),e,e1,e2,fac,
     &       zold(nsh2,3),znew(nsh2,3),cs(nsh2,3)
      
      INTEGER LM(NSH)
C
      INTEGER INQ,INP,INH,INF,NTSTEP,NLANCO,JTSTEP,JRUN,
     *  ISEED,I,J,ILANCO,m,k,nnstep,index
      INTEGER NUUR
      DATA INQ/ 2/, INP/ 12/, INH/ 3/, INF/ 4/
C
      NTSTEP=(NUUR*3)/2
      nstop=ntstep
      
      if(lforc) then
      write(*,*) 'ntstep',janstep
      ntstep=janstep
      nstop=ntstep
      end if

      if(lforcwr) then
        do i=1,nsh2
          do j=1,3
            do m=0,120
              forcingtraj(i,j,m)=0.0
            enddo
          enddo
        enddo 
      endif
    
      
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
      PIGR4=8.0D0*ASIN(1.0D0)
      DTT=DT*PIGR4*2.0D0
      CALL COPY(AS,DADNEW)

      PIGR4=8.0D0*ASIN(1.0D0)
      DTT=DT*PIGR4*2.0D0

C DADNEW is hier de foutenvektor op tijdstip T
      WRITE(*,*) 'START THE BACKWARD INTEGRATION'
C het terugintegreren met de adjoint
C iets anders bedenken voor ntstep <=1

      DO 200,JTSTEP=NTSTEP+1,1,-1	
	  
C *** QPRIME(=PSI-ORBIT) IS STORED
        CALL RETRIE(JTSTEP-1)
        CALL DDTAD(JTSTEP)
        CALL TIMIAD(JTSTEP)
	
	if(ltrajcost) then
	  k=jtstep-1
	  if((mod(k,nsamplestep).eq.0)
     &.and.(k.ge.((ntrajstart*3)/2)).and.(k.lt.ntstep)) then
	    index=( k-((ntrajstart*3)/2) )/nsamplestep


            fac=1.
	    call mult3(etraj(1,1,index),fac,bs)
	    call add3(dadnew,bs,dadnew)
	  endif
	endif
	
 200  CONTINUE

C
      CALL COPY(DADNEW,BS)
      
      write(*,*) 'lec=',lec
      if (lec) then
c see note in timeNL_forc_knmi.f
c see also OPK

        do i=1,nsh2
          do j=1,3
            cs(i,j)=0.0
          enddo
        enddo
	do kforce=1,numforce
	  call copy(cs,multiforcing(1,1,kforce))
	enddo
      
        do jtstep=ntstep+1,1,-1
          if (jtstep.eq.1) dtt=dt*pigr4*2.0d0/4.
          if (jtstep.eq.2) dtt=dt*pigr4*2.0d0/2.
          if (jtstep.gt.2) dtt=dt*pigr4*2.0d0
          fac=dtt
	  
	  call mult3(forcingtraj(1,1,jtstep-1),fac,
     &               forcingtraj(1,1,jtstep-1))
        enddo
	 
        do jtstep=ntstep+1,1,-1 
	  do kforce=1,numforce
	    if (  ( jtstep.le.(nforceswitch(kforce+1)))
     &    .and.(jtstep.gt.(nforceswitch(kforce))) ) then
              call add3(multiforcing(1,1,kforce),
     &     forcingtraj(1,1,jtstep-1),multiforcing(1,1,kforce))                    
            endif
	  enddo
	enddo 	
	
	call add3(multiforcing(1,1,1),
     &     forcingtraj(1,1,0),multiforcing(1,1,1))
          
      endif

      

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

      INTEGER NLMAX,NPERT,NSKIP,NOUTST,IDIF,NPRINT
      PARAMETER (nlmax=800)
      REAL*8 RMSV,DT,EPS,TDIS,ADDISL,ADDISH,TREL,TDIF,H0,
     *  RRDEF1,RRDEF2
      PARAMETER ( NPERT=0, NSKIP=0, RMSV=1.E-02 )
      PARAMETER ( NOUTST=36, DT=1./36., EPS=.03,
     *            TDIS=3.0, ADDISL=0.5, ADDISH=0.5,
     *            TREL=25., TDIF=2.0, IDIF=4,
     *            H0=9., RRDEF1=.110, RRDEF2=.070 )
C
C      INCLUDE COMFRA.F
C      INCLUDE COMTES.F
      include 'comfra.f'
      include 'comtes.f'

      REAL*8 X(NLAT),W(NLAT),RAD(NLAT),PIGR4,TREAL
      INTEGER LM(NSH)
C
      INTEGER INQ,INP,INH,INF,NTSTEP,NLANCO,JTSTEP,JRUN,
     *  I,J
      INTEGER NUUR
      logical lforc,lforc2,lforc3,lforc4,lforcnl
      common /addforc/ lforc,lforc2,lforc3,lforc4
      integer janstep
      common /stap/ janstep
      common /addforc2/ lforcnl
      

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
      PIGR4=8.0D0*ASIN(1.0D0)
      DTT=DT*PIGR4/2.0D0
      write(*,*) 'dtt =',dtt
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
      include 'comfra.f'

      REAL*8 AS(NSH2,3)
      INTEGER NUUR
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
