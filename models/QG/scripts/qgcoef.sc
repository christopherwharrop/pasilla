#!/bin/bash

# This script produces a legendre function files for the Molteni qg model
# for different triangular truncations

compiler='gfortran'
fflags="-O2  -ffixed-line-length-none -fdefault-real-8"

for ngp in 16 32 64 96 160;
do

cat > gausspoints.f <<==
c23456789012345678901234567890123456789012345678901234567890123456789012      
      PROGRAM  gauleg
      implicit none

      REAL*8     newv
      REAL*8     EPS, M_PI
      PARAMETER (EPS=3.0d-12)       	!EPS is the relative precision
      PARAMETER (M_PI=3.141592654d0)      ! Pi value

      INTEGER    i, j, m
      REAL*8     p1, p2, p3, pp, z, z1
      INTEGER    ngp            ! # of Gauss Points
      PARAMETER (ngp=${ngp})
      REAL*8     xabsc(ngp), weig(ngp)


	   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

	   do i = 1, m				! Loop over the desired roots */

     		z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON method   */
100     	p1 = 1.0d0
        	p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        	z1 = z
        	z = z1 - p1/pp             ! Newton Method  */

        	if (dabs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
      	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      enddo     ! i loop
      
      open(1,file='gausspoints.dat',access='append')
      write(1,'(2F18.10)') real(m),real(m-m)
      do i=m+1,2*m
        write(1,'(2F18.10)') xabsc(i),weig(i)
      enddo
      
      end
==

$compiler $fflags gausspoints.f -o gausspoints
./gausspoints

rm gausspoints.f gausspoints
done

for resol in T21 T42 T63 T106
do

cat > bve.inc <<==1
c *** MN-1 corresponds to the triangular spectral truncation
c *** 2*MG denotes the number of gridpoints in the meridional direction

      INTEGER     MN,MG
==1
if [ ${resol} == "T21" ]; then
cat >> bve.inc <<==2
      PARAMETER   (MN=22,MG=16)
==2
fi
if [ ${resol} == "T42" ]; then
cat >> bve.inc <<==3
      PARAMETER   (MN=43,MG=32)
==3
fi
if [ ${resol} == "T63" ]; then
cat >> bve.inc <<==4
      PARAMETER   (MN=64,MG=48)
==4
fi
if [ ${resol} == "T106" ]; then
cat >> bve.inc <<==5
      PARAMETER   (MN=107,MG=80)
==5
fi
cat >> bve.inc <<==
      REAL        GP,GW,PO,PW,PM,PS,AM,AD

      COMMON/INTG/GP(MG),GW(MG)
      COMMON/LEGN/PO(MG,MN,MN),PW(MG,MN,MN),PM(MG,MN,MN),PS(MG,MN,MN)
      COMMON/COEF/AM(MN),AD(MN)
==

cat > legendre.f <<==
      PROGRAM LEGENDRE
      
C**   READ GAUSS POINTS AND WEIGHTS AND CALCULATE LEGENDRE FUNCTIONS

      IMPLICIT NONE
      INCLUDE'bve.inc'

      CALL CALCULATE
      CALL OUTPUT
      
      END
      
      
c23456789012345678901234567890123456789012345678901234567890123456789012      
      SUBROUTINE CALCULATE
      
C**   THIS ROUTINE READS GAUSS POINTS AND WEIGHTS AND CALCULATES 
C**   LEGENDRE FUNCTIONS

      IMPLICIT NONE
      INCLUDE'bve.inc'
      INTEGER     MGT,I,J,M,N,M1,N1
      REAL        QO(MN,MN),QS(MN,MN)
      REAL        DUMMY,RMGT,SIT,WEIGHT,Z,W

C*    GAUSS POINTS AND WEIGHTS
      OPEN(UNIT=1,FILE='gausspoints.dat',STATUS='OLD')
    4 READ(1,101) RMGT,DUMMY
  101   FORMAT(1X,F18.10,F18.10)
        MGT=RMGT
        IF(MGT.EQ.MG) GOTO 2
        DO 3 J=1,MGT
          READ(1,102) SIT,WEIGHT
c          write(*,*) j,sit
  102     FORMAT(1X,F18.10,F18.10)
    3   CONTINUE
      GOTO 4
    2 DO 5 J=1,MG
        READ(1,103) SIT,WEIGHT
c          write(*,*) j,sit
  103   FORMAT(1X,F18.10,F18.10)
        GP(J)=SIT
        GW(J)=WEIGHT
    5 CONTINUE
      CLOSE(1)
      
      
C*    COEFFICIENTS
      DO 9 I=1,MN
        AM(I)=I-1
        AD(I)=(1-I)*I
    9 CONTINUE
C*    LEGENDRE FUNCTION VALUES
      DO 7 I=1,MG
        Z=GP(I)
        W=GW(I)
        CALL FLEGM(Z,QO,QS,MN)
        DO 8 M=1,MN
          DO 12 N=M,MN
            PO(I,M,N)=QO(M,N)
            PW(I,M,N)=(W*QO(M,N))/2.
            PM(I,M,N)=AM(M)*QO(M,N)
            PS(I,M,N)=QS(M,N)
   12     CONTINUE
    8   CONTINUE
    7 CONTINUE

      RETURN
      END 

c23456789012345678901234567890123456789012345678901234567890123456789012      
      SUBROUTINE FLEGM(X,F,G,MN)
C**   THIS ROUTINE CALCULATES ARRAYS F(I,J) AND G(I,J) OF LEGENDRE
C**   FUNCTIONS AND DERIVATIVES, WITH ARGUMENT X, ORDER M=I-1 AND
C**   DEGREE N=J-1 BY MEANS OF RECURRENCY RELATIONS (NORMALIZATION
C**   ACCORDING TO MACHENHAUER)
C**   WRITTEN BY WIM VERKLEY AND FRANK SELTEN
C**   VERSION 17 augustus 1992
      IMPLICIT NONE

      INTEGER  MN,MIJ,MMN,I,J,M,N,INFO,IUIT
      REAL     F(MN,MN),G(MN,MN),X1,X,AM,AN,F1,F2,F3,F4,F5
      REAL     T1,T2,T3,T4,T5,T6,T7,G1,G2,G3,G4,G5
      COMMON/DEBUG/INFO,IUIT

      X1=X
      MIJ=MN
      MMN=MN-1
C*    CALCULATION OF LEGENDRE FUNCTIONS
      DO 1 I=1,MIJ
        DO 10 J=1,MIJ
          F(I,J)=0.
   10   CONTINUE
    1 CONTINUE
      F(1,1)=1.
      DO 2 I=1,MMN
        M=I-1
        N=I-1
        AM=M
        AN=N
        F1=G1(AM,AN)
        F3=G3(AM,AN)
        F4=G4(AM,AN)
        T1=F1*(2.*AN+1.)*X1*F(I,I)
        F(I,I+1)=T1/(AN-AM+1.)
        T3=F3*(AN-AM+1.)*X1*F(I,I+1)
        T4=F4*(AN+AM+1.)*F(I,I)
        F(I+1,I+1)=(T3-T4)/(SQRT(1.-X1*X1))
        IF (I.EQ.MMN) GOTO 2
        DO 5 J=I+1,MMN
          N=J-1
          AN=N
          F1=G1(AM,AN)
          F2=G2(AM,AN)
          T1=F1*(2.*AN+1.)*X1*F(I,J)
          T2=F2*(AN+AM)*F(I,J-1)
          F(I,J+1)=(T1-T2)/(AN-AM+1.)
    5   CONTINUE
    2 CONTINUE
C*    CALCULATION OF DERIVATIVES
      DO 3 I=1,MIJ
        DO 7 J=1,MIJ
          G(I,J)=0.
    7   CONTINUE
    3 CONTINUE
      DO 4 I=1,MIJ
        M=I-1
        N=I-1
        AM=M
        AN=N
        T5=AN*X1*F(I,I)
        G(I,I)=(-1.*T5)/(1.-X1*X1)
        IF (I.EQ.MIJ) GOTO 4
        DO 6 J=I+1,MIJ
          N=J-1
          AN=N
          F5=G5(AM,AN)
          T6=AN*X1*F(I,J)
          T7=F5*(AN+AM)*F(I,J-1)
          G(I,J)=(T7-T6)/(1.-X1*X1)
    6   CONTINUE
    4 CONTINUE
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012      
      FUNCTION G1(AM,AN)
C**   THIS FUNCTION CALCULATES THE FACTOR F1 IN THE ROUTINE FLEGM
C**   VERSION 17 augustus 1992
      IMPLICIT NONE
      REAL     AM,AN,AM1,AN1,A2,A3,B2,B4,G1
      AM1=AM
      AN1=AN
      A2=2.*AN1+1.
      A3=2.*(AN1+1.)+1.
      B2=AN1-AM1+1.
      B4=AN1+AM1+1.
      G1=SQRT((A3*B2)/(A2*B4))
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012      
      FUNCTION G2(AM,AN)
C**   THIS FUNCTION CALCULATES THE FACTOR F2 IN THE ROUTINE FLEGM
C**   VERSION 17 augustus 1992
      IMPLICIT NONE
      REAL     AM,AN,AM1,AN1,A1,A3,B1,B2,B3,B4,G2
      AM1=AM
      AN1=AN
      A1=2.*(AN1-1.)+1.
      A3=2.*(AN1+1.)+1.
      B1=AN1-AM1
      B2=AN1-AM1+1.
      B3=AN1+AM1
      B4=AN1+AM1+1.
      G2=SQRT((A3*B2*B1)/(A1*B4*B3))
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012      
      FUNCTION G3(AM,AN)
C**   THIS FUNCTION CALCULATES THE FACTOR F3 IN THE ROUTINE FLEGM
C**   VERSION 17 augustus 1992
      IMPLICIT NONE 
      REAL     AM,AN,AM1,AN1,B2,B5,G3
      AM1=AM
      AN1=AN
      B2=AN1-AM1+1.
      B5=AN1+AM1+2.
      G3=-1.*SQRT(1./(B5*B2))
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012      
      FUNCTION G4(AM,AN)
C**   THIS FUNCTION CALCULATES THE FACTOR F4 IN THE ROUTINE FLEGM
C**   VERSION 17 augustus 1992
      IMPLICIT NONE
      REAL     AM,AN,AM1,AN1,A2,A3,B4,B5,G4
      AM1=AM
      AN1=AN
      A2=2.*AN1+1.
      A3=2.*(AN1+1.)+1.
      B4=AN1+AM1+1.
      B5=AN1+AM1+2.
      G4=-1.*SQRT(A3/(A2*B5*B4))
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012      
      FUNCTION G5(AM,AN)
C**   THIS FUNCTION CALCULATES THE FACTOR F5 IN THE ROUTINE FLEGM
C**   VERSION 17 augustus 1992
      IMPLICIT NONE
      REAL     AM,AN,AM1,AN1,A1,A2,B1,B3,G5
      AM1=AM
      AN1=AN
      A1=2.*AN1-1.
      A2=2.*AN1+1.
      B1=AN1-AM1
      B3=AN1+AM1
      G5=SQRT((A2*B1)/(A1*B3))
      RETURN
      END

c23456789012345678901234567890123456789012345678901234567890123456789012      
      SUBROUTINE OUTPUT
      IMPLICIT NONE
      INCLUDE'bve.inc'

      INTEGER      M,N,I,J
      REAL         teken0,tekenS
      CHARACTER*2  ft2
      CHARACTER*3  ft3
      
      if (mn.lt.100) then
        write(ft2,'(I2.2)') mn-1
        OPEN(11,FILE='qgcoefT'//ft2//'.dat',FORM='FORMATTED')
      else
        write(ft3,'(I3.3)') mn-1
        OPEN(11,FILE='qgcoefT'//ft3//'.dat',FORM='FORMATTED')
      endif
      
c      m=2
c      n=3
c      do i=MG,1,-1
c        write(*,'(I3,F12.5,2I3,3F12.5)') 
c     *  i+MG,gp(i),m,n,PO(i,m+1,n+1),PS(i,m+1,n+1),PW(i,m+1,n+1)
c      enddo
c      do i=1,MG
c        write(*,'(I3,F12.5,2I3,3F12.5)') 
c     *  MG+1-i,gp(i),m,n,-PO(i,m+1,n+1),PS(i,m+1,n+1),-PW(i,m+1,n+1)
c      enddo
c      write(*,*)
c      m=1
c      n=3
c      do i=MG,1,-1
c        write(*,'(I3,F12.5,2I3,3F12.5)') 
c     *  i+MG,gp(i),m,n,PO(i,m+1,n+1),PS(i,m+1,n+1),PW(i,m+1,n+1)
c      enddo
c      do i=1,MG
c        write(*,'(I3,F12.5,2I3,3F12.5)') 
c     *    MG+1-i,gp(i),m,n,PO(i,m+1,n+1),-PS(i,m+1,n+1),PW(i,m+1,n+1)
c      enddo
      
      
      do i=MN,1,-1
        write(11,*) i
      enddo
      
      do m=0,MN-1
        do n=m,MN-1
          write(11,*) n
        enddo
      enddo
      
      do m=0,MN-1
        do n=m,MN-1
          if (mod(n-m,2).eq.1) then
            teken0=-1d0
            tekenS=1d0
          else
            teken0=1d0
            tekenS=-1d0
          endif
          do j=MG,1,-1
            write(11,*) teken0*PO(j,m+1,n+1)
          enddo
          do j=1,MG
            write(11,*) PO(j,m+1,n+1)
          enddo
        enddo
      enddo
      do m=0,MN-1
        do n=m,MN-1
          if (mod(n-m,2).eq.1) then
            teken0=-1d0
            tekenS=1d0
          else
            teken0=1d0
            tekenS=-1d0
          endif
          do j=MG,1,-1
            write(11,*) tekenS*PS(j,m+1,n+1)
          enddo
          do j=1,MG
            write(11,*) PS(j,m+1,n+1)
          enddo
        enddo
      enddo
      do m=0,MN-1
        do n=m,MN-1
          if (mod(n-m,2).eq.1) then
            teken0=-1d0
            tekenS=1d0
          else
            teken0=1d0
            tekenS=-1d0
          endif
          do j=MG,1,-1
            write(11,*) teken0*PW(j,m+1,n+1)
          enddo
          do j=1,MG
            write(11,*) PW(j,m+1,n+1)
          enddo
        enddo
      enddo
      
      close(11)
      
      END
      
==

$compiler $fflags legendre.f -o legendre
./legendre

rm legendre legendre.f bve.inc

done

rm gausspoints.dat