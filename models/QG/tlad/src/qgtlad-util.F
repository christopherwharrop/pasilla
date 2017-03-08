CC NOTE - SWITCHED VARIABLE "GJACOB" TO "ININAG", AS THEY MATCH

      SUBROUTINE ADDAB(AA,BB)
C
C***  computes the sum of AA and BB, 
C***  the result is put in AA.
C
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
      REAL*8 AA(NSH2),BB(NSH2)
      INTEGER K
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
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
      REAL*8 AA(NSH2),BB(NSH2)
      INTEGER K
C
      DO 10 K=1,NSH2
        AA(K)=AA(K)-BB(K)
   10 CONTINUE
C
      RETURN
      END
      
      SUBROUTINE NEG(AA)
C *** Negate AA
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
      REAL*8 AA(NSH2)
      INTEGER K
      
      DO 10 K=1,NSH2
        AA(K)=-AA(K)
 10   CONTINUE
 
      RETURN
      END
      
      
      SUBROUTINE INPRSQ(XS,YS,E)
C Compute the squared norm inner product
      IMPLICIT NONE 
#include "truncation.h"
#include "comqg.h"
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

      SUBROUTINE ADD(AA,BB,CC)
C
C***  computes the sum of AA and BB, 
C***  the result is put in CC.
C
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2),BB(NSH2),CC(NSH2)
      INTEGER K
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2,3),BB(NSH2,3),CC(NSH2,3)
      INTEGER K,L
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2),BB(NSH2),CC(NSH2)
      INTEGER K
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2,3),BB(NSH2,3),CC(NSH2,3)
      INTEGER K,L
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2),BB(NSH2),C
      INTEGER K
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2,3),BB(NSH2,3),C
      INTEGER K,L
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2),BB(NSH2)
      INTEGER K
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2,3),BB(NSH2,3)
      INTEGER K,L
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2,3),BB(NSH2,3),E
      INTEGER K,L
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 AA(NSH2,3),BB(NSH2,3),E
      INTEGER K,L
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 XS(NSH2,3),E
      INTEGER K,L
C
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
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 XS(NSH2),E
      INTEGER K,L
C
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

      
      SUBROUTINE LAP3(XS,YS)
C
C***  computes the Laplacian of XS, the result is put in YS.
C
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 XS(NSH2,3),YS(NSH2,3)
      INTEGER K,L
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

      SUBROUTINE LAPINV3(XS,YS)
C
C***  computes the inverse Laplacian of XS
C***  the result is put in YS.
C
      IMPLICIT NONE      
#include "truncation.h"
#include "comqg.h"  
      REAL*8 XS(NSH2,3),YS(NSH2,3)
      INTEGER K,L
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

      SUBROUTINE RETRIE(IPLACE)
C RETRIEVE ELEMENT IPLACE FROM THE STORED ORBIT
C PUT THE DATA IN PSI
C DATE: 25-3-1991
C STATUS: TESTED
cc    IMPLICIT REAL*8 (A-H,O-Z)
cc    include 'comfra.f'
cc    include 'comtes.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      INTEGER I,J,IPLACE
      print *,"REAL RETRIE",IPLACE
      IF (IPLACE.LT.0) THEN
        WRITE(*,*)  'RETRIE: FATAL ERROR IPLACE LT 0',IPLACE
        STOP
      ELSE IF (IPLACE.GT.NTADJ+1) THEN
        WRITE(*,*)  'RETRIE: FATAL ERROR IPLACE GT NTADJ+1',
     &    IPLACE
      ENDIF
      DO 1,I=1,3
       DO 2,J=1,NSH2
        PSI(J,I)=PORBIT(J,I,IPLACE)
 2     CONTINUE
 1    CONTINUE
 
      RETURN
      END
      

      SUBROUTINE STOREP(IPLACE)
C STORE ELEMENT IPLACE
C PUT THE DATA FROM QPRIME INTO QORBIT
cc    IMPLICIT REAL*8 (A-H,O-Z)
cc    include 'comfra.f'
cc    include 'comtes.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"

      INTEGER I,J,IPLACE
c     print *,"REAL STOREP",IPLACE 
      IF (IPLACE.LT.0) THEN
        WRITE(*,*) 'STOREQ: FATAL ERROR IPLACE LT 0',
     &    IPLACE
        STOP
      ELSE IF (IPLACE.GT.NTADJ+1) THEN
        WRITE(*,*) 'STOREQ: FATAL ERROR IPLACE GT NTSTEP+1: ',
     &    IPLACE
        STOP
      ENDIF
      DO 1,I=1,3
       DO 2,J=1,NSH2
        PORBIT(J,I,IPLACE)=PSI(J,I)
 2     CONTINUE
 1    CONTINUE
 
cc    if((0.eq.1).and.(ltangwr)) then
cc       if(IPLACE.eq.0) then
cc         write(*,*) 'WRITE NL TRAJ to 51 !!!'
cc       endif
cc       if(mod(IPLACE,9).eq.0) then
cc         k=(IPLACE/3)*2
cc         call tometview(51,psi,k,138)
cc       endif
cc    endif
      
      RETURN
      END


      SUBROUTINE JACOBP (PSILOC,PVOR,SJACOB)
C
C *** Advection of potential vorticity and dissipation on gaussian grid
C *** Adapted version of JACOBD for testing the tangent linear equations
cc    IMPLICIT REAL*8 (A-H,O-Z)
cc    COMMON /CLINEA/LITEST(10)
cc    LOGICAL LITEST
cc    include 'comfra.f'
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
#include "comqgtlad.h"
      REAL*8 PSILOC(NSH2), PVOR(NSH2), SJACOB(NSH2)
      INTEGER I,J,K,L,M,N
C
C *** Isotropic dissipation term
C
cc    IF (LITEST(10)) THEN
        DO 50 K=1,NSH2
   50   VV(K)=DISS(K,2)*PSILOC(K)
C
        CALL SPTOGG (VV,ININAG,PP)
cc    ELSE
cc      DO 150 K=1,NGP
cc150     ININAG(K)=0.0D0        
cc    ENDIF
C
C *** Space derivatives of potential vorticity
C
cc    IF (LITEST(7).OR.LITEST(8)) THEN
        CALL DDL (PVOR,VV)
        CALL SPTOGG (VV,DVORDL,PP)
        CALL SPTOGG (PVOR,DVORDM,PD)
cc    ENDIF
C
C *** Space derivatives of streamfunction
C
      CALL DDL (PSILOC,VV)
      CALL SPTOGG (VV,DPSIDL,PP)
      CALL SPTOGG (PSILOC,DPSIDM,PD)
C
C *** Jacobian and unisotropic dissipation terms
C
      
      DO 100 I=1,NLAT 
        DO 100 J=1,NLON
cc      IF (LITEST(10)) THEN
          ININAG(I,J)=-DPSIDM(I,J)*DDISDY(I,J)
     *                -DPSIDL(I,J)*DDISDX(I,J)
     *                +RDISS(I,J)*ININAG(I,J)
cc      ENDIF
cc      IF (LITEST(7).OR.LITEST(8)) THEN
          ININAG(I,J)=DPSIDM(I,J)*DVORDL(I,J)-
     *                DPSIDL(I,J)*DVORDM(I,J)+
     *                ININAG(I,J)
cc      ENDIF
  100 CONTINUE
C
      CALL GGTOSP (ININAG,SJACOB)
C
C *** Planetary vorticity advection
C
cc    IF (LITEST(7).OR.LITEST(8)) THEN
        DO 200 K=1,NSH2
  200   SJACOB(K)=SJACOB(K)-VV(K)
cc    ENDIF
C
      RETURN
      END

      SUBROUTINE QFORTOPSI(FORQ,FORP)
      IMPLICIT NONE
#include "truncation.h"
#include "comqg.h"
      integer k
      real*8  r3
      real*8  forp(nsh2,nvl),forq(nsh2,nvl)
      real*8  work(nsh2),    fort(nsh2,ntl)

      do k=1,nsh2
        work(k)=forq(k,1)+forq(k,3)
        forp(k,1)=rinhel(k,1)*(work(k)+forq(k,2))
        forp(k,2)=work(k)-2.d0*forq(k,2)
        forp(k,3)=forq(k,1)-forq(k,3)
      enddo

      do k=1,nsh2
        fort(k,1)=rinhel(k,2)*forp(k,2)+rinhel(k,3)*forp(k,3)
        fort(k,2)=rinhel(k,4)*forp(k,2)+rinhel(k,5)*forp(k,3)
      enddo

      r3=1./3.
      do k=1,nsh2
        forp(k,2)=r3*(forp(k,1)-fort(k,1)+fort(k,2))
        forp(k,1)=forp(k,2)+fort(k,1)
        forp(k,3)=forp(k,2)-fort(k,2)
      enddo

      RETURN
      END


      SUBROUTINE ALOMA3SP(EPSIN,EPSOUT)
C adjoint of tmalo3sp wrt KE innerproduct

      IMPLICIT NONE
cc    include 'comfra.f'
cc    include 'comtes.f'
#include "truncation.h"
#include "comqg.h"
      real*8 epsin(nsh2,3),epsout(nsh2,3),xs(nsh2,3),
     &       agg(nlat,nlon)
      integer i,j,k,l,ilat,ilon

      do l=1,3
        call lap(epsin(1,l),xs(1,l))
        call sptogg(xs(1,l),agg,pp)
cc      do ilat=1,nlat
cc        do ilon=1,nlon
cc          agg(ilat,ilon)=agg(ilat,ilon)*zmask(ilat,ilon)
cc        enddo
cc      enddo
        call ggtosp(agg,epsout(1,l))
        call lapinv(epsout(1,l),epsout(1,l))
      enddo

      return
      end

      SUBROUTINE TMALO3SP(EPSIN,EPSOUT)
C localization of errorfield EPSIN 
C masking array ZMASK must be precomputed.
      IMPLICIT NONE
cc    include 'comfra.f'
cc    include 'comtes.f'
#include "truncation.h"
#include "comqg.h"
      real*8 epsin(nsh2,3),epsout(nsh2,3),
     &       agg(nlat,nlon),e,fac
      integer l,ilat,ilon

      do l=1,3
        call sptogg(epsin(1,l),agg,pp)
cc      do ilat=1,nlat
cc        do ilon=1,nlon
cc          agg(ilat,ilon)=agg(ilat,ilon)*zmask(ilat,ilon)
cc        enddo
cc      enddo
        call ggtosp(agg,epsout(1,l))
      enddo

      return
      end
 
