c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comqgtlad.h                                                    
c *** Contents: Common declarations for three level qg model tl and ad
c-----------------------------------------------------------------------
c *** COMMON /qphist/ PORBIT
c     porbit: the trajectories of streamfunction, for use in the TL/AD
c
c *** COMMON /backin/ DADOLD,DADNEW,DDADNE
c     dadold:
c     dadnew:
c     ddadne:
c
c *** COMMON /backf/  FOROLD,FORNEW,DFORNE
c     forold:
c     fornew:
c     dforne:
c
c *** COMMON /gfield/ DPSIDL,DPSIDM,DVORDL,DVORDM,ININAG,VV
c     DPSIDL
c     DPSIDM
c     DVORDL
c     DVORDM
c     ININAG
c     VV
c
c *** COMMON /etraj/  ETRAJ
c     etraj:
c-----------------------------------------------------------------------

      REAL*8 PORBIT(NSH2,3,0:NTADJ+1)
      REAL*8 DADOLD(NSH2,3),DADNEW(NSH2,3),DDADNE(NSH2,3)
      REAL*8 FOROLD(NSH2,3),FORNEW(NSH2,3),DFORNE(NSH2,3)
      COMMON /QPHIST/ PORBIT
      COMMON /BACKIN/ DADOLD,DADNEW,DDADNE
      COMMON /BACKF/ FOROLD,FORNEW,DFORNE
      REAL*8 DPSIDL,DPSIDM,DVORDL,DVORDM,ININAG,VV
      COMMON /GFIELD/ DPSIDL(NLAT,NLON), DPSIDM(NLAT,NLON),
     *                DVORDL(NLAT,NLON), DVORDM(NLAT,NLON),
     *                ININAG(NLAT,NLON), VV(NSH2)
      real*8 etraj
      common etraj(nsh2,3,0:100)
