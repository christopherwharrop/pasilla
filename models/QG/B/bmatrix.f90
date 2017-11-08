      PROGRAM BMATRIX

      use netcdf

      IMPLICIT NONE
      INTEGER X,Y,Z
      PARAMETER(X=64, Y=32,Z=3)
      REAL LAT(X,Y)
      REAL LON(X,Y)
      REAL BIN(X,Y)
!      REAL BBB(X,Y,Z)
      REAL, allocatable :: BBB(:,:,:)
!      REAL FBBB(X,Y,Z,X*Y*Z)
!      REAL FBBB(X * Y * Z, X * Y * Z)
      REAL, allocatable :: FBBB(:,:)
      INTEGER XX,YY,ZZ,XL,YL
      INTEGER*8 MASKS 
      REAL A,C,DLAT,DLON,TLAT,TLAX,DTOR,ERAD,TMPA,TMPC,DIST,LFAC
      REAL PARM1, PARM2
      REAL MASK(X,Y)
      INTEGER*8 nrow
      REAL ratio
      INTEGER SUMM(3,3),GOOD(3,3),ZERO(3,3),SYMM(3,3)
    integer :: ierr          ! return value of function

    ! General netCDF variables
    integer :: ncFileID      ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: BVarDimID, BVarID


    allocate(bbb(x,y,z))
    allocate(fbbb(X * Y * Z, X * Y * Z))

      DTOR=ATAN(1.0)/45.0
      ERAD=6372.8
      PARM1=3.0*111.1
      PARM2=4.0*PARM1
      MASKS=0

      PRINT *,"SIZE OF OUTPUT: ",SIZE(FBBB)
      OPEN (20,FILE="latlon.dat",FORM="UNFORMATTED")
      READ (20) LAT
      PRINT *,"LAT",LAT(1,1),LAT(64,32)
      READ (20) LON
      PRINT *,"LON",LON(1,1),LON(64,32)
      CLOSE (20)

      OPEN (40,FILE="bmatrix.rows.dat",FORM="UNFORMATTED")

      nrow=1
      DO  ZZ=1,Z

        IF(ZZ.EQ.1) THEN
          OPEN (31,FILE="psidiff.z1.dat",FORM="UNFORMATTED")
        ENDIF
        IF(ZZ.EQ.2) THEN
          OPEN (31,FILE="psidiff.z2.dat",FORM="UNFORMATTED")
        ENDIF
        IF(ZZ.EQ.3) THEN
          OPEN (31,FILE="psidiff.z3.dat",FORM="UNFORMATTED")
        ENDIF

        DO YY=1,Y
          DO XX=1,X

            DO YL=1,Y
              DO XL=1,X
                DLAT=DTOR*(LAT(XL,YL)-LAT(XX,YY)) 
                DLON=DTOR*(LON(XL,YL)-LON(XX,YY)) 
                TLAT=DTOR*LAT(XL,YL)
                TLAX=DTOR*LAT(XX,YY)
                A=(SIN(DLAT/2.0))**2 + COS(TLAT)*COS(TLAX)*(SIN(DLON/2.0))**2
                C=2.0*ASIN(SQRT(A))
                DIST=ERAD*C
!               LFAC=0.5+0.5*ABS(COS(TLAX))
                LFAC=1.0
                MASK(XL,YL)=1.0
          !     PRINT *,XX,YY,XL,YL
          !     PRINT *,LAT(XX,YY),LON(XX,YY),LAT(XL,YL),LON(XL,YL)
          !     PRINT *,DLAT,DLON
          !     PRINT *,SIN(DLAT/2.0),COS(TLAT),COS(TLAX),SIN(DLON/2.0)
          !     PRINT *,A,C,DIST,LFAC
                IF(DIST.GT.(PARM1/LFAC)) MASK(XL,YL)=EXP((LFAC/PARM1)*(PARM1/LFAC-DIST/LFAC))
                IF(DIST.GT.(PARM2/LFAC)) MASK(XL,YL)=0.0
          !     PRINT *,"MASK: ",MASK(XL,YL)
          !     PRINT *,"-------------------------"
!               MASK(XL,YL)=0.0
!               IF((XX.EQ.XL).AND.(YY.EQ.YL)) MASK(XL,YL)=1.0
              end do
            end do

            READ (31) BIN
!           BIN(:,:)=SQRT(ABS(BIN(:,:)))*MASK(:,:) 
            BIN(:,:)=BIN(:,:)*MASK(:,:)
            IF(ZZ.NE.1) BIN(:,:)=0.0
            BBB(:,:,1)=BIN(:,:) 

            READ (31) BIN
!           BIN(:,:)=SQRT(ABS(BIN(:,:)))*MASK(:,:) 
            BIN(:,:)=BIN(:,:)*MASK(:,:)
            IF(ZZ.NE.2) BIN(:,:)=0.0          
            BBB(:,:,2)=BIN(:,:)

            READ (31) BIN
!           BIN(:,:)=SQRT(ABS(BIN(:,:)))*MASK(:,:) 
            BIN(:,:)=BIN(:,:)*MASK(:,:)
            IF(ZZ.NE.3) BIN(:,:)=0.0  
            BBB(:,:,3)=BIN(:,:)

!            WRITE (40) BBB
            FBBB(nrow,:) = reshape(BBB,(/x * y * z/))
            MASKS=MASKS+1
      !     PRINT *,XX,YY,ZZ,MASKS

            nrow = nrow + 1

          end do
        end do

        CLOSE(31)

      end do

      WRITE (40) FBBB
      CLOSE(40)

    ! assume normal termination
    ierr = 0 

    ! Open new file, overwriting previous contents
    call nc_check(nf90_create("b_sqr.nc", OR(NF90_CLOBBER, NF90_64BIT_OFFSET), ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Define the matrix size
    call nc_check(nf90_def_dim(ncid=ncFileID, name="BDim", &
                               len=x*y*z, dimid = BVarDimID))

    ! Define the actual Matrix data 
    call nc_check(nf90_def_var(ncid=ncFileID, name="B", xtype=nf90_real, &
               dimids=(/BVarDimID, BVarDimID/), varid=BVarID))
    call nc_check(nf90_put_att(ncFileID, BVarID, "long_name", "B Covariance Matrix"))
    call nc_check(nf90_put_att(ncFileID, BVarID, "units", "Nondimensional"))

    ! Leave define mode so we can fill
    call nc_check(nf90_enddef(ncfileID))

    ! Fill the b variable
    call nc_check(nf90_put_var(ncFileID, BVarID, FBBB))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    PRINT *,"TOTAL NUMBER OF MASKS: ",MASKS

    PRINT *,"SIZE!",SIZE(FBBB)
     MASKS=0
     SYMM(:,:)=0
     GOOD(:,:)=0
     ZERO(:,:)=0
     SUMM(:,:)=0.0
     Do 100 XX=1, (X*Y*Z)
!    print *, "STEAD AT: ",XX 
     MASKS=MASKS+X*Y*Z
     Do 100 YY=1, (X*Y*Z)
        XL=1+((XX-1)/(X*Y))
        YL=1+((YY-1)/(X*Y))
        SUMM(XL,YL)=SUMM(XL,YL)+1
        ratio=(abs(fbbb(XX,YY))+1.0)/(abs(fbbb(YY,XX))+1.0)
        if ((ratio.gt.1.01).or.(ratio.lt.0.99)) then 
            if ((fbbb(XX,YY).gt.0.1).and.(fbbb(XX,YY).gt.0.1)) then 
!           print *,XX,YY,fbbb(xx,yy),fbbb(yy,xx),ratio
            symm(XL,YL)=symm(XL,YL)+1
            else
            zero(XL,YL)=zero(XL,YL)+1
            endif
        else
        good(XL,YL)=good(XL,YL)+1
        endif   
100   continue

      Do 101 XX=1,3
      Do 101 YY=1,3
        print *,XX,YY,SUMM(XX,YY),GOOD(XX,YY),SYMM(XX,YY),ZERO(XX,YY)
101   continue
      PRINT *,"TOTAL NUMBER OF MASKS, TWO: ",MASKS
!     PRINT *,"TOTAL NUMBER OF VALUES: ",MASKS*X*Y*Z
!     PRINT *,"EXPECTED SIZE OF FILE: ",(MASKS*X*Y*Z*4)+MASKS*8
      END
      
  !------------------------------------------------------------------
  ! nc_check
  ! 
  ! Checks return status from a NetCDF API call.  If an error was
  ! returned, print the message and abort the program.
  !------------------------------------------------------------------
  subroutine nc_check(istatus)

    use netcdf

    integer, intent (in)                   :: istatus
  
    character(len=512) :: error_msg
  
    ! if no error, nothing to do here.  we are done.
    if( istatus == nf90_noerr) return

    error_msg = nf90_strerror(istatus)
  
    print *,error_msg
    stop  

  end subroutine nc_check
