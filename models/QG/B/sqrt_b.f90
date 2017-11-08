      program sqrtb

      use netcdf

      ! General variables
      implicit none
      include "mkl.fi"
      ! General B-matrix variables
      integer                      ::      bkg_len
      integer                      ::      LWORK
      real(KIND=4), allocatable    ::      bkg_cov(:,:)
      ! General LAPACK variables
      integer                      ::      I,N,M
      real(KIND=4), allocatable    ::      A(:,:),W(:),Z(:,:),Y(:,:)
      real(KIND=4), allocatable    ::      WORK(:)
      ! General netCDF variables
      integer :: ncFileID  ! netCDF file identifier
      integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
      integer :: BVarDimID, CoordinatesVarID, LocationVarID, BVarID
      integer :: info, ierr

      bkg_len=64*32*3
!      bkg_len=128*64*3
!     bkg_len=3
      N=bkg_len
      allocate (bkg_cov(bkg_len,bkg_len))
      allocate (W(bkg_len))
      allocate (Z(bkg_len,bkg_len))
      allocate (Y(bkg_len,bkg_len))
      allocate (A(bkg_len,bkg_len))
!     print *,"GET_BKG_COV_MAT"

      ! Open file for read only
      call nc_check(nf90_open("b_sqr.nc", NF90_NOWRITE, ncFileID))
      call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables,   & 
      nAttributes, unlimitedDimID))
      ! Get the B variable ID
      call nc_check(nf90_inq_varid(ncFileID, "B", BVarID))
      ! Get the B variable
      call nc_check(nf90_get_var(ncFileID, BVarID, bkg_cov(:,:)))
      ! Flush buffers
      call nc_check(nf90_sync(ncFileID))
      ! Close the NetCDF file
      call nc_check(nf90_close(ncFileID))
!     bkg_cov(1,1)=3.0 
!     bkg_cov(2,1)=2.0 
!     bkg_cov(3,1)=1.0 
!     bkg_cov(1,2)=2.0
!     bkg_cov(2,2)=2.0
!     bkg_cov(3,2)=2.0
!     bkg_cov(1,3)=1.0
!     bkg_cov(2,3)=2.0
!     bkg_cov(3,3)=3.0
!     write (*,'(A10,9F6.2)')  "B BEFORE = ",bkg_cov
      bkg_cov(:,:)=36.0*bkg_cov(:,:)
      print *,"GET_BKG_COV_MAT COMPLETE"


      print *,"SQRT OF BKG_COV"
      LWORK = -1
      allocate(WORK(1))
      CALL ssyev('V','U',N,bkg_cov,N,W,WORK,LWORK,INFO)
      LWORK = INT(WORK(1))
      print *,"LWORK = ",LWORK
      deallocate(WORK)
      allocate (WORK(LWORK))
      call ssyev('V','U',N,bkg_cov,N,W,WORK,LWORK,INFO)
      deallocate(WORK)
      print *, "EIG-VALS: ",W(1:3) 

      A(:,:)=0.0
      M=0
      do 100 i=1,bkg_len
        if(W(i).gt.0.0) M=M+1
!       if(W(i).gt.0.0) A(i,i)=sqrt(W(i)) 
!       if(W(i).gt.0.0) A(i,i)=(W(i))
        if(W(i).gt.0.0) A(i,i)=sqrt(sqrt(W(i)))
        if(W(i).le.0.0) A(i,i)=sqrt(sqrt(-W(i)))
100   continue
      call sgemm("N","N",N,N,N,1.0,bkg_cov,N,A,N,0.0,Z,N)

      A(:,:)=0.0
      do 110 i=1,bkg_len
!       if(W(i).gt.0.0) A(i,i)=sqrt(W(i))
!       if(W(i).gt.0.0) A(i,i)=(W(i))
        if(W(i).gt.0.0) A(i,i)=sqrt(sqrt(W(i)))
        if(W(i).le.0.0) A(i,i)=-sqrt(sqrt(-W(i)))
110   continue
      call sgemm("N","N",N,N,N,1.0,Z,N,A,N,0.0,Y,N)  
      call sgemm("N","T",N,N,N,1.0,Y,N,bkg_cov,N,0.0,A,N)
      bkg_cov(:,:)=A(:,:)

      Do 200 m=1,bkg_len
        Do 200 n=1,bkg_len
           if (bkg_cov(m,n).lt.1.0) bkg_cov(m,n)=0.0
200   continue

!     write (*,'(A10,9F6.2)') "B AFTER = ",bkg_cov 
      print *,"SQRT OF BKG_COV COMPLETE"


      print *,"PUT_SQRT_BKG_COV"
      ! assume normal termination
      ierr = 0 
      ! Open new file, overwriting previous contents
      call nc_check(nf90_create("b.nc", OR(NF90_CLOBBER,         &
        NF90_64BIT_OFFSET), ncFileID))
      call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables,   &
        nAttributes, unlimitedDimID))
      ! Define the matrix size
      call nc_check(nf90_def_dim(ncid=ncFileID, name="BDim",          &
        len=bkg_len, dimid = BVarDimID))
      ! Define the actual Matrix data 
      call nc_check(nf90_def_var(ncid=ncFileID, name="B",             &
        xtype=nf90_real,dimids=(/BVarDimID, BVarDimID/), varid=BVarID))
      call nc_check(nf90_put_att(ncFileID, BVarID, "long_name", "B    &
        Covariance Matrix"))
      call nc_check(nf90_put_att(ncFileID, BVarID, "units",           &
        "Nondimensional"))
      ! Leave define mode so we can fill
      call nc_check(nf90_enddef(ncfileID))
      ! Fill the b variable
      call nc_check(nf90_put_var(ncFileID, BVarID, bkg_cov))
      ! Flush buffers
      call nc_check(nf90_sync(ncFileID))
      ! Close the NetCDF file
      call nc_check(nf90_close(ncFileID))

      print *,"PUT_SQRT_BKG_COV COMPLETE"



      end

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

