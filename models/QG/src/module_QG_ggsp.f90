module QG_GGSP

  use kind

  implicit none

  private

  public :: qg_ggsp_type

  type :: qg_ggsp_type
    private
    integer :: nm                            ! The truncation is of type T(riangular) nm
    integer :: nlon                          ! Number of longitude points of the Gaussian grid
    integer :: nlat                          ! Number of latitude  points of the Gaussian grid
    integer :: nsh
    integer, allocatable :: nshm(:)          ! Contains numbers 22 down to 1 for index 0 to 21
    real(r8kind), allocatable  :: pp(:,:)    ! Legendre polynomials defined at Gausian latitudes
    real(r8kind), allocatable  :: pd(:,:)    ! Mu derivative of Legendre polynomials
    real(r8kind), allocatable  :: pw(:,:)    ! Weights for Legendre integrals
    real(r8kind), allocatable  :: trigd(:,:) ! Trigonometric coefficients used by the nag version of the fft
    real(r8kind), allocatable  :: trigi(:,:) ! Trigonometric coefficients used by the nag version of the fft
    real(r8kind), allocatable  :: rm(:)      ! contains zonal wavenumber m of each spherical harmonic of the corresponding index for zonal derivative operator
  contains
    final :: destructor_qg_ggsp
    procedure          :: sptogg_pp
    procedure          :: sptogg_pd
    procedure          :: sptogg_pw
    procedure, private :: sptogg
    procedure          :: ggtosp
    procedure          :: ddl
    procedure          :: ddl_b
  end type qg_ggsp_type

  interface qg_ggsp_type
    procedure constructor_qg_ggsp
  end interface


contains


  !-------------------------------------------------------------------------------
  ! constructor_qg_ggsp
  !-------------------------------------------------------------------------------
  function constructor_qg_ggsp(nm, nlat, nlon, nshm, pp, pd, pw) result (qg_ggsp)

    ! Input
    integer, intent(in) :: nm            ! The truncation is of type T(riangular) nm
    integer, intent(in) :: nlon          ! Number of longitude points of the Gaussian grid
    integer, intent(in) :: nlat          ! Number of latitude  points of the Gaussian grid
    integer, intent(in) :: nshm(0:nm)    ! Contains numbers 22 down to 1 for index 0 to 21
    real(r8kind), intent(in) :: pp(:,:)  ! Legendre polynomials defined at Gausian latitudes
    real(r8kind), intent(in) :: pd(:,:)  ! Mu derivative of Legendre polynomials
    real(r8kind), intent(in) :: pw(:,:)  ! Weights for Legendre integrals

    ! Return value
    type(qg_ggsp_type) :: qg_ggsp

    ! Local data
    integer :: i, j, k, k1, k2, m, ifail
    real(r8kind), allocatable :: ininag(:,:)  ! FFT initialization field
    real(r8kind), allocatable :: tmp(:,:)     ! Work space used by the nag version of the FFT

    qg_ggsp%nm = nm
    qg_ggsp%nlon = nlon
    qg_ggsp%nlat = nlat
    qg_ggsp%nsh = ((nm + 1) * (nm + 2)) / 2

    allocate(qg_ggsp%nshm(0:nm))
    qg_ggsp%nshm(:) = nshm(:)

    allocate(qg_ggsp%pp(nlat,qg_ggsp%nsh))
    qg_ggsp%pp(:,:) = pp(:,:)

    allocate(qg_ggsp%pd(nlat,qg_ggsp%nsh))
    qg_ggsp%pd(:,:) = pd(:,:)

    allocate(qg_ggsp%pw(nlat,qg_ggsp%nsh))
    qg_ggsp%pw(:,:) = pw(:,:)

    ! Initialization of trigonometric coefficients for fft
    allocate(ininag(nlat, nlon))
    allocate(qg_ggsp%trigd(nlon,2), qg_ggsp%trigi(nlon,2))
    allocate(tmp(nlat,nlon))
    do j = 1, nlon
      do i = 1, nlat
        ininag(i, j) = 1.0d0
      enddo
    enddo 
    ifail = 0
    call c06fpf(nlat, nlon, ininag, 'i', qg_ggsp%trigd, tmp, ifail)
    ifail = 0
    call c06fqf(nlat, nlon, ininag, 'i', qg_ggsp%trigi, tmp, ifail)


    ! Precalculate zonal derivative operator
    allocate(qg_ggsp%rm(qg_ggsp%nsh))
    k2 = 0
    do m = 0, nm
      k1 = k2 + 1
      k2 = k2 + nshm(m)
      do k = k1, k2
        qg_ggsp%rm(k) = dble(m)
      enddo
    enddo


  end function constructor_qg_ggsp


  !------------------------------------------------------------------
  ! destructor_qg_ggsp
  !
  ! Deallocates pointers used by a qg_ggsp_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_qg_ggsp(this)

    type(qg_ggsp_type), intent(inout) :: this

    ! No pointers in qg_model_type object so we do nothing

  end subroutine destructor_qg_ggsp


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid using
  ! legendre polynomials
  !
  ! input  spectral field as
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg_pp (this, as) result(agg)

    ! Input
    class(qg_ggsp_type), intent(inout) :: this
    real(r8kind),        intent(   in) :: as(this%nsh, 2)

    ! Return value
    real(r8kind) :: agg(this%nlat, this%nlon)

    agg = this%sptogg(as, this%pp)

  end function sptogg_pp


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid using
  ! derivatives with respect to sin(fi) .
  !
  ! input  spectral field as
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg_pd (this, as) result(agg)

    ! Input
    class(qg_ggsp_type), intent(inout) :: this
    real(r8kind),        intent(   in) :: as(this%nsh, 2)

    ! Return value
    real(r8kind) :: agg(this%nlat, this%nlon)

    agg = this%sptogg(as, this%pd)

  end function sptogg_pd


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid using
  ! derivatives with respect to pw
  !
  ! input  spectral field as
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg_pw (this, as) result(agg)

    ! Input
    class(qg_ggsp_type), intent(inout) :: this
    real(r8kind),        intent(   in) :: as(this%nsh, 2)

    ! Return value
    real(r8kind) :: agg(this%nlat, this%nlon)

    agg = this%sptogg(as, this%pw)

  end function sptogg_pw


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid
  ! input  spectral field as,  legendre polynomials pploc (pp or pd) 
  !        where pp are legendre polynomials and pd derivatives with
  !        respect to sin(fi)
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg (this, as, pploc) result(agg)

    ! Input 
    class(qg_ggsp_type), intent(inout) :: this
    real(r8kind),        intent(   in) :: as(:,:)
    real(r8kind),        intent(   in) :: pploc(:,:)

    ! Return value
    real(r8kind) :: agg(this%nlat, this%nlon)

    ! Local data
    real(r8kind), allocatable :: tmp(:,:)   ! Work space used by the nag version of the fft    
    integer :: i, j, k, k1, k2, m, mi, mr, nlon1
    integer :: ifail

    ! inverse legendre transform
    do j = 1, this%nlon
      do i = 1, this%nlat
        agg(i, j) = 0.0d0
      enddo
    enddo

    nlon1 = this%nlon + 1
    k2 = this%nshm(0)

    do k = 1, k2
      do i = 1, this%nlat
        agg(i, 1) = agg(i, 1) + as(k, 1) * pploc(i, k)
      enddo
    enddo

    do m = 1, this%nm
      mr = m + 1
      mi = nlon1 - m
      k1 = k2 + 1
      k2 = k2 + this%nshm(m)
      do k = k1, k2
        do i = 1, this%nlat
          agg(i, mr) = agg(i, mr) + as(k, 1) * pploc(i, k)
        enddo
        do i = 1, this%nlat
          agg(i, mi) = agg(i, mi) - as(k, 2) * pploc(i, k)
        enddo
      enddo
    enddo

    ! inverse fourier transform
    ifail = 0
    allocate(tmp(this%nlat,this%nlon))
    call c06fqf (this%nlat, this%nlon, agg, 'r', this%trigi, tmp, ifail)

  end function sptogg
 

  !-----------------------------------------------------------------------
  ! conversion from gaussian grid (agg) to spectral coefficients (as)
  ! input gaussian grid field agg
  ! output as contains spectral coefficients
  !-----------------------------------------------------------------------
  function ggtosp (this, agg) result (as)

    ! Input
    class(qg_ggsp_type), intent(inout) :: this
    real(r8kind),        intent(   in) :: agg(:,:)

    ! Return value
    real(r8kind) :: as(this%nsh, 2)

    ! Local data
    real(r8kind), allocatable :: agg_copy(:,:)  ! Copy of input gaussian grid field
    real(r8kind), allocatable :: tmp(:,:)       ! Work space used by the nag version of the fft    
    integer :: i, k, k1, k2, m, mi, mr, nlon1
    integer :: ifail

    ! Make a local copy of agg so it is not destroyed by c06fpf
    allocate(agg_copy(this%nlat,this%nlon))
    agg_copy(:,:) = agg(:,:)

    ! fourier transform
    ifail = 0
    allocate(tmp(this%nlat,this%nlon))
    call c06fpf (this%nlat, this%nlon, agg_copy, 'r', this%trigd, tmp, ifail)

    ! legendre transform
    do i = 1, 2
      do k = 1, this%nsh
        as(k, i) = 0.0d0
      enddo
    enddo

    nlon1 = this%nlon + 1

    k2 = this%nshm(0)

    do k = 1, k2
      do i = 1, this%nlat
        as(k, 1) = as(k, 1) + agg_copy(i, 1) * this%pw(i, k)
      enddo
    enddo

    do m = 1, this%nm
      mr = m + 1
      mi = nlon1 - m
      k1 = k2 + 1
      k2 = k2 + this%nshm(m)
      do k = k1, k2
        do i = 1, this%nlat
          as(k, 1) = as(k, 1) + agg_copy(i, mr) * this%pw(i, k)
          as(k, 2) = as(k, 2) + agg_copy(i, mi) * this%pw(i, k)
        enddo
      enddo
    enddo

  end function ggtosp


  !-----------------------------------------------------------------------
  ! zonal derivative in spectral space
  ! input spectral field as
  ! output spectral field dadl which is as differentiated wrt lambda
  !-----------------------------------------------------------------------
  pure function ddl (this,as) result(dadl)

    implicit none

    class(qg_ggsp_type), intent(in) :: this
    real(r8kind),        intent(in) :: as(this%nsh, 2)
    real(r8kind)                    :: dadl(this%nsh, 2)

    integer :: k

    do k = 1, this%nsh
      dadl(k, 1) = -this%rm(k) * as(k, 2)
      dadl(k, 2) =  this%rm(k) * as(k, 1)
    enddo

    return

  end function ddl


  !-----------------------------------------------------------------------
  !  Differentiation of ddl in reverse (adjoint) mode:
  !   gradient     of useful results: as dadl
  !   with respect to varying inputs: as
  !-----------------------------------------------------------------------
  ! zonal derivative in spectral space
  ! input spectral field as
  ! output spectral field dadl which is as differentiated wrt lambda
  !-----------------------------------------------------------------------
  subroutine ddl_b(this, as, asb, dadlb)

    class(qg_ggsp_type), intent( in) :: this
    real(r8kind),        intent( in) :: as(this%nsh, 2)
    real(r8kind),        intent(out) :: asb(this%nsh, 2)
    real(r8kind),        intent(out) :: dadlb(this%nsh, 2)

    integer :: k

    do k=this%nsh,1,-1
      asb(k, 1) = asb(k, 1) + this%rm(k) * dadlb(k, 2)
      dadlb(k, 2) = 0.0_r8kind
      asb(k, 2) = asb(k, 2) - this%rm(k) * dadlb(k, 1)
      dadlb(k, 1) = 0.0_r8kind
    end do

  end subroutine ddl_b


end module QG_GGSP
