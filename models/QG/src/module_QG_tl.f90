module QG_Model_TL

  use kind
  use Abstract_Model, only : abstract_model_type
  use QG_Config
  use QG_GGSP
  use QG_Util

  implicit none

  private

  public :: qg_tl_type

  type, extends(abstract_model_type) :: qg_tl_type
    private
    ! Model configuration
    type(qg_config_type) :: config

    ! Model grid dimensions
    character(len=3) ::  ft   ! Character string containing resolution
    integer :: nm             ! The truncation is of type T(riangular) nm
    integer :: nlon           ! Number of longitude points of the Gaussian grid
    integer :: nlat           ! Number of latitude  points of the Gaussian grid
    integer :: nvl            ! Number of vorticity levels in the vertical (should be set to 3)
    integer :: ntl            ! Number of temperature levels in the vertical (equal to nvl-1)
    integer :: nsh            ! Half of nsh2
    integer :: nsh2           ! Number of coefficients needed to define one level of the T nm model
    integer :: nsh23          ! Number of coefficients needed to define three levels of the T nm model

    ! Current model step and simulaiton clock
    integer      :: step
    real(r8kind) :: clock

    ! Model time step
    real(r8kind) :: dtt                      ! dimensionless time step

    ! Model state
    real(r8kind), allocatable :: psi(:,:)        ! Stream function at the nvl levels
    real(r8kind), allocatable :: trajectory(:,:) ! Tangent Linear trajectory
    real(r8kind), allocatable :: psit(:,:)       ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprime(:,:)     ! Potential vorticity
    real(r8kind), allocatable :: psitd(:,:)       ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprimed(:,:)     ! Potential vorticity
    
    ! Model Forcing
    real(r8kind), allocatable :: for(:,:)    ! Constant potential vorticity forcing at the nvl levels

    ! Spectral Coefficients
    integer,      allocatable :: nshm(:)     ! Contains numbers 22 down to 1 for index 0 to 21
    integer,      allocatable :: ll(:)       ! Contains total wavenumber n of each spherical harmonic of the corresponding index
    real(r8kind), allocatable :: pp(:,:)     ! Legendre polynomials defined at Gausian latitudes
    real(r8kind), allocatable :: pd(:,:)     ! Mu derivative of Legendre polynomials
    real(r8kind), allocatable :: pw(:,:)     ! Weights for Legendre integrals
!    real(r8kind), allocatable :: trigd(:,:)  ! Trigonometric coefficients used by the nag version of the fft
!    real(r8kind), allocatable :: trigi(:,:)  ! Trigonometric coefficients used by the nag version of the fft
!    real(r8kind), allocatable :: rm(:)       ! contains zonal wavenumber m of each spherical harmonic of the corresponding index for zonal derivative operator

    ! Grid conversion object
    type(qg_ggsp_type) :: ggsp

    ! Laplace/Helmholtz direct and inverse operators
    real(r8kind), allocatable :: rinhel(:,:) ! Laplace and Helmholtz operator for Q-PSI inversion
    real(r8kind), allocatable :: diss(:,:)   ! Dissipation coefficients for each spherical harmonic
                                             !   diss(k,1) : Hyperviscosity at the three levels (tdif sets timescale)
                                             !   diss(k,2) : Ekman friction at lower level (tdis sets timescale)
    real(r8kind)              :: rl1         ! One over Rossby rad. of def. squared of 200-500 thickness
    real(r8kind)              :: rl2         ! One over Rossby rad. of def. squared of 500-800 thickness
    real(r8kind)              :: relt1       ! Nondimensional relaxation coefficient of 200-500 thickness
    real(r8kind)              :: relt2       ! Nondimensional relaxation coefficient of 500-800 thickness

    ! Orography
    real(r8kind), allocatable  :: phi(:)      ! Gauss points in radians
    real(r8kind), allocatable  :: sinfi(:)    ! Sine of phi
    real(r8kind), allocatable  :: cosfi(:)    ! Cosine of phi
    logical                    :: lgdiss      ! If .true. then orography and land-sea mask dependent friction at the lower level plus Ekman friction, else only Ekman friction
    real(r8kind), allocatable  :: orog(:)     ! Orography
    real(r8kind), allocatable  :: dorodl(:,:) ! Derivative of orog wrt lambda
    real(r8kind), allocatable  :: dorodm(:,:) ! Derivative of orag wrt sin(fi)
    real(r8kind), allocatable  :: rdiss(:,:)  ! Landsea-mask/orography dependent friction
    real(r8kind), allocatable  :: ddisdx(:,:) ! Landsea-mask/orography dependent friction
    real(r8kind), allocatable  :: ddisdy(:,:) ! Landsea-mask/orography dependent friction

  contains
    final :: destructor_qg_tl
    procedure :: adv_nsteps
    procedure :: tang
    procedure :: tang_reverse
    procedure :: gridfields
    procedure :: get_config
    procedure :: get_step
    procedure :: get_clock
    procedure :: get_nlat
    procedure :: get_nlon
    procedure :: get_nsh2
    procedure :: get_nvl
    procedure :: get_psi
    procedure :: get_trajectory
    procedure :: get_qprimed
    procedure :: get_psig
    procedure :: get_lat_lon_grid
    procedure :: get_for
    procedure :: get_state_vector
    procedure :: get_location_vector
    procedure :: get_interpolation_weights
    procedure, private :: allocate_comqg
    procedure, private :: init_spectral_coeff
    procedure, private :: init_laplace_helmholtz
    procedure, private :: init_orography
    procedure, private :: init_forcing
    procedure, private :: init_state
    procedure, private :: artiforc
    procedure, private :: dqdt
    procedure, private :: dqdt_d
    procedure, private :: ddt
    procedure, private :: ddt_d
    procedure, private :: ddtl
    procedure, private :: timinp
    procedure, private :: jacob
    procedure, private :: jacob_d
    procedure, private :: jacobd
    procedure, private :: jacobd_d
    procedure, private :: jacobp
    procedure, private :: psitoq
    procedure, private :: psiq
    procedure, private :: qpsi
    procedure, private :: qtopsi
    procedure, private :: lap
    procedure, private :: lapinv
    procedure, private :: fmtofs
    procedure, private :: fstofm
  end type qg_tl_type

  interface qg_tl_type
    procedure :: constructor_qg_tl
  end interface

  ! Mathematical and physical constants
  real(r8kind), parameter :: pi = 4d0 * atan(1d0)     ! value of pi


contains


  !-------------------------------------------------------------------------------
  ! constructor_qg_tl
  !-------------------------------------------------------------------------------
  function constructor_qg_tl(config, state, state_vector, trajectory, trajectory_vector, for, step) result (qg_tl)

    type(qg_config_type),   intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:,:)
    real(r8kind), optional, intent(in) :: state_vector(:)
    real(r8kind), optional, intent(in) :: trajectory(:,:)
    real(r8kind), optional, intent(in) :: trajectory_vector(:)
    real(r8kind), optional, intent(in) :: for(:,:)
    integer,      optional, intent(in) :: step

    type(qg_tl_type)          :: qg_tl
    real(r8kind), allocatable :: psig3d(:,:,:), psisp(:,:)
    real(r8kind), allocatable :: trajg3d(:,:,:)

    real(r8kind) :: nsteps_per_day

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    ! Set the model config
    qg_tl%config = config

    ! Initialize model step and clock
    if (present(step)) then
      qg_tl%step = step
    else
      qg_tl%step = 0
    end if
    qg_tl%clock = qg_tl%step * config%get_time_step()

    ! Initialize time step of the model:
    nsteps_per_day = 24.0d0 * 3600.0d0 / real(config%get_time_step())
    qg_tl%dtt = (1d0 / nsteps_per_day) * pi * 4d0

    ! Set model resolution dependent parameters and allocate model state variables
    call qg_tl%allocate_comqg(config%get_resolution())

    ! Read model input from qgcoefT*
    call qg_tl%init_spectral_coeff()

    ! Instantiate a ggsp grid conversion object
    qg_tl%ggsp = qg_ggsp_type(qg_tl%nm, qg_tl%nlat, qg_tl%nlon, qg_tl%nshm, qg_tl%pp, qg_tl%pd, qg_tl%pw)

    ! Initialize laplace and helmholtz operators
    call qg_tl%init_laplace_helmholtz(config%get_rrdef1(), config%get_rrdef2(), config%get_trel(), config%get_tdis(), config%get_tdif(), config%get_idif())

    ! Initialize orography
    call qg_tl%init_orography(config%get_h0(), config%get_addisl(), config%get_addish())

    ! Initialize model forcing
    if (present(for)) then
      call qg_tl%init_forcing(config%get_obsf(), for=for)
    else
      call qg_tl%init_forcing(config%get_obsf())
    end if

    ! Initialize streamfunction
    if (present(state)) then
      call qg_tl%init_state(psi=state)
    else if (present(state_vector)) then
      allocate(psig3d(qg_tl%nlon, qg_tl%nlat, qg_tl%nvl))
      allocate(psisp(qg_tl%nsh2, qg_tl%nvl))
      psig3d = reshape(state_vector,(/qg_tl%nlon, qg_tl%nlat, qg_tl%nvl/))
      psig3d(:,:,:) = psig3d(:,:,:) / facsf
      psisp(:,1) = reshape(qg_tl%ggsp%ggtosp(transpose(psig3d(:,:,1))), (/qg_tl%nsh2/))
      psisp(:,2) = reshape(qg_tl%ggsp%ggtosp(transpose(psig3d(:,:,2))), (/qg_tl%nsh2/))
      psisp(:,3) = reshape(qg_tl%ggsp%ggtosp(transpose(psig3d(:,:,3))), (/qg_tl%nsh2/))
      call qg_tl%init_state(psi=psisp)
    else
      call qg_tl%init_state()
    end if

    ! Initialize trajectory
    if (present(trajectory)) then
      allocate(qg_tl%trajectory, source = trajectory)
    else if (present(trajectory_vector)) then
      allocate(trajg3d(qg_tl%nlon, qg_tl%nlat, qg_tl%nvl))
      allocate(qg_tl%trajectory(qg_tl%nsh2,qg_tl%nvl))
      trajg3d = reshape(trajectory_vector,(/qg_tl%nlon, qg_tl%nlat, qg_tl%nvl/))
      trajg3d(:,:,:) = trajg3d(:,:,:) / facsf
      qg_tl%trajectory(:,1) = reshape(qg_tl%ggsp%ggtosp(transpose(trajg3d(:,:,1))), (/qg_tl%nsh2/))
      qg_tl%trajectory(:,2) = reshape(qg_tl%ggsp%ggtosp(transpose(trajg3d(:,:,2))), (/qg_tl%nsh2/))
      qg_tl%trajectory(:,3) = reshape(qg_tl%ggsp%ggtosp(transpose(trajg3d(:,:,3))), (/qg_tl%nsh2/))
    else
      allocate(qg_tl%trajectory, source = qg_tl%psi)
    end if
    allocate(qg_tl%psitd(qg_tl%nsh2,qg_tl%ntl))
    allocate(qg_tl%qprimed(qg_tl%nsh2,qg_tl%nvl))
    call qg_tl%psitoq(qg_tl%trajectory, qg_tl%psitd, qg_tl%qprimed)

  end function constructor_qg_tl


  !------------------------------------------------------------------
  ! destructor_qg_tl
  !
  ! Deallocates pointers used by a qg_tl_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_qg_tl(this)


    type(qg_tl_type), intent(inout) :: this

    ! No pointers in qg_tl_type object so we do nothing

  end subroutine destructor_qg_tl


  !-------------------------------------------------------------------------------
  ! allocate_comqg
  !-------------------------------------------------------------------------------
  subroutine allocate_comqg(this, resolution)

    class(qg_tl_type) :: this
    integer, intent(in) :: resolution

    select case (resolution)

      case(106)
        this%nm = 106
        this%nlon = 320
        this%nlat = 160
        this%ft = "106"

      case(63)
        this%nm = 63
        this%nlon = 192
        this%nlat = 96
        this%ft = "63"

      case(42)
        this%nm = 42
        this%nlon = 128
        this%nlat = 64
        this%ft = "42"

      case(21)
        this%nm = 21
        this%nlon = 64
        this%nlat = 32
        this%ft = "21"

      case DEFAULT
        stop 'ERROR: Unsupported resolution'

    end select

    this%nvl = 3
    this%ntl = this%nvl - 1
    this%nsh = ((this%nm + 1) * (this%nm + 2)) / 2
    this%nsh2 = 2 * this%nsh
    this%nsh23 = this%nsh2 * this%nvl

  end subroutine allocate_comqg


  !-------------------------------------------------------------------------------
  ! init_spectral_coeff
  !-------------------------------------------------------------------------------
  subroutine init_spectral_coeff(this)

    class(qg_tl_type), intent(inout) :: this

    integer      :: i, j, k
    real(r8kind) :: sqn, rsqn

    ! Allocate spectral coefficient arrays
    allocate(this%nshm(0:this%nm))
    allocate(this%ll(this%nsh))
    allocate(this%pp(this%nlat,this%nsh))
    allocate(this%pd(this%nlat,this%nsh))
    allocate(this%pw(this%nlat,this%nsh))

    ! Read spectral coefficients from qgcoefT*
    open(11, file = './qgcoefT' // trim(this%ft) // '.dat', form = 'formatted')
    do i = 0, this%nm
      read(11, *) this%nshm(i)
    enddo
    do i = 1, this%nsh
      read(11, *) this%ll(i)
    enddo
    do k = 1, this%nsh
      do j = 1, this%nlat
        read(11, *) this%pp(j, k)
      enddo
    enddo
    do k = 1, this%nsh
      do j = 1, this%nlat
        read(11, *) this%pd(j, k)
      enddo
    enddo
    do k = 1, this%nsh
      do j = 1, this%nlat
        read(11, *) this%pw(j, k)
      enddo
    enddo
    close(11)

    ! compensation for normalization in nag fft routines
    sqn = sqrt(dble(this%nlon))
    rsqn = 1d0 / sqn
    do k = 1, this%nsh
      do i = 1, this%nlat
        this%pp(i, k) = this%pp(i, k) * sqn
        this%pd(i, k) = this%pd(i, k) * sqn
        this%pw(i, k) = this%pw(i, k) * rsqn
      enddo
    enddo

  end subroutine init_spectral_coeff


  !-------------------------------------------------------------------------------
  ! init_laplace_helmholtz
  !-------------------------------------------------------------------------------
  subroutine init_laplace_helmholtz(this, rrdef1, rrdef2, trel, tdis, tdif, idif)

    class(qg_tl_type), intent(inout) :: this
    real(r8kind) :: rrdef1
    real(r8kind) :: rrdef2
    real(r8kind) :: trel
    real(r8kind) :: tdis
    real(r8kind) :: tdif
    integer      :: idif

    real(r8kind) :: pigr4, dis, rll, dif
    real(r8kind) :: r1, a, b, c, d, e
    integer      :: j, k

    pigr4 = 4.d0 * pi
    this%rl1 = 1.0d0 / rrdef1**2
    this%rl2 = 1.0d0 / rrdef2**2
    this%relt1 = max(0.0d0, this%rl1 / (trel * pigr4))
    this%relt2 = max(0.0d0, this%rl2 / (trel * pigr4))
    dis = max(0.0d0, 1.0d0 / (tdis * pigr4))
    rll = dble(this%ll(this%nsh))
    dif = max(0.0d0, 1.0d0 / (tdif * pigr4 * (rll * (rll + 1))**idif))

    ! Allocate laplace helmholtz arrays
    allocate(this%rinhel(this%nsh2,0:5))
    allocate(this%diss(this%nsh2,2))

    ! laplace/helmholtz direct and inverse operators
    do j = 0, 5
      this%rinhel(1, j) = 0.0d0
    enddo
    this%diss(1, 1) = 0.0d0
    this%diss(1, 2) = 0.0d0
    do k = 2, this%nsh
      r1 = dble(this%ll(k) * (this%ll(k) + 1))
      a = -r1 - 3.0d0 * this%rl1
      b = -r1 - 3.0d0 * this%rl2
      c = -r1 - this%rl1
      d = -r1 - this%rl2
      e = a * d + b * c
      this%rinhel(k, 0) = -r1
      this%rinhel(k, 1) = -1.0d0 / r1
      this%rinhel(k, 2) =  d / e
      this%rinhel(k, 3) =  b / e
      this%rinhel(k, 4) = -c / e
      this%rinhel(k, 5) =  a / e
      this%diss(k, 2) = dis * r1
      this%diss(k, 1) = -dif * r1**idif
    enddo
    do j = 0, 5
      do k = 1, this%nsh
        this%rinhel(k + this%nsh, j) = this%rinhel(k, j)
      enddo
    enddo
    do j = 1, 2
      do k = 1, this%nsh
        this%diss(k + this%nsh, j) = this%diss(k, j)
      enddo
    enddo

  end subroutine init_laplace_helmholtz


  !-------------------------------------------------------------------------------
  ! init_orography
  !-------------------------------------------------------------------------------
  subroutine init_orography(this, h0, addisl, addish)

    class(qg_tl_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: h0
    real(r8kind),         intent(   in) :: addisl
    real(r8kind),         intent(   in) :: addish

    real(r8kind), allocatable :: fmu(:,:)
    real(r8kind), allocatable :: agg(:,:), agg1(:,:), agg2(:,:)
!    real(r8kind), allocatable :: orog(:)
    real(r8kind), allocatable :: ws(:), wsx(:)
    real(r8kind)              :: rnorm, dlon, rh0, dd
    integer                   :: i, j

    ! orography and dissipation terms
    ! fmu(i, 1): sin(phi(i))
    ! fmu(i, 2): 1 - sin**2(phi(i))      
    allocate(fmu(this%nlat, 2))
    allocate(this%phi(this%nlat))
    allocate(this%sinfi(this%nlat))
    allocate(this%cosfi(this%nlat))
    rnorm = 1.0d0 / sqrt(3.0d0 * this%nlon)
    do i = 1, this%nlat
      fmu(i, 1) = rnorm * this%pp(i, 2)
      fmu(i, 2) = 1.d0 - fmu(i, 1)**2
      this%sinfi(i) = fmu(i, 1)
      this%phi(i) = asin(this%sinfi(i))
      this%cosfi(i) = cos(this%phi(i))
      this%phi(i) = 180d0 * this%phi(i) / pi
    enddo
    dlon = 360d0 / real(this%nlon)

    ! height of orography in meters
    allocate(agg(this%nlat, this%nlon))
    allocate(agg1(this%nlat, this%nlon))
    open(13, file = './qgbergT' // trim(this%ft) // '.dat', form = 'formatted')
    do i = 1, this%nlon
      do j = 1, this%nlat
        read(13, *) agg1(j, i)
      enddo
    enddo
    rh0 = max(0.0d0, 0.001d0 / h0)
    do j = 1, this%nlon
      do i = 1, this%nlat
        agg(i, j) = fmu(i, 1) * agg1(i, j) * rh0
    !      agg(i, j) = agg1(i, j) * rh0
      enddo
    enddo

    ! surface dependent friction
    allocate(this%orog(this%nsh2))
    allocate(ws(this%nsh2))
    allocate(this%dorodl(this%nlat,this%nlon))
    allocate(this%dorodm(this%nlat,this%nlon))
    allocate(agg2(this%nlat, this%nlon))
    this%lgdiss = ((addisl .gt. 0.0) .or. (addish .gt. 0.0))
    this%orog = reshape(this%ggsp%ggtosp (agg), (/this%nsh2/))
    ws = reshape(this%ggsp%ddl (this%orog), (/this%nsh2/))
    this%dorodl = this%ggsp%sptogg_pp (ws)
    this%dorodm = this%ggsp%sptogg_pd (this%orog)
    if (this%lgdiss) then
      do i = 1, this%nlon
        do j = 1, this%nlat
          read(13, *) agg2(j, i)
        enddo
      enddo
      do j = 1, this%nlon
        do i = 1, this%nlat
          agg(i, j) = 1.0d0 + addisl * agg2(i, j) + addish * (1.0d0 - exp(-0.001d0 * agg1(i, j)))
        enddo
      enddo

      allocate(wsx(this%nsh2))
      allocate(this%rdiss(this%nlat,this%nlon))
      allocate(this%ddisdx(this%nlat,this%nlon))
      allocate(this%ddisdy(this%nlat,this%nlon))
      ws = reshape(this%ggsp%ggtosp (agg), (/this%nsh2/))
      wsx = reshape(this%ggsp%ddl (ws), (/this%nsh2/))
      this%rdiss = this%ggsp%sptogg_pp (ws)
      this%ddisdx = this%ggsp%sptogg_pp (wsx)
      this%ddisdy = this%ggsp%sptogg_pd (ws)
      dd = 0.5d0 * this%diss(2, 2)
      do j = 1, this%nlon
        do i = 1, this%nlat
          this%ddisdx(i, j) = dd * this%ddisdx(i, j) / fmu(i, 2)
          this%ddisdy(i, j) = dd * this%ddisdy(i, j) * fmu(i, 2)
        enddo
      enddo

    endif
    close(13)

    open(13, file = 'qgbergT' // trim(this%ft) // '.grads', form = 'unformatted')
    write(13) ((real(agg1(j, i)), i = 1, this%nlon), j = 1, this%nlat)
    write(13) ((real(agg2(j, i)), i = 1, this%nlon), j = 1, this%nlat)
    close(13)

    open(50, file = 'qgbergT' // trim(this%ft) // '.ctl', form = 'formatted')
    write(50, '(A)') 'dset ^qgbergT' // trim(this%ft) // '.grads'
    write(50, '(A)') 'undef 9.99e+10'
    write(50, '(A)') 'options sequential big_endian'
    write(50, '(A)') 'title three level QG model'
    write(50, '(A)') '*'
    write(50, '(A, i4, A, F19.14)') 'xdef ', this%nlon, ' linear  0.000 ', dlon
    write(50, '(A)') '*'
    write(50, '(A, I4, A, 1F19.14)') 'ydef ', this%nlat, ' levels ', this%phi(1)
    write(50, '(F19.14)') (this%phi(j), j = 2, this%nlat)
    write(50, '(A)') '*'
    write(50, '(A)') 'zdef  1 levels 1000'
    write(50, '(A)') '*'
    write(50, '(A)') 'tdef 1 linear 1jan0001 1dy'
    write(50, '(A)') '*'
    write(50, '(A)') 'vars  2'
    write(50, '(A)') 'oro    1  99 orography [m]'
    write(50, '(A)') 'friction    1  99 friction mask'
    write(50, '(A)') 'endvars'
    close(50)

  end subroutine init_orography


  !-------------------------------------------------------------------------------
  ! init_forcing
  !-------------------------------------------------------------------------------
  subroutine init_forcing(this, obsf, for)

    use netcdf
    use netcdf_utilities

    class(qg_tl_type),   intent(inout) :: this
    logical,                intent(   in) :: obsf
    real(r8kind), optional, intent(   in) :: for(:,:)

    ! General netCDF variables
    integer :: ncFileID      ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: ForVarID

    ! Initialize forcing
    if (present(for)) then

      ! Initialize forcing to optional input
      allocate(this%for(this%nsh2,this%nvl), source=for)

    elseif (obsf) then
      ! Calculate an artificial forcing from observations if requested
      allocate(this%for(this%nsh2,this%nvl))
      this%for = this%artiforc()
    else

      ! Initialize the forcing from a bootstrap file
      allocate(this%for(this%nsh2,this%nvl))

      ! Open file for read only
      call nc_check(nf90_open('qgforcingT' // trim(this%ft) // '.nc', NF90_NOWRITE, ncFileID))
      call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

      ! Get the streamfunction ID
      call nc_check(nf90_inq_varid(ncFileID, "For", ForVarID))

      ! Get the streamfunction variable
      call nc_check(nf90_get_var(ncFileID, ForVarID, this%for))

      ! Flush buffers
      call nc_check(nf90_sync(ncFileID))

      ! Close the NetCDF file
      call nc_check(nf90_close(ncFileID))

!      do l = 1, this%nvl
!        do k = 1, this%nsh2
!          this%for(k, l) = 0d0
!        enddo
!      enddo

   endif

  end subroutine init_forcing


  !-------------------------------------------------------------------------------
  ! init_state
  !-------------------------------------------------------------------------------
  subroutine init_state(this, psi)

    use netcdf
    use netcdf_utilities

    class(qg_tl_type),   intent(inout) :: this
    real(r8kind), optional, intent(   in) :: psi(:,:)

    ! General netCDF variables
    integer :: ncFileID      ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: PsiVarID

    if (present(psi)) then

      ! Initialize streamfunction to optional input
      allocate(this%psi(this%nsh2,this%nvl), source=psi)

    else

      ! Read initial streamfunction from a bootstrap file
      allocate(this%psi(this%nsh2,this%nvl))

      ! Open file for read only
      call nc_check(nf90_open('qginitT' // trim(this%ft) // '.nc', NF90_NOWRITE, ncFileID))
      call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

      ! Get the streamfunction ID
      call nc_check(nf90_inq_varid(ncFileID, "Psi", PsiVarID))

      ! Get the streamfunction variable
      call nc_check(nf90_get_var(ncFileID, PsiVarID, this%psi))

      ! Flush buffers
      call nc_check(nf90_sync(ncFileID))

      ! Close the NetCDF file
      call nc_check(nf90_close(ncFileID))

!      do l = 1, this%nvl
!        do k = 1, this%nsh2
!          this%psi(k, l) = 0d0
!        enddo
!      enddo

    endif

    ! Allocate streamfunction and layer thickness arrays
    allocate(this%psit(this%nsh2,this%ntl))

    ! Allocate potential vorticity array
    allocate(this%qprime(this%nsh2,this%nvl))

    ! Initialize potential vorticity from streamfunction
    call this%psitoq(this%psi, this%psit, this%qprime)

  end subroutine init_state


  !-----------------------------------------------------------------------
  ! performs a fourth order runge kutta time step at truncation nm
  ! with time step dt
  ! dqdt calculates the time derivative
  ! input  qprime at current time
  ! output qprime at current time plus dt
  !-----------------------------------------------------------------------
  subroutine adv_nsteps(this, nsteps)

    class(qg_tl_type) :: this
    integer           :: nsteps

    integer :: step, k, l, nvar
    real(r8kind) :: dt2, dt6
    real(r8kind) :: y(this%nsh2, this%nvl), dydt(this%nsh2, this%nvl), yt(this%nsh2, this%nvl)
    real(r8kind) :: dyt(this%nsh2, this%nvl), dym(this%nsh2, this%nvl)
    real(r8kind) :: yd(this%nsh2, this%nvl), dydtd(this%nsh2, this%nvl), ytd(this%nsh2, this%nvl)
    real(r8kind) :: dytd(this%nsh2, this%nvl), dymd(this%nsh2, this%nvl)

    if (nsteps > 0) then

      nvar = (this%nm + 2) * this%nm
      dt2 = this%dtt * 0.5d0
      dt6 = this%dtt / 6d0

      ! Advance the model forward in time n steps
      do step = 1, nsteps

        yd = this%fmtofs(this%qprimed)
        y = this%fmtofs(this%qprime)
        CALL this%DQDT_D(y, yd, dydt, dydtd)
        ytd = 0.0_8
        DO l = 1, this%nvl
          DO k = 1, nvar
            ytd(k, l) = yd(k, l) + dt2 * dydtd(k, l)
            yt(k, l) = y(k, l) + dt2 * dydt(k, l)
          END DO
        END DO
        CALL this%DQDT_D(yt, ytd, dyt, dytd)
        DO l = 1, this%nvl
          DO k = 1, nvar
            ytd(k, l) = yd(k, l) + dt2 * dytd(k, l)
            yt(k, l) = y(k, l) + dt2 * dyt(k, l)
          END DO
        END DO
        CALL this%DQDT_D(yt, ytd, dym, dymd)
        DO l = 1, this%nvl
          DO k = 1, nvar
            ytd(k, l) = yd(k, l) + this%dtt * dymd(k, l)
             yt(k, l) = y(k, l) + this%dtt * dym(k, l)
            dymd(k, l) = dytd(k, l) + dymd(k, l)
             dym(k, l) = dyt(k, l) + dym(k, l)
          END DO
        END DO
        CALL this%DQDT_D(yt, ytd, dyt, dytd)
        DO l = 1, this%nvl
          DO k = 1, nvar
            yd(k, l) = yd(k, l) + dt6 * (dydtd(k, l) + dytd(k, l) + 2. * dymd(k, l))
            y(k, l)  = y(k, l) + dt6 * ( dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
          END DO
        END DO

        ! Update model state with original trajectory
        this%qprime = this%qprime + this%qprimed

        ! Update trajectory
        this%qprimed = this%fstofm(yd, this%nm)

        ! Inrement the step count
        this%step = this%step + 1

      end do

      ! Make stream function consistent with potential vorticity
      call this%qtopsi(this%qprimed, this%trajectory, this%psitd)
      call this%qtopsi(this%qprime, this%psi, this%psit)

    end if

  end subroutine adv_nsteps


  SUBROUTINE DQDT_D(this, y, yd, dydt, dydtd)

    class(qg_tl_type) :: this
    real(r8kind), INTENT( IN) :: y(:, :)
    real(r8kind), INTENT( IN) :: yd(:, :)
    real(r8kind), INTENT(OUT) :: dydt(:, :)
    real(r8kind), INTENT(OUT) :: dydtd(:, :)

    real(r8kind) :: qprime(this%nsh2,this%nvl) ! qprime
    real(r8kind) :: qprimed(this%nsh2,this%nvl) ! qprime
    real(r8kind) :: psi(this%nsh2,this%nvl)    ! psi
    real(r8kind) :: psid(this%nsh2,this%nvl)    ! psi
    real(r8kind) :: psit(this%nsh2,this%ntl)   ! psit
    real(r8kind) :: psitd(this%nsh2,this%ntl)   ! psit
    real(r8kind) :: dqprdt(this%nsh2, this%nvl)  ! time derivative of qprime
    real(r8kind) :: dqprdtd(this%nsh2, this%nvl) ! time derivative of qprime

    qprimed = this%fstofm(yd, this%nm)
    qprime = this%fstofm(y, this%nm)

    call this%qtopsi(qprimed, psid, psitd) ! qprime --> psi and psit
    call this%qtopsi(qprime, psi, psit)          ! qprime --> psi and psit

    ! psi, psit, qprime, for, diss --> dqprdt
    call this%DDT_D(psi, psid, psit, psitd, qprime, qprimed, this%for, dqprdt, dqprdtd)

    dydtd = this%fmtofs(dqprdtd)
    dydt = this%fmtofs(dqprdt)


  END SUBROUTINE DQDT_D

  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  subroutine dqdt(this, y, dydt)

    class(qg_tl_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: y(:,:)
    real(r8kind),         intent(  out) :: dydt(:,:)

    real(r8kind) :: qprime(this%nsh2,this%nvl) ! qprime
    real(r8kind) :: psi(this%nsh2,this%nvl)    ! psi
    real(r8kind) :: psit(this%nsh2,this%nvl)   ! psit
    real(r8kind) :: dqprdt(this%nsh2,this%nvl) ! time derivative of qprime

    qprime = this%fstofm(y, this%nm)
    call this%qtopsi(qprime, psi, psit)
    dqprdt = this%ddt(psi, psit, qprime, this%for) ! psi, psit, qprime, for, diss --> dqprdt
    dydt = this%fmtofs(dqprdt)

    return

  end subroutine dqdt


  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  subroutine DDT_D(this, psi, psid, psit, psitd, qprime, qprimed, for, dqprdt, dqprdtd)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),      intent( in) :: psi(this%nsh2, this%nvl)    ! stream function at the nvl levels
    real(r8kind),      intent( in) :: psid(this%nsh2, this%nvl)    ! stream function at the nvl levels
    real(r8kind),      intent( in) :: psit(this%nsh2, this%ntl)   ! thickness at the ntl levels
    real(r8kind),      intent( in) :: psitd(this%nsh2, this%ntl)   ! thickness at the ntl levels
    real(r8kind),      intent( in) :: qprime(this%nsh2, this%nvl) ! potential vorticity
    real(r8kind),      intent( in) :: qprimed(this%nsh2, this%nvl) ! potential vorticity
    real(r8kind),      intent( in) :: for(this%nsh2, this%nvl)    ! constant potential vorticity forcing at the nvl levels
    real(r8kind),      intent(out) :: dqprdt(this%nsh2, this%nvl)
    real(r8kind),      intent(out) :: dqprdtd(this%nsh2, this%nvl)

    integer      :: k, l, i, j
    real(r8kind) :: dum1, dum2
    real(r8kind) :: dum1d, dum2d

    dqprdtd = 0.0_8

    ! advection of potential vorticity at upper level
    call this%JACOB_D(psi(:, 1), psid(:, 1), qprime(:, 1), qprimed(:, 1), dqprdt(:, 1), dqprdtd(:, 1))

    ! advection of potential vorticity at middle level
    call this%JACOB_D(psi(:, 2), psid(:, 2), qprime(:, 2), qprimed(:, 2), dqprdt(:, 2), dqprdtd(:, 2))

    ! advection of potential vorticity and dissipation at lower level
    call this%JACOBD_D(psi(:, 3), psid(:, 3), qprime(:, 3), qprimed(:, 3), dqprdt(:, 3), dqprdtd(:, 3))

    ! relaxation of temperature and forcing
    DO k = 1, this%nsh2
      dum1d = this%relt1 * psitd(k, 1)
      dum1 = this%relt1 * psit(k, 1)
      dum2d = this%relt2 * psitd(k, 2)
      dum2 = this%relt2 * psit(k, 2)
      dqprdtd(k, 1) = dqprdtd(k, 1) + dum1d
      dqprdt(k, 1) = dqprdt(k, 1) + dum1 + for(k, 1)
      dqprdtd(k, 2) = dqprdtd(k, 2) - dum1d + dum2d
      dqprdt(k, 2) = dqprdt(k, 2) - dum1 + dum2 + for(k, 2)
      dqprdtd(k, 3) = dqprdtd(k, 3) - dum2d
      dqprdt(k, 3) = dqprdt(k, 3) - dum2 + for(k, 3)
    END DO

    ! explicit horizontal diffusion
    DO l = 1, this%nvl
      DO k = 1, this%nsh2
        dqprdtd(k, l) = dqprdtd(k, l) + this%diss(k, 1) * qprimed(k, l)
        dqprdt(k, l) = dqprdt(k, l) + this%diss(k, 1) * qprime(k, l)
      END DO
    END DO

  END subroutine DDT_D


  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  function ddt(this, psi, psit, qprime, for) result(dqprdt)

    class(qg_tl_type), intent(in) :: this
    real(r8kind),      intent(in) :: psi(this%nsh2,this%nvl)    ! stream function at the nvl levels
    real(r8kind),      intent(in) :: psit(this%nsh2,this%ntl)   ! thickness at the ntl levels
    real(r8kind),      intent(in) :: qprime(this%nsh2,this%nvl) ! potential vorticity
    real(r8kind),      intent(in) :: for(this%nsh2,this%nvl)    ! constant potential vorticity forcing at the nvl levels
    real(r8kind)                  :: dqprdt(this%nsh2,this%nvl)

    integer :: k, l, i, j
    real(r8kind) :: dum1, dum2

    ! advection of potential vorticity at upper level
    dqprdt(:, 1) = reshape(this%jacob (psi(:, 1), qprime(:, 1)), (/this%nsh2/))

    ! advection of potential vorticity at middle level
    dqprdt(:, 2) = reshape(this%jacob (psi(:, 2), qprime(:, 2)), (/this%nsh2/))

    ! advection of potential vorticity and dissipation at lower level
    dqprdt(:, 3) = reshape(this%jacobd (psi(:, 3), qprime(:, 3)), (/this%nsh2/))

    ! relaxation of temperature and forcing
    do k = 1, this%nsh2
      dum1 = this%relt1 * psit(k, 1)
      dum2 = this%relt2 * psit(k, 2)
      dqprdt(k, 1) = dqprdt(k, 1) + dum1        + for(k, 1)
      dqprdt(k, 2) = dqprdt(k, 2) - dum1 + dum2 + for(k, 2)
      dqprdt(k, 3) = dqprdt(k, 3)        - dum2 + for(k, 3)
    enddo

    ! explicit horizontal diffusion
    do l = 1, 3
      do k = 1, this%nsh2
        dqprdt(k, l) = dqprdt(k, l) + this%diss(k, 1) * qprime(k, l)
      enddo
    enddo

    return

  end function ddt


  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  subroutine JACOB_D(this, psiloc, psilocd, pvor, pvord, sjacob, sjacobd)

    class(qg_tl_type), intent(in) :: this
    real(r8kind), intent( in) :: psiloc(this%nsh2)
    real(r8kind), intent( in) :: psilocd(this%nsh2)
    real(r8kind), intent( in) :: pvor(this%nsh2)
    real(r8kind), intent( in) :: pvord(this%nsh2)
    real(r8kind), intent(out) :: sjacob(this%nsh2)
    real(r8kind), intent(out) :: sjacobd(this%nsh2)

    integer      :: i, j, k
    real(r8kind) :: vv(this%nsh2)
    real(r8kind) :: vvd(this%nsh2)
    real(r8kind) :: dpsidl(this%nlat, this%nlon),  dpsidm(this%nlat, this%nlon),  dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidld(this%nlat, this%nlon),  dpsidmd(this%nlat, this%nlon),  dvordld(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon),  gjacob(this%nlat, this%nlon),  dpsidls(this%nsh2)
    real(r8kind) :: dvordmd(this%nlat, this%nlon),  gjacobd(this%nlat, this%nlon),  dpsidlsd(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vvd = reshape(ggsp%DDL(pvord), (/this%nsh2/))
    vv = reshape(ggsp%DDL(pvor), (/this%nsh2/))
    dvordld = ggsp%SPTOGG_PP(vvd)
    dvordl = ggsp%SPTOGG_PP(vv)
    dvordmd = ggsp%SPTOGG_PD(pvord)
    dvordm = ggsp%SPTOGG_PD(pvor)

    ! space derivatives of streamfunction
    dpsidlsd = reshape(ggsp%DDL(psilocd), (/this%nsh2/))
    dpsidls = reshape(ggsp%DDL(psiloc), (/this%nsh2/))
    dpsidld = ggsp%SPTOGG_PP(dpsidlsd)
    dpsidl = ggsp%SPTOGG_PP(dpsidls)
    dpsidmd = ggsp%SPTOGG_PD(psilocd)
    dpsidm = ggsp%SPTOGG_PD(psiloc)

    gjacobd = 0.0_8

    ! jacobian term
    DO j = 1, this%nlon
      DO i = 1, this%nlat
        gjacobd(i, j) = dpsidmd(i, j) * dvordl(i, j) + dpsidm(i, j) * dvordld(i, j) - dpsidld(i, j) * dvordm(i, j) - dpsidl(i, j) * dvordmd(i, j)
        gjacob(i, j) = dpsidm(i, j) * dvordl(i, j) - dpsidl(i, j) * dvordm(i, j)
      END DO
    END DO
    sjacobd = reshape(ggsp%GGTOSP(gjacobd), (/this%nsh2/))
    sjacob = reshape(ggsp%GGTOSP(gjacob), (/this%nsh2/))

    ! planetary vorticity advection
    DO k = 1, this%nsh2
      sjacobd(k) = sjacobd(k) - dpsidlsd(k)
      sjacob(k) = sjacob(k) - dpsidls(k)
    END DO

  END subroutine JACOB_D


  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacob (this, psiloc, pvor) result(sjacob)

    implicit none

    class(qg_tl_type), intent(in) :: this
    real(r8kind), intent( in) :: psiloc(this%nsh2)
    real(r8kind), intent( in) :: pvor(this%nsh2)
    real(r8kind)              :: sjacob(this%nsh2)

    integer      :: i, j, k
    real(r8kind) :: vv(this%nsh2)
    real(r8kind) :: dpsidl(this%nlat, this%nlon),  dpsidm(this%nlat, this%nlon),  dvordl(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon),  gjacob(this%nlat, this%nlon),  dpsidls(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vv = reshape(ggsp%ddl (pvor), (/this%nsh2/))
    dvordl = ggsp%sptogg_pp (vv)
    dvordm = ggsp%sptogg_pd (pvor)

    ! space derivatives of streamfunction
    dpsidls = reshape(ggsp%ddl (psiloc), (/this%nsh2/))
    dpsidl = ggsp%sptogg_pp (dpsidls)
    dpsidm = ggsp%sptogg_pd (psiloc)

    ! jacobian term
    do j = 1, this%nlon
      do i = 1, this%nlat
        gjacob(i, j) = dpsidm(i, j) * dvordl(i, j) - dpsidl(i, j) * dvordm(i, j)
      enddo
    enddo

    sjacob = reshape(ggsp%ggtosp (gjacob), (/this%nsh2/))

    ! planetary vorticity advection
    do k = 1, this%nsh2
      sjacob(k) = sjacob(k) - dpsidls(k)
    enddo

    return

  end function jacob


  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  subroutine JACOBD_D(this, psiloc, psilocd, pvor, pvord, sjacob, sjacobd)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),      intent( in) :: psiloc(this%nsh2)
    real(r8kind),      intent( in) :: psilocd(this%nsh2)
    real(r8kind),      intent( in) :: pvor(this%nsh2)
    real(r8kind),      intent( in) :: pvord(this%nsh2)
    real(r8kind),      intent(out) :: sjacob(this%nsh2)
    real(r8kind),      intent(out) :: sjacobd(this%nsh2)

    integer      :: i, j, k
    real(r8kind) :: dpsidl(this%nlat, this%nlon),  dpsidm(this%nlat, this%nlon),  dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidld(this%nlat, this%nlon),  dpsidmd(this%nlat, this%nlon),  dvordld(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon),  gjacob(this%nlat, this%nlon),  vv(this%nsh2)
    real(r8kind) :: dvordmd(this%nlat, this%nlon),  gjacobd(this%nlat, this%nlon),  vvd(this%nsh2)
    real(r8kind) :: azeta(this%nlat, this%nlon), dpsidls(this%nsh2)
    real(r8kind) :: azetad(this%nlat, this%nlon), dpsidlsd(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vvd = reshape(ggsp%DDL(pvord), (/this%nsh2/))
    vv = reshape(ggsp%DDL(pvor), (/this%nsh2/))
    dvordld = ggsp%SPTOGG_PP(vvd)
    dvordl = ggsp%SPTOGG_PP(vv)
    dvordmd = ggsp%SPTOGG_PD(pvord)
    dvordm = ggsp%SPTOGG_PD(pvor)

    ! space derivatives of streamfunction
    dpsidlsd = reshape(ggsp%DDL(psilocd), (/this%nsh2/))
    dpsidls = reshape(ggsp%DDL(psiloc), (/this%nsh2/))
    dpsidld = ggsp%SPTOGG_PP(dpsidlsd)
    dpsidl = ggsp%SPTOGG_PP(dpsidls)
    dpsidmd = ggsp%SPTOGG_PD(psilocd)
    dpsidm = ggsp%SPTOGG_PD(psiloc)

    gjacobd = 0.0_8

    ! jacobian term + orographic forcing
    DO j = 1, this%nlon
      DO i = 1, this%nlat
        gjacobd(i, j) = dpsidmd(i, j) * (dvordl(i, j) + this%sinfi(i) * this%dorodl(i, j)) + dpsidm(i, j) * dvordld(i, j) - &
                     &  dpsidld(i, j) * (dvordm(i, j) + this%sinfi(i) * this%dorodm(i, j)) - dpsidl(i, j) * dvordmd(i, j)
        gjacob(i, j) = dpsidm(i, j) * (dvordl(i, j) + this%sinfi(i) * this%dorodl(i, j)) - dpsidl(i, j) * (dvordm(i, j) + &
                     & this%sinfi(i) * this%dorodm(i, j))
      END DO
    END DO

    ! dissipation 
    IF (this%lgdiss) THEN

      !   spatially varying dissipation 
      DO k = 1, this%nsh2
        vvd(k) = this%diss(k, 2) * psilocd(k)
        vv(k) = this%diss(k, 2) * psiloc(k)
      END DO

      azetad = ggsp%SPTOGG_PP(vvd)
      azeta = ggsp%SPTOGG_PP(vv)

      DO j = 1, this%nlon
        DO i = 1, this%nlat
          gjacobd(i, j) = gjacobd(i, j) - this%ddisdy(i, j) * dpsidmd(i, j) - this%ddisdx(i, j) * dpsidld(i, j) + this%rdiss(i, j) * azetad(i, j)
          gjacob(i, j) = gjacob(i, j) - dpsidm(i, j) * this%ddisdy(i, j) - dpsidl(i, j) * this%ddisdx(i, j) + this%rdiss(i, j) * azeta(i, j)
        END DO
      END DO

      sjacobd = reshape(ggsp%GGTOSP(gjacobd), (/this%nsh2/))
      sjacob = reshape(ggsp%GGTOSP(gjacob), (/this%nsh2/))

    ELSE

      !   uniform dissipation
      sjacobd = reshape(ggsp%GGTOSP(gjacobd), (/this%nsh2/))
      sjacob = reshape(ggsp%GGTOSP(gjacob), (/this%nsh2/))
      DO k = 1, this%nsh2
        sjacobd(k) = sjacobd(k) + this%diss(k, 2) * psilocd(k)
        sjacob(k) = sjacob(k) + this%diss(k, 2) * psiloc(k)
      END DO

    END IF

    ! planetary vorticity advection
    DO k = 1, this%nsh2
      sjacobd(k) = sjacobd(k) - dpsidlsd(k)
      sjacob(k) = sjacob(k) - dpsidls(k)
    END DO

  END subroutine JACOBD_D


  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacobd (this, psiloc, pvor) result(sjacob)

    class(qg_tl_type), intent(in) :: this
    real(r8kind),         intent(in) :: psiloc(this%nsh2)
    real(r8kind),         intent(in) :: pvor(this%nsh2)
    real(r8kind)                     :: sjacob(this%nsh2)

    integer      :: i, j, k
    real(r8kind) :: dpsidl(this%nlat, this%nlon),  dpsidm(this%nlat, this%nlon),  dvordl(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon),  gjacob(this%nlat, this%nlon),  vv(this%nsh2)
    real(r8kind) :: azeta(this%nlat, this%nlon), dpsidls(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity 
    vv = reshape(ggsp%ddl (pvor), (/this%nsh2/))
    dvordl = ggsp%sptogg_pp (vv)
    dvordm = ggsp%sptogg_pd (pvor)

    ! space derivatives of streamfunction
    dpsidls = reshape(ggsp%ddl (psiloc), (/this%nsh2/))
    dpsidl = ggsp%sptogg_pp (dpsidls)
    dpsidm = ggsp%sptogg_pd (psiloc)

    ! jacobian term + orographic forcing
    do j = 1, this%nlon
      do i = 1, this%nlat
        gjacob(i, j) = dpsidm(i, j) * (dvordl(i, j) + this%sinfi(i) * this%dorodl(i, j)) -  &
   &                   dpsidl(i, j) * (dvordm(i, j) + this%sinfi(i) * this%dorodm(i, j))
      enddo
    enddo

    ! dissipation 
    if (this%lgdiss) then

      !   spatially varying dissipation 
      do k = 1, this%nsh2
        vv(k) = this%diss(k, 2) * psiloc(k)
      enddo

      azeta = ggsp%sptogg_pp (vv)

      do j = 1, this%nlon
        do i = 1, this%nlat
          gjacob(i, j) = gjacob(i, j) - dpsidm(i, j)     * this%ddisdy(i, j) &
   &                                  - dpsidl(i, j)     * this%ddisdx(i, j) &
   &                                  + this%rdiss(i, j) * azeta(i, j)
        enddo
      enddo

      sjacob = reshape(ggsp%ggtosp (gjacob), (/this%nsh2/))

    else

      !   uniform dissipation
      sjacob = reshape(ggsp%ggtosp (gjacob), (/this%nsh2/))

      do k = 1, this%nsh2
        sjacob(k) = sjacob(k) + this%diss(k, 2) * psiloc(k)
      enddo

    endif

    ! planetary vorticity advection
    do k = 1, this%nsh2
      sjacob(k) = sjacob(k) - dpsidls(k)
    enddo

    return

  end function jacobd


  !----------------------------------------------------------------------
  ! Advection of potential vorticity and dissipation on gaussian grid
  ! Adapted version of JACOBD for testing the tangent linear equations
  !----------------------------------------------------------------------
  function jacobp(this, psiloc, pvor) result(sjacob)

    class(qg_tl_type), intent(inout) :: this
    real(r8kind),      intent(   in) :: psiloc(this%nsh2)
    real(r8kind),      intent(   in) :: pvor(this%nsh2)

    real(r8kind) :: sjacob(this%nsh2)

    real(r8kind) :: vv(this%nsh2)
    real(r8kind) :: ininag(this%nlat, this%nlon)
    real(r8kind) :: dpsidl(this%nlat, this%nlon)
    real(r8kind) :: dvordl(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon)
    real(r8kind) :: dpsidm(this%nlat, this%nlon)

    integer      :: i, j, k

    ! Isotropic dissipation term
    do k = 1, this%nsh2
      vv(k) = this%diss(k, 2) * psiloc(k)
    end do

    ininag = this%ggsp%sptogg_pp(vv)

    ! Space derivatives of potential vorticity
    vv = reshape(this%ggsp%ddl(pvor), (/this%nsh2/))
    dvordl = this%ggsp%sptogg_pp(vv)
    dvordm = this%ggsp%sptogg_pd(pvor)

    ! Space derivatives of streamfunction
    vv = reshape(this%ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidl = this%ggsp%sptogg_pp(vv)
    dpsidm = this%ggsp%sptogg_pd(psiloc)

    ! Jacobian and unisotropic dissipation terms
    do i = 1, this%nlat
      do j = 1, this%nlon
        ininag(i, j) = -dpsidm(i, j) * this%ddisdy(i, j) - dpsidl(i, j) * this%ddisdx(i, j) + this%rdiss(i, j) * ininag(i, j)
        ininag(i, j) =  dpsidm(i, j) * dvordl(i, j)      - dpsidl(i, j) * dvordm(i, j)      + ininag(i, j)
      end do
    end do

    sjacob = reshape(this%ggsp%ggtosp(ininag), (/this%nsh2/))

    ! Planetary vorticity advection
    do k = 1, this%nsh2
      sjacob(k) = sjacob(k) - vv(k)
    end do

  end function


  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  subroutine qtopsi(this, qprime, psi, psit)

    class(qg_tl_type),    intent( in) :: this
    real(r8kind),         intent( in) :: qprime(:,:) ! potential vorticity
    real(r8kind),         intent(out) :: psi(:,:)    ! stream function at the nvl levels
    real(r8kind),         intent(out) :: psit(:,:)   ! thickness at the ntl levels

    integer :: k
    real(r8kind) :: r3
    real(r8kind) :: ws(this%nsh2)       ! only used as portable workspace

    do k = 1, size(psi,1)
      ws(k) = qprime(k, 1) + qprime(k, 3)
      psi(k, 1) = this%rinhel(k, 1) * (ws(k) + qprime(k, 2))
      psi(k, 2) = ws(k) - 2.d0 * qprime(k, 2)
      psi(k, 3) = qprime(k, 1) - qprime(k, 3)
    enddo

    do k = 1, size(psit,1)
      psit(k, 1) = this%rinhel(k, 2) * psi(k, 2) + this%rinhel(k, 3) * psi(k, 3)
      psit(k, 2) = this%rinhel(k, 4) * psi(k, 2) + this%rinhel(k, 5) * psi(k, 3)
    enddo

    r3 = 1. / 3.
    do k = 1, size(psi,1)
      psi(k, 2) = r3 * (psi(k, 1) - psit(k, 1) + psit(k, 2))
      psi(k, 1) = psi(k, 2) + psit(k, 1)
      psi(k, 3) = psi(k, 2) - psit(k, 2)
    enddo

    return

  end subroutine qtopsi


  !-----------------------------------------------------------------------
  ! computation of potential vorticity from stream function
  ! input psi streamfunction
  ! output qprime,  the potential vorticity and psit,  the layer thick.
  !-----------------------------------------------------------------------
  subroutine psitoq(this, psi, psit, qprime)
      
    class(qg_tl_type), intent( in) :: this
    real(r8kind),         intent( in) :: psi(:,:)    ! stream function at the nvl levels
    real(r8kind),         intent(out) :: psit(:,:)   ! thickness at the ntl levels
    real(r8kind),         intent(out) :: qprime(:,:) ! potential vorticity

    integer :: k

    do k = 1, size(psit,1)
      psit(k, 1) = psi(k, 1) - psi(k, 2)
      psit(k, 2) = psi(k, 2) - psi(k, 3)
      qprime(k, 1) = this%rinhel(k, 0) * psi(k, 1) - this%rl1 * psit(k, 1)
      qprime(k, 2) = this%rinhel(k, 0) * psi(k, 2) + this%rl1 * psit(k, 1) - this%rl2 * psit(k, 2)
      qprime(k, 3) = this%rinhel(k, 0) * psi(k, 3) + this%rl2 * psit(k, 2)
    enddo

    return

  end subroutine psitoq


  !-----------------------------------------------------------------------
  ! Computation of potential vorticity qout from stream function sfin
  !-----------------------------------------------------------------------
  pure function psiq(this, sfin) result(qout)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),      intent( in) :: sfin(this%nsh2, this%nvl)

    real(r8kind) :: qout(this%nsh2, this%nvl)

    real(r8kind) :: tus(this%nsh2)
    integer      :: k

    do k = 1, this%nsh2
      tus(k) = this%rl1 * sfin(k, 1) - this%rl1 * sfin(k, 2)
    enddo

    do k = 1, this%nsh2
      qout(k, 1) = this%rinhel(k, 0) * sfin(k, 1) - tus(k)
      qout(k, 2) = this%rinhel(k, 0) * sfin(k, 2) + tus(k)
    enddo

    do k = 1, this%nsh2
      tus(k) = this%rl2 * sfin(k, 2) - this%rl2 * sfin(k, 3)
    enddo

    do k = 1, this%nsh2
      qout(k, 2) = qout(k, 2) - tus(k)
      qout(k, 3) = this%rinhel(k, 0) * sfin(k, 3) + tus(k)
    enddo

  end function psiq


  !-----------------------------------------------------------------------
  ! Computation of streamfunction bb from potential vorticity qin
  !-----------------------------------------------------------------------
  pure function qpsi(this, qin) result(sfout)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),      intent( in) :: qin(this%nsh2, this%nvl)

    real(r8kind) :: sfout(this%nsh2, this%nvl)

    real(r8kind) :: ws(this%nsh2)
    real(r8kind) :: tus(this%nsh2, this%ntl), r3
    integer      :: k

    do k = 1, this%nsh2
      ws(k) = qin(k, 1) + qin(k, 3)
      sfout(k, 1) = this%rinhel(k, 1) * (ws(k) + qin(k, 2))
      sfout(k, 2) = ws(k) - 2. * qin(k, 2)
      sfout(k, 3) = qin(k, 1) - qin(k, 3)
    enddo

    do k = 1, this%nsh2
      tus(k, 1) = this%rinhel(k, 2) * sfout(k, 2) + this%rinhel(k, 3) * sfout(k, 3)
      tus(k, 2) = this%rinhel(k, 4) * sfout(k, 2) + this%rinhel(k, 5) * sfout(k, 3)
    enddo

    r3 = 1. / 3
    do k = 1, this%nsh2
      sfout(k, 2) =r3 * (sfout(k, 1) - tus(k, 1) + tus(k, 2))
      sfout(k, 1) = sfout(k, 2) + tus(k, 1)
      sfout(k, 3) = sfout(k, 2) - tus(k, 2)
    enddo
    
  end function qpsi


  !-----------------------------------------------------------------------
  ! transforms francos format to the french format for global fields
  ! input  y spectral coefficients in francos format
  ! output z spectral coefficients in french format
  ! fm format:
  ! k       m  n
  ! 1       0  0
  ! 2       0  1
  ! 3       0  2
  ! :       :  :
  ! nm+1    0  nm
  ! nm+2    1  1 --> real part
  ! nm+3    1  2 --> real part
  ! :       :  :
  ! nm+nm+1 1  nm --> real part
  ! :       :  :
  ! :       nm nm --> real part
  !  repeat for imaginary part
  !  disadvantage: 0 0 mode and imaginary parts of m = 0 modes are obsolete
  ! fs format stores all m for every n first and has no obsolete indices
  ! 
  ! k       m  n
  ! 1       0  1
  ! 2       1  1 --> real part
  ! 3       1  1 --> imaginary part: k = 1-3 is T1 truncation
  ! 4       0  2
  ! 5       1  2 --> real part
  ! 6       1  2 --> imaginary part
  ! 7       2  2 --> real part
  ! 8       2  2 --> imaginary part: k = 1-8 is T2 truncation
  ! etcetera
  !-----------------------------------------------------------------------
  pure function fmtofs (this, y) result(z)

    class(qg_tl_type), intent(in) :: this
    real(r8kind), intent( in)        :: y(:,:)

    real(r8kind), dimension(size(y,1),size(y,2)) :: z

    integer ::  m, n, k, indx, l

    do l = 1, size(y,2)
      k = 1
      do m = 0, this%nm
        do n = max(m, 1), this%nm
          k = k + 1
          if (m .eq. 0) then
            indx = n**2
          else
            indx = n**2 + 2 * m - 1
          end if
          z(indx, l) = y(k, l)
          if (m .ne. 0) z(indx + 1, l) = y(k + this%nsh, l)
        enddo
      enddo
    enddo

    return

  end function fmtofs


  !-----------------------------------------------------------------------
  ! transforms the french format to francos format for global fields
  ! input  y spectral coef. in french format,  ntr is truncation limit
  ! output z spectral coefficients in francos format
  ! fm format:
  ! k       m  n
  ! 1       0  0
  ! 2       0  1
  ! 3       0  2
  ! :       :  :
  ! nm+1    0  nm
  ! nm+2    1  1 --> real part
  ! nm+3    1  2 --> real part
  ! :       :  :
  ! nm+nm+1 1  nm --> real part
  ! :       :  :
  ! :       nm nm --> real part
  !  repeat for imaginary part
  !  disadvantage: 0 0 mode and imaginary parts of m = 0 modes are obsolete
  ! fs format stores all m for every n first and has no obsolete indices
  ! 
  ! k       m  n
  ! 1       0  1
  ! 2       1  1 --> real part
  ! 3       1  1 --> imaginary part: k = 1-3 is T1 truncation
  ! 4       0  2
  ! 5       1  2 --> real part
  ! 6       1  2 --> imaginary part
  ! 7       2  2 --> real part
  ! 8       2  2 --> imaginary part: k = 1-8 is T2 truncation
  ! etcetera
  !-----------------------------------------------------------------------
  pure function fstofm (this, y, ntr) result(z)

    class(qg_tl_type), intent(in) :: this
    real(r8kind),         intent(in) :: y(:,:)
    integer,              intent(in) :: ntr

    real(r8kind), dimension(size(y,1),size(y,2)) :: z

    integer :: m, n, k, indx, i, l

    do l = 1, size(y,2)
      do i = 1, size(y,1)
        z(i, l) = 0d0
      enddo
      k = 1
      do m = 0, this%nm
        do n = max(m, 1), this%nm
          k = k + 1
          if ((m .le. ntr).and.(n .le. ntr)) then
            if (m .eq. 0) then
              indx = n**2
            else
              indx = n**2 + 2 * m - 1
            end if
            z(k, l) = y(indx, l)
            if (m .ne. 0) z(k + this%nsh, l) = y(indx + 1, l)
          endif
        enddo
      enddo
    enddo

    return

  end function fstofm


  !-----------------------------------------------------------------------
  ! computation of geostrophic winds at all levels
  ! computes geopotential height in [m2 / s2[ = [gz] from streamfunction 
  ! by solving the linear balance equation: 
  ! del phi = (1 - mu**2 ) d psi / dmu + mu del psi
  ! the global mean value is not determined and set to zero
  !-----------------------------------------------------------------------
  subroutine gridfields(this, lat, lon, lvl, geopg, psig, forg, qgpv, ug, vg)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),         intent(out) :: lat(:)      ! Grid latitude
    real(r8kind),         intent(out) :: lon(:)      ! Grid longitude
    real(r8kind),         intent(out) :: lvl(:)      ! Grid level
    real(r8kind),         intent(out) :: geopg(:,:,:)! Geopotential on the grid
    real(r8kind),         intent(out) :: psig(:,:,:) ! Grid values of dimensional streamfunction at the three levels
    real(r8kind),         intent(out) :: forg(:,:,:) ! Grid values of dimensional forcing at the three levels
    real(r8kind),         intent(out) :: qgpv(:,:,:) ! Grid values of dimensional pv at the three levels
    real(r8kind),         intent(out) :: ug(:,:,:)   ! Grid values of zonal velocity at the three levels in m/s
    real(r8kind),         intent(out) :: vg(:,:,:)   ! Grid values of meridional velocity at the three levels in m/s

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)

    integer :: i, j, k, l
    real(r8kind) :: dlon
    real(r8kind) :: facwind, facsf, facgp, facpv
    real(r8kind) :: dpsdl(this%nlat, this%nlon), dpsdm(this%nlat, this%nlon), psik(this%nsh2), vv(this%nsh2)
    real(r8kind) :: fmu(this%nlat)
    real(r8kind) :: delpsis(this%nsh2), delpsig(this%nlat, this%nlon)
    real(r8kind) :: dmupsig(this%nlat, this%nlon), delgeog(this%nlat, this%nlon)
    real(r8kind) :: delgeos(this%nsh2), geos(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! Get longitude increment
    dlon = 360d0 / real(this%nlon)

    ! space derivatives of streamfunction
    facwind = radius * om
    facsf = om * (radius)**2
    facgp = (om**2) * (radius**2)
    facpv = om

    do i = 1, this%nlat
      fmu(i) = 1 - sin(pi * this%phi(i) / 180d0)**2
    enddo

    ! Calculate in the latitude and longitude fields
    do i = 1, this%nlat
      lat(i) = this%phi(i)
    end do
    do j = 1, this%nlon
      lon(j) = (j - 1) * dlon
    end do

    ! Calculate the level field
    lvl = (/200.0, 500.0, 800.0/)

    ! Calculate remaining fields
    do l = 1, this%nvl

      psig(:, :, l) = ggsp%sptogg_pp(this%psi(:, l))
      forg(:,:,l) = ggsp%sptogg_pp(this%for(:, l))
      qgpv(:, :, l) = ggsp%sptogg_pp(this%qprime(:, l))

      do k = 1, this%nsh2
        psik(k) = this%psi(k, l)
      enddo

      vv = reshape(ggsp%ddl (psik), (/this%nsh2/))
      dpsdl = ggsp%sptogg_pp (vv)
      dpsdm = ggsp%sptogg_pd (psik)

      ! solve linear balance equation
      delpsis = this%lap(this%psi(:, l))
      delpsig = ggsp%sptogg_pp(delpsis)
      dmupsig = ggsp%sptogg_pd(this%psi(:, l))

      do j = 1, this%nlon
        do i = 1, this%nlat
          psig(i, j, l) = facsf * psig(i, j, l)
          qgpv(i, j, l) = facpv * qgpv(i, j, l)
          ug(i, j, l) = -facwind * dpsdm(i, j) * this%cosfi(i)
          vg(i, j, l) = +facwind * dpsdl(i, j) / this%cosfi(i)
          delgeog(i, j) = fmu(i) * dmupsig(i, j) + this%sinfi(i) * delpsig(i, j)
        enddo
      enddo

      delgeos = reshape(ggsp%ggtosp(delgeog), (/this%nsh2/))
      geos = this%lapinv(delgeos)
      geos(1) = 0.d0
      geopg(:, :, l) = ggsp%sptogg_pp(geos)

      do j = 1, this%nlon
        do i = 1, this%nlat
          geopg(i, j, l) = facgp * geopg(i, j, l)
        enddo
      enddo

    enddo

  end subroutine gridfields

      
  !-----------------------------------------------------------------------
  ! Return the lat/lon grid values
  !-----------------------------------------------------------------------
  function get_lat_lon_grid(this) result(grid)

    class(qg_tl_type), intent( in) :: this
    real(r8kind), dimension(3) :: grid(this%nlon, this%nlat, this%nvl, 3)

    integer      :: i, j, l
    real(r8kind) :: dlon, lat, lvl
    real(r8kind) :: levels(this%nvl)      ! Grid level

    ! Set the levels
    levels = (/200.0, 500.0, 800.0/)

    ! Get longitude increment
    dlon = 360d0 / real(this%nlon)

    ! Calculate the grid locations
    do l = 1, this%nvl
      lvl = levels(l)
      do i = 1, this%nlat
        lat = this%phi(i)
        do j = 1, this%nlon
          grid(j,i,l,:) = (/(j - 1) * dlon, lat, lvl/)
        end do
      end do
    enddo

  end function get_lat_lon_grid


  !-----------------------------------------------------------------------
  ! Return the streamfunction on the gaussian grid
  !-----------------------------------------------------------------------
  function get_psig(this) result(psig)

    class(qg_tl_type), intent( in) :: this
    real(r8kind), dimension(this%nlon, this%nlat, this%nvl) :: psig ! Grid values of dimensional streamfunction at the three levels

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)

    integer :: i, j, l
    real(r8kind) :: facsf
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of streamfunction
    facsf = om * (radius)**2

    ! Calculate remaining fields
    do l = 1,this%nvl
      psig(:, :, l) = transpose(ggsp%sptogg_pp(this%psi(:, l)))
      do j = 1, this%nlat
        do i = 1, this%nlon
          psig(i, j, l) = facsf * psig(i, j, l)
        enddo
      enddo

    enddo

  end function get_psig


  !-----------------------------------------------------------------------
  ! computation of laplace operator in spectral domain
  ! input  xs  field in spectral form
  ! output xsl laplace of xs in spectral form
  !-----------------------------------------------------------------------
  function lap(this, xs) result(xsl)

    class(qg_tl_type),  intent(in) :: this
    real(r8kind),          intent(in) :: xs(:)
    real(r8kind), dimension(size(xs)) :: xsl

    integer :: k

    do k = 1, size(xs)
      xsl(k) = xs(k) * this%rinhel(k, 0)
    enddo

    return

  end function lap

      
  !-----------------------------------------------------------------------
  ! computation of laplace operator in spectral domain
  ! input  xsl field in spectral form
  ! output xs  inverse laplace of xs in spectral form
  !-----------------------------------------------------------------------
  function lapinv(this, xsl) result(xs)

    class(qg_tl_type),   intent(in) :: this
    real(r8kind),           intent(in) :: xsl(:)
    real(r8kind), dimension(size(xsl)) :: xs

    integer :: k

    do k = 1, size(xsl)
      xs(k) = xsl(k) * this%rinhel(k, 1)
    enddo

    return

  end function lapinv


  !-----------------------------------------------------------------------
  ! computation of artifical forcing according to roads(1987)
  ! the forcing is computed from file obsfile
  !------------------------------------------------------------------------
  function artiforc(this) result(for)

    class(qg_tl_type), intent(in) :: this

    integer :: i, j, k, l, iday, nvar
    real(r4kind) :: psi4(this%nsh2, 3)
    real(r8kind) :: sum(this%nsh2, 3), forg(this%nlat, this%nlon), dlon, psifs(this%nsh2, 3)
    real(r8kind) :: psi(this%nsh2,this%nvl)    ! stream function at the nvl levels
    real(r8kind) :: psit(this%nsh2,this%ntl)   ! thickness at the ntl levels
    real(r8kind) :: qprime(this%nsh2,this%nvl) ! potential vorticity
    real(r8kind) :: for(this%nsh2,this%nvl)    ! constant potential vorticity forcing at the nvl levels
    real(r8kind) :: dqprdt(this%nsh2,this%nvl) ! time derivative of qprime

    ! Only used in artiforc
    character(len=32)  :: obsfile  ! Name of observation file
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    obsfile = this%config%get_obsfile()

    nvar = (this%nm + 2) * this%nm
    dlon = 360d0 / real(this%nlon)

    write(*, '(A, A)') "Calculating forcing from ", obsfile

    do l = 1, this%nvl
      do k = 1, this%nsh2
        sum(k, l) = 0.
        for(k, l) = 0d0
      enddo
    enddo

    ! calculate the mean tendency
    open(unit = 46, file = './' // obsfile, status = 'old', form = 'unformatted')
    iday = 0
 10 continue
    do l = 1, this%nvl
      read(46, end = 20) (psi4(k, this%nvl - l + 1), k = 1, nvar)
    enddo

    iday = iday + 1
    do l = this%nvl, 1, -1
      do k = 1, nvar
        psifs(k, l) = psi4(k, l)
      enddo
    enddo
    psi = this%fstofm(psifs, this%nm)
    call this%psitoq(psi, psit, qprime)
    dqprdt = this%ddt(psi, psit, qprime, for)
    do l = 1, this%nvl
      do k = 1, this%nsh2
        sum(k, l) = sum(k, l) + dqprdt(k, l)
      enddo
    enddo
    goto 10
 20 continue
    close(46)

    do l = 1, this%nvl
      do k = 1, this%nsh2
        for(k, l) = -sum(k, l) / real(iday)
      enddo
    enddo

    write(*, '(A, I6)') "Number of states used to calculate forcing: ", iday

  end function artiforc


  !------------------------------------------------------------------
  ! timinp
  !
  ! VERSION OF TIMINT WITH ARGUMENTS
  !
  ! *** JTSTEP = 0 : Forward integration, half time step
  ! *** JTSTEP = 1 : Leapfrog integration, half time step
  ! *** JTSTEP > 1 : Leapfrog integration + implicit diffusion
  ! ***              and time filter, full time step
  !
  !------------------------------------------------------------------
  subroutine timinp (this, jtstep, xold, xnew, dxnew)

    class(qg_tl_type), intent( in) :: this
    integer,           intent( in) :: jtstep
    real(r8kind),      intent(out) :: xold(this%nsh2, 3), xnew(this%nsh2, 3)
    real(r8kind),      intent( in) :: dxnew(this%nsh2, 3)

    integer      :: jstep
    integer      :: i, j, k, l, m, n
    real(r8kind) :: as(this%nsh2, 3)
    real(r8kind) :: ws(this%nsh2)
    real(r8kind) :: eps1, eps2

    eps2 = 0.03
    eps1 = 1.0D0-2.0D0 * eps2 

    if (jtstep > 1) then
      do l = 1, 3
        do k = 1, this%nsh2
          ws(k) = xnew(k, l)
          xnew(k, l) = this%diss(k, 1) * (xold(k, l) + this%dtt * dxnew(k, l))
          xold(k, l) = eps1 * ws(k) + eps2 * (xnew(k, l) + xold(k, l))
        end do
      end do
    else
      if (jtstep==0) then
        do l = 1, 3
          do k = 1, this%nsh2
            xold(k, l) = xnew(k, l)
          end do
        end do
      endif
      do l = 1, 3
        do k = 1, this%nsh2
          xnew(k, l) = xold(k, l) + this%dtt * dxnew(k, l)
        end do
      end do

!     this%dtt = 2. * this%dtt

    endif

  end subroutine timinp

 
  !------------------------------------------------------------------
  ! COMPUTES THE TANGENT LINEAR EQUATIONS
  !------------------------------------------------------------------
  subroutine ddtl(this, ds, depsdt, jstep)

    class(qg_tl_type), intent(inout) :: this
    real(r8kind),      intent(   in) :: ds(this%nsh2, 3)
    real(r8kind),      intent(  out) :: depsdt(this%nsh2, 3)
    integer,           intent(   in) :: jstep

    integer      :: nnnstep
    real(r8kind) :: cp(this%nsh2), wspace(this%nsh2), as(this%nsh2, 3), bs(this%nsh2, 3)
    real(r8kind) :: dummy, dfeps(this%nsh2, 3), delpst(this%nsh2, 2)
    real(r8kind) :: delpsi(this%nsh2, 3), lpsi(this%nsh2, 3)

    integer i, j, k, l, m, n

    real(r8kind) :: xs(this%nsh2, 3), fac
!    integer      :: nstop, janstep
!    common /ntime/ nstop

    ! *** Initializations needed for the computation of c1, c2 and c3
    do l = 1, 3
      do k = 1, this%nsh2
        dfeps(k, l) = 0.0D0
      end do
    end do
    do k = 1, this%nsh2
      cp(k) = 0.0D0
    end do
    cp(2) = 1.0D0 / sqrt(3.0D0)

    call copy(ds, delpsi)

    ! MATRIX L TOEPASSEN:
    as = this%psiq(delpsi)
    lpsi = this%psiq(this%psi)
    do k = 1, this%nsh2
      delpst(k, 1) = this%relt1 * (delpsi(k, 1) - delpsi(k, 2))
      delpst(k, 2) = this%relt2 * (delpsi(k, 2) - delpsi(k, 3))
    end do

    call sub (as(:, 1), cp, bs(:, 1))
    wspace(:) = this%jacob(this%psi(:,1), bs(:,1))
    call addab (dfeps(:, 1), wspace(:))

    bs(:,1) = this%jacob(delpsi(:,1), lpsi(:,1))
    call addab (dfeps(:, 1), bs(:, 1))

    call addab (dfeps(:, 1), delpst(:, 1))

    call sub (as(:, 2), cp, bs(:, 2))
    wspace(:) = this%jacob(this%psi(:, 2), bs(:, 2))
    call addab (dfeps(:, 2), wspace(:))

    bs(:,2) = this%jacob(delpsi(:, 2), lpsi(:, 2))
    call addab (dfeps(:, 2), bs(:, 2))

    call subab (dfeps(:, 2), delpst(:, 1))
    call addab (dfeps(:, 2), delpst(:, 2))

    call sub (as(:, 3), cp, bs(:, 3))  
    wspace(:) = this%jacob(this%psi(:, 3), bs(:, 3))
    call addab (dfeps(:, 3), wspace(:))

    call add(lpsi(:, 3), this%orog, wspace(:))
    bs(:,3) = this%jacobp(delpsi(:, 3), wspace(:))
    call addab (dfeps(:, 3), bs(:, 3))

    call subab (dfeps(:, 3), delpst(:, 2))

    depsdt = this%qpsi(dfeps)

  end subroutine ddtl


  !-------------------------------------------------------------------------------
  ! This subroutine integrates the tangent linear equation with initial condition
  ! AS(NSH2,3) and integration time NUUR. It uses the orbit stored by STOREP.
  ! The result is BS(NSH2,3).
  !-------------------------------------------------------------------------------
!  subroutine tang(this, as, bs, nuur)
  function tang(this, as, nuur) result(bs)

    class(qg_tl_type), intent(inout) :: this
    real(r8kind),      intent(   in) :: as(this%nsh2, 3)
    integer,           intent(   in) :: nuur

    real(r8kind) :: bs(this%nsh2, 3)
    real(r8kind) :: dum(this%nsh23), dumold(this%nsh23), dumnew(this%nsh23)
    real(r8kind) :: etraj(this%nsh2, 3, 0:100)
!    real(r8kind) :: e, fac

!    integer      :: janstep
    integer      :: k, n, index
!    integer      :: nstop
!    common /stap/ janstep
!    common /ntime/ nstop

    real(r8kind) :: zt(this%nlat, this%nlon), dummy, sum
    integer      :: ix, iy, ilevel
    integer      :: ntstep, jtstep, i

    ntstep = (nuur * 3) / 2
!    nstop=ntstep

!    dumnew = reshape(this%trajectory, (/this%nsh23/))
    dumnew = reshape(as, (/this%nsh23/))

    do i = 1, this%nsh23
      dumold(i) = 0.0D0
    end do

    write(*,*) 'TANGENT MODEL'
    print *, "IN TANGENT, DTT=", this%dtt, ntstep

    do jtstep = 0, ntstep
      ! call retrie(jtstep)
!CWH      call retrie(0)
      
      call this%ddtl(dumnew, dum, jtstep)
      call this%timinp(jtstep, dumold, dumnew, dum)

      if (0 == 1) then
        k = jtstep
        if( mod(k,9) == 0) then
          write(*,*) 'write-out linear evo', k
          if (k == 0) then
            write(66) as
!            write(66) this%trajectory
          else
            write(66) dumnew
          endif
        endif
      endif

      if (0 == 1) then
        k = jtstep
        if (k == 0) then
          write(*,*) 'write LIN TRAJ to 51 !!!'
        endif
        if (mod(k, 9) == 0) then
          k = (jtstep / 3) * 2
          write(66) dumnew
        endif
      endif

      if (0 == 1) then 
        ! fill array ETRAJ to determine costfunction
        k = jtstep
        if (0 == 1) then
          call copy(reshape(dumnew, (/this%nsh2, 3/)), etraj(:, :, index))

          ! one gridpoint (vorticity)  --------------------------------
          if (0 == 1) then
            etraj(:, 3, index) = this%lap(etraj(:, 3, index))
            zt = this%ggsp%sptogg_pp(etraj(:, 3, index))
            dummy = zt(24, 30)
            write(*,*) 'etraj(index)', dummy
            if (jtstep == 0) sum = 0
            sum = sum + dummy**2
            if (jtstep == ntstep) write(*,*) 'costf = ', sum
            do iy = 1, this%nlat
              do ix = 1, this%nlon
                zt(iy, ix) = 0.0
              enddo
            enddo
            zt(24, 30) = dummy
            etraj(:, 3, index) = reshape(this%ggsp%ggtosp(zt), (/this%nsh2/))
            etraj(:, 3, index) = this%lap(etraj(:, 3, index))
            do ilevel = 1, 2
              do n = 1, this%nsh2
                etraj(n, ilevel, index) = 0.0
              enddo
            enddo
          endif
        endif
      endif

    end do

!    call copy(reshape(dumnew, (/this%nsh2, 3/)), this%trajectory)
    call copy(reshape(dumnew, (/this%nsh2, 3/)), bs)

  end function


  !-------------------------------------------------------------------------------
  ! tang_reverse
  !-------------------------------------------------------------------------------
!  subroutine tang_reverse(this, as, bs, nuur)
  subroutine tang_reverse(this, nuur)

    class(qg_tl_type), intent(inout) :: this
!    real(r8kind),      intent(   in) :: as(this%nsh2, 3)
!    real(r8kind),      intent(  out) :: bs(this%nsh2, 3)
    integer,           intent(   in) :: nuur
    
!    this%dtt = -1.0D0 * this%dtt

!    call this%tang(as, bs, nuur)
!    call this%tang(nuur)

!    this%dtt = -1.0D0 * this%dtt

  end subroutine tang_reverse


  !-------------------------------------------------------------------------------
  ! get_config
  !-------------------------------------------------------------------------------
  function get_config(this) result(config)

    class(qg_tl_type), intent(in) :: this
    type(qg_config_type)          :: config

    config = this%config

  end function get_config


  !-------------------------------------------------------------------------------
  ! get_step
  !-------------------------------------------------------------------------------
  function get_step(this) result(step)

    class(qg_tl_type), intent(in) :: this
    integer                          :: step

    step = this%step

  end function get_step


  !-------------------------------------------------------------------------------
  ! get_clock
  !-------------------------------------------------------------------------------
  function get_clock(this) result(clock)

    class(qg_tl_type), intent(in) :: this
    real(r8kind)                     :: clock

    clock = this%clock

  end function get_clock


  !-------------------------------------------------------------------------------
  ! get_nlat
  !-------------------------------------------------------------------------------
  function get_nlat(this) result(nlat)

    class(qg_tl_type), intent(in) :: this
    integer                          :: nlat

    nlat = this%nlat

  end function get_nlat


  !-------------------------------------------------------------------------------
  ! get_nlon
  !-------------------------------------------------------------------------------
  function get_nlon(this) result(nlon)

    class(qg_tl_type), intent(in) :: this
    integer                          :: nlon

    nlon = this%nlon

  end function get_nlon


  !-------------------------------------------------------------------------------
  ! get_nsh2
  !-------------------------------------------------------------------------------
  function get_nsh2(this) result(nsh2)

    class(qg_tl_type), intent(in) :: this
    integer                          :: nsh2

    nsh2 = this%nsh2

  end function get_nsh2

  !-------------------------------------------------------------------------------
  ! get_nvl
  !-------------------------------------------------------------------------------
  function get_nvl(this) result(nvl)

    class(qg_tl_type), intent(in) :: this
    integer                          :: nvl

    nvl = this%nvl

  end function get_nvl


  !-------------------------------------------------------------------------------
  ! get_psi
  !-------------------------------------------------------------------------------
  function get_psi(this) result(psi)

    class(qg_tl_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: psi

    psi = this%psi

  end function get_psi


  !-------------------------------------------------------------------------------
  ! get_trajectory
  !-------------------------------------------------------------------------------
  function get_trajectory(this) result(trajectory)

    class(qg_tl_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: trajectory

    trajectory = this%trajectory

  end function get_trajectory


  !-------------------------------------------------------------------------------
  ! get_qprimed
  !-------------------------------------------------------------------------------
  function get_qprimed(this) result(qprimed)

    class(qg_tl_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: qprimed

    qprimed = this%qprimed

  end function get_qprimed


  !-------------------------------------------------------------------------------
  ! get_for
  !-------------------------------------------------------------------------------
  function get_for(this) result(for)

    class(qg_tl_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: for

    for = this%for

  end function get_for


  !-------------------------------------------------------------------------------
  ! get_state_vector
  !-------------------------------------------------------------------------------
  function get_state_vector(this) result(state_vector)

    class(qg_tl_type),            intent(in) :: this
    real(r8kind), dimension(this%nlat * this%nlon * this%nvl) :: state_vector

    state_vector = reshape(this%get_psig(),(/this%nlat * this%nlon * this%nvl/))

  end function get_state_vector


  !-------------------------------------------------------------------------------
  ! get_location_vector
  !-------------------------------------------------------------------------------
  function get_location_vector(this) result(location_vector)

    class(qg_tl_type),            intent(in) :: this
    real(r8kind), dimension(this%nlat * this%nlon * this%nvl, 3) :: location_vector

    location_vector = reshape(this%get_lat_lon_grid(),(/this%nlat * this%nlon * this%nvl, 3/))

  end function get_location_vector


  !------------------------------------------------------------------
  ! get_interpolation_weights
  !------------------------------------------------------------------
  subroutine get_interpolation_weights(this, latx, lonx, lvl, NW_index, NE_index, SW_index, SE_index, NW_weight, NE_weight, SW_weight, SE_weight)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),         intent( in) :: latx, lonx, lvl
    integer,              intent(out) :: NW_index, NE_index, SW_index, SE_index
    real(r8kind),         intent(out) :: NW_weight, NE_weight, SW_weight, SE_weight

    real(r8kind) :: lat, lon, lat_north, lat_south, lon_west, lon_east, sum_weight, levels(this%nvl)
    integer      :: ilat_north, ilat_south, ilon_west, ilon_east, ilvl

    ! Get the bounding latitudes
    ilat_south = 1
    do while (this%phi(ilat_south) < latx)
      ilat_south = ilat_south + 1
    end do
    ilat_south = ilat_south - 1
    ilat_north = ilat_south + 1
    lat_north = this%phi(ilat_north)
    lat_south = this%phi(ilat_south)

    ! Get the bounding longitudes
    ilon_west = int(lonx / (360d0 / real(this%nlon))) + 1
    ilon_east = ilon_west + 1
    if (ilon_east > this%nlon) ilon_east = 1
    lon_west = (ilon_west - 1) * 360d0 / real(this%nlon)
    lon_east = (ilon_east - 1) * 360d0 / real(this%nlon)

    ! Get the level index
    ilvl = 1
    levels = (/200.0, 500.0, 800.0/)
    do while (levels(ilvl) < lvl)
      ilvl = ilvl + 1
    end do

    ! Compute the indices of the bounding box after it is mapped from 3D array into a 1D vector
    NW_index = this%nlat * this%nlon * (ilvl - 1) + this%nlon * (ilat_north - 1) + ilon_west
    NE_index = this%nlat * this%nlon * (ilvl - 1) + this%nlon * (ilat_north - 1) + ilon_east
    SW_index = this%nlat * this%nlon * (ilvl - 1) + this%nlon * (ilat_south - 1) + ilon_west
    SE_index = this%nlat * this%nlon * (ilvl - 1) + this%nlon * (ilat_south - 1) + ilon_east

    ! Compute the distances to bounding box vertices
    NW_weight = 1.0 / distance(latx, lonx, lat_north, lon_west)
    NE_weight = 1.0 / distance(latx, lonx, lat_north, lon_east)
    SW_weight = 1.0 / distance(latx, lonx, lat_south, lon_west)
    SE_weight = 1.0 / distance(latx, lonx, lat_south, lon_east)

    ! Compute weights by normalizing distances
    sum_weight = NW_weight + NE_weight + SW_weight + SE_weight
    NW_weight = NW_weight / sum_weight
    NE_weight = NE_weight / sum_weight
    SW_weight = SW_weight / sum_weight
    SE_weight = SE_weight / sum_weight

  end subroutine get_interpolation_weights


  function distance(lat, lon, latx, lonx) result(d)

    real(r8kind), intent(in) :: lat, lon, latx, lonx
    real(r8kind)             :: d

    real(r8kind) :: rlat, rlatx, dlat, dlon
    real(r8kind) :: a, c

    real(r8kind), parameter :: dtor = atan(1.0) / 45.0
    real(r8kind), parameter :: erad = 6372.8

    rlat = dtor * lat
    rlatx = dtor * latx

    dlat = dtor * (latx - lat)
    dlon = dtor * (lonx - lon)

    a = (sin(dlat / 2.0))**2 + cos(rlat) * cos(rlatx) * (sin(dlon / 2.0))**2
    c = 2.0 * asin(sqrt(a))
    d = erad * c

  end function distance

end module QG_Model_TL
