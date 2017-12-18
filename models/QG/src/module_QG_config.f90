module QG_Config

  use kind,            only : r4kind, r8kind
  use Abstract_Config, only : abstract_config_type
  use QG_GGSP

  implicit none

  private

  public :: qg_config_type

  type, extends(abstract_config_type) :: qg_config_type
    private
    ! Namelist params variables
    integer           :: resolution ! Model resolution
    integer           :: time_step  ! Model time step (in seconds)
    character(len=32) :: obsfile    ! Name of observation file
    logical           :: inf        ! If .true. then artificial PV forcing read from file
    logical           :: obsf       ! If .true. PV forcing is calculated from observations in routine artiforc
    real(r8kind)      :: tdis       ! Ekman dissipation timescale in days at lower level
    real(r8kind)      :: addisl     ! Parameter used in the computation of the dissipation timescale at the lower level over land
    real(r8kind)      :: addish     ! Parameter used in the computation of the dissipation timescale at the lower level as a function of topography
    real(r8kind)      :: trel       ! Relaxation time scale in days of the temperature
    real(r8kind)      :: tdif       ! Dissipation timescale of scale-selective diffusion in days for wavenumber nm
    integer           :: idif       ! Determines scale-selectivity of hyperviscosity; power of laplace operator
    real(r8kind)      :: h0         ! scale factor for the topographically induced upward motion at the lower level
    real(r8kind)      :: rrdef1     ! Rossby radius of deformation of 200-500 thickness
    real(r8kind)      :: rrdef2     ! Rossby radius of deformation of 500-800 thickness

    ! Model grid dimensions
    character(len=3) ::  ft   ! Character string containing resolution
    integer :: nm             ! The truncation is of type T(riangular) nm
    integer :: nlon           ! Number of longitude points of the Gaussian grid
    integer :: nlat           ! Number of latitude  points of the Gaussian grid
    integer :: nvl            ! Number of vorticity levels in the vertical (should be set to 3)
    integer :: ntl            ! Number of temperature levels in the vertical (equal to nvl-1)
    integer :: nsh            ! Half of nsh2
    integer :: nsh2           ! Number of coefficients needed to define one level of the T nm model

    ! Spectral coefficients
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
    real(r8kind), allocatable  :: dorodl(:,:) ! Derivative of orog wrt lambda
    real(r8kind), allocatable  :: dorodm(:,:) ! Derivative of orag wrt sin(fi)
    real(r8kind), allocatable  :: rdiss(:,:)  ! Landsea-mask/orography dependent friction
    real(r8kind), allocatable  :: ddisdx(:,:) ! Landsea-mask/orography dependent friction
    real(r8kind), allocatable  :: ddisdy(:,:) ! Landsea-mask/orography dependent friction

    ! Model Forcing
    real(r8kind), allocatable :: for(:,:)    ! Constant potential vorticity forcing at the nvl levels

  contains
    final :: destructor_qg_config
    procedure :: init_grid
    procedure :: init_spectral_coeff
    procedure :: init_laplace_helmholtz
    procedure :: init_orography
    procedure :: init_forcing
    procedure :: artiforc
    procedure :: get_resolution
    procedure :: get_time_step
    procedure :: get_obsfile
    procedure :: get_inf
    procedure :: get_obsf
    procedure :: get_tdis
    procedure :: get_addisl
    procedure :: get_addish
    procedure :: get_trel
    procedure :: get_tdif
    procedure :: get_idif
    procedure :: get_h0
    procedure :: get_rrdef1
    procedure :: get_rrdef2

    procedure :: get_ggsp
    procedure :: get_for
    procedure :: get_nm
    procedure :: get_nsh
    procedure :: get_nsh2
    procedure :: get_nvl
    procedure :: get_ntl
    procedure :: get_nlat
    procedure :: get_nlon
    procedure :: get_ft
    procedure :: get_ll
    procedure :: get_relt1
    procedure :: get_relt2
    procedure :: get_diss
    procedure :: get_phi
    procedure :: get_sinfi
    procedure :: get_cosfi
    procedure :: get_dorodl
    procedure :: get_dorodm
    procedure :: get_lgdiss
    procedure :: get_rdiss
    procedure :: get_ddisdx
    procedure :: get_ddisdy
    procedure :: get_rinhel
    procedure :: get_rl1
    procedure :: get_rl2

  end type qg_config_type

  interface qg_config_type
    procedure constructor_arglist
    procedure constructor_namelist_file
    procedure constructor_namelist_unit
  end interface

  ! Mathematical and physical constants
  real(r8kind), parameter :: pi = 4d0 * atan(1d0)     ! value of pi

contains


  !------------------------------------------------------------------
  ! constructor_arglist
  !
  ! Returns an initialized qg_config_type object
  !------------------------------------------------------------------
  function constructor_arglist(resolution, time_step, obsfile, inf, obsf, tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2, for) result(qg_config)

    integer,                intent(in) :: resolution ! Model resolution
    integer,                intent(in) :: time_step  ! Model time step (in seconds)
    character(len=32),      intent(in) :: obsfile    ! Name of observation file
    logical,                intent(in) :: inf        ! If .true. then artificial PV forcing read from file
    logical,                intent(in) :: obsf       ! If .true. PV forcing is calculated from observations in routine artiforc
    real(r8kind),           intent(in) :: tdis       ! Ekman dissipation timescale in days at lower level
    real(r8kind),           intent(in) :: addisl     ! Parameter used in the computation of the dissipation timescale at the lower level over land
    real(r8kind),           intent(in) :: addish     ! Parameter used in the computation of the dissipation timescale at the lower level as a function of topography
    real(r8kind),           intent(in) :: trel       ! Relaxation time scale in days of the temperature
    real(r8kind),           intent(in) :: tdif       ! Dissipation timescale of scale-selective diffusion in days for wavenumber nm
    integer,                intent(in) :: idif       ! Determines scale-selectivity of hyperviscosity; power of laplace operator
    real(r8kind),           intent(in) :: h0         ! Scale factor for the topographically induced upward motion at the lower level
    real(r8kind),           intent(in) :: rrdef1     ! Rossby radius of deformation of 200-500 thickness
    real(r8kind),           intent(in) :: rrdef2     ! Rossby radius of deformation of 500-800 thickness
    real(r8kind), optional, intent(in) :: for(:,:)
    type(qg_config_type)               :: qg_config

    ! Initialize model parameters    
    qg_config%resolution = resolution
    qg_config%time_step = time_step
    qg_config%obsfile = obsfile
    qg_config%inf = inf
    qg_config%obsf = obsf
    qg_config%tdis = tdis
    qg_config%addisl = addisl
    qg_config%addish = addish
    qg_config%trel = trel
    qg_config%tdif = tdif
    qg_config%idif = idif
    qg_config%h0 = h0
    qg_config%rrdef1 = rrdef1
    qg_config%rrdef2 = rrdef2
   
    ! Initialize the grid
    call qg_config%init_grid(qg_config%resolution)

    ! Read and initialize spectral coefficients
    call qg_config%init_spectral_coeff()

    ! Instantiate a ggsp grid conversion object
    qg_config%ggsp = qg_ggsp_type(qg_config%nm, qg_config%nlat, qg_config%nlon, qg_config%nshm, qg_config%pp, qg_config%pd, qg_config%pw)

    ! Initialize laplace and helmholtz operators
    call qg_config%init_laplace_helmholtz(qg_config%rrdef1, qg_config%rrdef2, qg_config%trel, qg_config%tdis, qg_config%tdif, qg_config%idif)

    ! Initialize orography
    call qg_config%init_orography(qg_config%h0, qg_config%addisl, qg_config%addish)

    ! Initialize model forcing
    if (present(for)) then
      call qg_config%init_forcing(qg_config%obsf, for=for)
    else
      call qg_config%init_forcing(qg_config%obsf)
    end if

  end function constructor_arglist


  !------------------------------------------------------------------
  ! constructor_namelist_file
  !
  ! Returns an initialized qg_config_type object
  !------------------------------------------------------------------
  function constructor_namelist_file(nl_filename) result(qg_config)

    ! Namelist filename
    character(len=*), intent(in) :: nl_filename
    type(qg_config_type)         :: qg_config

    ! Namelist file descriptor
    integer :: nl_unit

    ! Namelist param variables
    integer           :: resolution = 21            ! Model resolution
    integer           :: time_step = 1200           ! Model time step (in seconds)
    character(len=32) :: obsfile ='sf7910T106.shfs' ! Name of observation file
    logical           :: inf = .false.              ! If .true. then artificial PV forcing read from file
    logical           :: obsf = .false.             ! If .true. PV forcing is calculated from observations in routine artiforc
    real(r8kind)      :: tdis = 3.0                 ! Ekman dissipation timescale in days at lower level
    real(r8kind)      :: addisl = 0.5               ! parameter used in the computation of the dissipation timescale at the lower level over land
    real(r8kind)      :: addish = 0.5               ! parameter used in the computation of the dissipation timescale at the lower level as a function of topography
    real(r8kind)      :: trel = 25.0                ! relaxation time scale in days of the temperature
    real(r8kind)      :: tdif = 3.0                 ! dissipation timescale of scale-selective diffusion in days for wavenumber nm
    integer           :: idif = 4                   ! determines scale-selectivity of hyperviscosity; power of laplace operator
    real(r8kind)      :: h0 = 3.0                   ! scale factor for the topographically induced upward motion at the lower level
    real(r8kind)      :: rrdef1 = 0.110             ! Rossby radius of deformation of 200-500 thickness
    real(r8kind)      :: rrdef2 = 0.070             ! Rossby radius of deformation of 500-800 thickness

    namelist /params/ resolution, time_step, obsfile, inf, obsf, tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2

    ! Open the namelist file
    open(newunit=nl_unit, file=trim(nl_filename), form='formatted', status='old')

    ! Read the configuration
    read(nl_unit, nml=params)

    ! Close the namelist
    close(nl_unit)

    qg_config = qg_config_type(resolution, time_step, obsfile, inf, obsf, tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2)

  end function constructor_namelist_file


  !------------------------------------------------------------------
  ! constructor_namelist_unit
  !
  ! Returns an initialized qg_config_type object
  !------------------------------------------------------------------
  function constructor_namelist_unit(nl_unit) result(qg_config)

    ! Namelist unit number (must be an open file)
    integer, intent(in)  :: nl_unit
    type(qg_config_type) :: qg_config

    ! Namelist param variables
    integer           :: resolution = 21            ! Model resolution
    integer           :: time_step = 1200           ! Model time step (in seconds)
    character(len=32) :: obsfile ='sf7910T106.shfs' ! Name of observation file
    logical           :: inf = .false.              ! If .true. then artificial PV forcing read from file
    logical           :: obsf = .false.             ! If .true. PV forcing is calculated from observations in routine artiforc
    real(r8kind)      :: tdis = 3.0                 ! Ekman dissipation timescale in days at lower level
    real(r8kind)      :: addisl = 0.5               ! parameter used in the computation of the dissipation timescale at the lower level over land
    real(r8kind)      :: addish = 0.5               ! parameter used in the computation of the dissipation timescale at the lower level as a function of topography
    real(r8kind)      :: trel = 25.0                ! relaxation time scale in days of the temperature
    real(r8kind)      :: tdif = 3.0                 ! dissipation timescale of scale-selective diffusion in days for wavenumber nm
    integer           :: idif = 4                   ! determines scale-selectivity of hyperviscosity; power of laplace operator
    real(r8kind)      :: h0 = 3.0                   ! scale factor for the topographically induced upward motion at the lower level
    real(r8kind)      :: rrdef1 = 0.110             ! Rossby radius of deformation of 200-500 thickness
    real(r8kind)      :: rrdef2 = 0.070             ! Rossby radius of deformation of 500-800 thickness

    namelist /params/ resolution, time_step, obsfile, inf, obsf, tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2

    ! Read the configuration
    read(nl_unit, nml=params)

    ! Rewind the namelist
    rewind(nl_unit)

    ! Initialize model parameters    
    qg_config = qg_config_type(resolution, time_step, obsfile, inf, obsf, tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2)

  end function constructor_namelist_unit


  !------------------------------------------------------------------
  ! destructor_qg_config
  !
  ! Deallocates pointers used by a qg_config_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_qg_config(this)

    type(qg_config_type), intent(inout) :: this

    ! No pointers in qg_config_type object so we do nothing

  end subroutine destructor_qg_config


  !-------------------------------------------------------------------------------
  ! init_grid
  !-------------------------------------------------------------------------------
  subroutine init_grid(this, resolution)

    class(qg_config_type) :: this
    integer, intent(in)   :: resolution

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

  end subroutine init_grid


  !-------------------------------------------------------------------------------
  ! init_spectral_coeff
  !-------------------------------------------------------------------------------
  subroutine init_spectral_coeff(this)

    class(qg_config_type), intent(inout) :: this

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

    class(qg_config_type), intent(inout) :: this
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

    class(qg_config_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: h0
    real(r8kind),         intent(   in) :: addisl
    real(r8kind),         intent(   in) :: addish

    real(r8kind), allocatable :: fmu(:,:)
    real(r8kind), allocatable :: agg(:,:), agg1(:,:), agg2(:,:)
    real(r8kind), allocatable :: orog(:)
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
    allocate(orog(this%nsh2))
    allocate(ws(this%nsh2))
    allocate(this%dorodl(this%nlat,this%nlon))
    allocate(this%dorodm(this%nlat,this%nlon))
    allocate(agg2(this%nlat, this%nlon))
    this%lgdiss = ((addisl .gt. 0.0) .or. (addish .gt. 0.0))
    orog = reshape(this%ggsp%ggtosp (agg), (/this%nsh2/))
    ws = reshape(this%ggsp%ddl (orog), (/this%nsh2/))
    this%dorodl = this%ggsp%sptogg_pp (ws)
    this%dorodm = this%ggsp%sptogg_pd (orog)
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

    class(qg_config_type),   intent(inout) :: this
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


  !-----------------------------------------------------------------------
  ! computation of artifical forcing according to roads(1987)
  ! the forcing is computed from file obsfile
  !------------------------------------------------------------------------
  function artiforc(this) result(for)

    class(qg_config_type), intent(in) :: this

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

    obsfile = this%obsfile

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
! Put this back in   psi = this%fstofm(psifs, this%nm)
! Put this back in   call this%psitoq(psi, psit, qprime)
! Put this back in   dqprdt = this%ddt(psi, psit, qprime, for)
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
  ! get_resolution
  !
  ! Retun config resolution
  !------------------------------------------------------------------
  function get_resolution(this) result(resolution)

    class(qg_config_type), intent(in) :: this
    integer                           :: resolution

    ! Return resolution
    resolution = this%resolution

  end function get_resolution


  !------------------------------------------------------------------
  ! get_time_step
  !
  ! Retun config time_step
  !------------------------------------------------------------------
  function get_time_step(this) result(time_step)

    class(qg_config_type), intent(in) :: this
    integer                           :: time_step

    ! Return time_step
    time_step = this%time_step

  end function get_time_step


  !------------------------------------------------------------------
  ! get_obsfile
  !
  ! Retun config obsfile
  !------------------------------------------------------------------
  function get_obsfile(this) result(obsfile)

    class(qg_config_type), intent(in) :: this
    character(len=32)                 :: obsfile

    ! Return inf
    obsfile = this%obsfile

  end function get_obsfile


  !------------------------------------------------------------------
  ! get_inf
  !
  ! Retun config inf
  !------------------------------------------------------------------
  function get_inf(this) result(inf)

    class(qg_config_type), intent(in) :: this
    logical                           :: inf

    ! Return inf
    inf = this%inf

  end function get_inf


  !------------------------------------------------------------------
  ! get_obsf
  !
  ! Retun config obsf
  !------------------------------------------------------------------
  function get_obsf(this) result(obsf)

    class(qg_config_type), intent(in) :: this
    logical                           :: obsf

    ! Return obsf
    obsf = this%obsf

  end function get_obsf


  !------------------------------------------------------------------
  ! get_tdis
  !
  ! Retun config tdis
  !------------------------------------------------------------------
  function get_tdis(this) result(tdis)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: tdis

    ! Return tdis
    tdis = this%tdis

  end function get_tdis


  !------------------------------------------------------------------
  ! get_addisl
  !
  ! Retun config addisl
  !------------------------------------------------------------------
  function get_addisl(this) result(addisl)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: addisl

    ! Return addisl
    addisl = this%addisl

  end function get_addisl


  !------------------------------------------------------------------
  ! get_addish
  !
  ! Retun config addish
  !------------------------------------------------------------------
  function get_addish(this) result(addish)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: addish

    ! Return addish
    addish = this%addish

  end function get_addish


  !------------------------------------------------------------------
  ! get_trel
  !
  ! Retun config trel
  !------------------------------------------------------------------
  function get_trel(this) result(trel)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: trel

    ! Return trel
    trel = this%trel

  end function get_trel


  !------------------------------------------------------------------
  ! get_tdif
  !
  ! Retun config tdif
  !------------------------------------------------------------------
  function get_tdif(this) result(tdif)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: tdif

    ! Return tdif
    tdif = this%tdif

  end function get_tdif


  !------------------------------------------------------------------
  ! get_idif
  !
  ! Retun config idif
  !------------------------------------------------------------------
  function get_idif(this) result(idif)

    class(qg_config_type), intent(in) :: this
    integer                           :: idif

    ! Return idif
    idif = this%idif

  end function get_idif


  !------------------------------------------------------------------
  ! get_h0
  !
  ! Retun config h0
  !------------------------------------------------------------------
  function get_h0(this) result(h0)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: h0

    ! Return h0
    h0 = this%h0

  end function get_h0


  !------------------------------------------------------------------
  ! get_rrdef1
  !
  ! Retun config rrdef1
  !------------------------------------------------------------------
  function get_rrdef1(this) result(rrdef1)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: rrdef1

    ! Return rrdef1
    rrdef1 = this%rrdef1

  end function get_rrdef1


  !------------------------------------------------------------------
  ! get_rrdef2
  !
  ! Retun config rrdef2
  !------------------------------------------------------------------
  function get_rrdef2(this) result(rrdef2)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: rrdef2

    ! Return rrdef2
    rrdef2 = this%rrdef2

  end function get_rrdef2


  !------------------------------------------------------------------
  ! get_ggsp
  !
  ! Retun config ggsp
  !------------------------------------------------------------------
  function get_ggsp(this) result(ggsp)

    class(qg_config_type), intent(in) :: this
    type(qg_ggsp_type)                :: ggsp

    ! Return ggsp
    ggsp = this%ggsp

  end function get_ggsp


  !------------------------------------------------------------------
  ! get_for
  !
  ! Retun config for
  !------------------------------------------------------------------
  function get_for(this) result(for)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: for(this%nsh2, this%nvl)

    ! Return for
    for = this%for

  end function get_for


  !------------------------------------------------------------------
  ! get_nm
  !
  ! Retun config nm
  !------------------------------------------------------------------
  function get_nm(this) result(nm)

    class(qg_config_type), intent(in) :: this
    integer                           :: nm

    ! Return nm
    nm = this%nm

  end function get_nm

  !------------------------------------------------------------------
  ! get_nsh2
  !
  ! Retun config nsh2
  !------------------------------------------------------------------
  function get_nsh2(this) result(nsh2)

    class(qg_config_type), intent(in) :: this
    integer                           :: nsh2

    ! Return nsh2
    nsh2 = this%nsh2

  end function get_nsh2


  !------------------------------------------------------------------
  ! get_nsh
  !
  ! Retun config nsh
  !------------------------------------------------------------------
  function get_nsh(this) result(nsh)

    class(qg_config_type), intent(in) :: this
    integer                           :: nsh

    ! Return nsh
    nsh = this%nsh

  end function get_nsh

  !------------------------------------------------------------------
  ! get_nvl
  !
  ! Retun config nvl
  !------------------------------------------------------------------
  function get_nvl(this) result(nvl)

    class(qg_config_type), intent(in) :: this
    integer                           :: nvl

    ! Return nvl
    nvl = this%nvl

  end function get_nvl

  !------------------------------------------------------------------
  ! get_ntl
  !
  ! Retun config ntl
  !------------------------------------------------------------------
  function get_ntl(this) result(ntl)

    class(qg_config_type), intent(in) :: this
    integer                           :: ntl

    ! Return ntl
    ntl = this%ntl

  end function get_ntl


  !------------------------------------------------------------------
  ! get_nlat
  !
  ! Retun config nlat
  !------------------------------------------------------------------
  function get_nlat(this) result(nlat)

    class(qg_config_type), intent(in) :: this
    integer                           :: nlat

    ! Return nlat
    nlat = this%nlat

  end function get_nlat


  !------------------------------------------------------------------
  ! get_nlon
  !
  ! Retun config nlon
  !------------------------------------------------------------------
  function get_nlon(this) result(nlon)

    class(qg_config_type), intent(in) :: this
    integer                           :: nlon

    ! Return nlon
    nlon = this%nlon

  end function get_nlon


  !------------------------------------------------------------------
  ! get_ft
  !
  ! Retun config ft
  !------------------------------------------------------------------
  function get_ft(this) result(ft)

    class(qg_config_type), intent(in) :: this
    character(len=3)                  :: ft

    ! Return ft
    ft = this%ft

  end function get_ft


  !------------------------------------------------------------------
  ! get_ll
  !
  ! Retun config ll
  !------------------------------------------------------------------
  function get_ll(this) result(ll)

    class(qg_config_type), intent(in) :: this
    integer                           :: ll(this%nsh)

    ! Return ll
    ll = this%ll

  end function get_ll


  !------------------------------------------------------------------
  ! get_relt1
  !
  ! Retun config relt1
  !------------------------------------------------------------------
  function get_relt1(this) result(relt1)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: relt1

    ! Return relt1
    relt1 = this%relt1

  end function get_relt1


  !------------------------------------------------------------------
  ! get_relt2
  !
  ! Retun config relt2
  !------------------------------------------------------------------
  function get_relt2(this) result(relt2)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: relt2

    ! Return relt2
    relt2 = this%relt2

  end function get_relt2


  !------------------------------------------------------------------
  ! get_diss
  !
  ! Retun config diss
  !------------------------------------------------------------------
  function get_diss(this) result(diss)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: diss(this%nsh2, 2)

    ! Return diss
    diss = this%diss

  end function get_diss


  !------------------------------------------------------------------
  ! get_phi
  !
  ! Retun config phi
  !------------------------------------------------------------------
  function get_phi(this) result(phi)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: phi(this%nlat)

    ! Return phi
    phi = this%phi

  end function get_phi


  !------------------------------------------------------------------
  ! get_sinfi
  !
  ! Retun config sinfi
  !------------------------------------------------------------------
  function get_sinfi(this) result(sinfi)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: sinfi(this%nlat)

    ! Return sinfi
    sinfi = this%sinfi

  end function get_sinfi


  !------------------------------------------------------------------
  ! get_cosfi
  !
  ! Retun config cosfi
  !------------------------------------------------------------------
  function get_cosfi(this) result(cosfi)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: cosfi(this%nlat)

    ! Return cosfi
    cosfi = this%cosfi

  end function get_cosfi


  !------------------------------------------------------------------
  ! get_dorodl
  !
  ! Retun config dorodl
  !------------------------------------------------------------------
  function get_dorodl(this) result(dorodl)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: dorodl(this%nlat, this%nlon)

    ! Return dorodl
    dorodl = this%dorodl

  end function get_dorodl


  !------------------------------------------------------------------
  ! get_dorodm
  !
  ! Retun config dorodm
  !------------------------------------------------------------------
  function get_dorodm(this) result(dorodm)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: dorodm(this%nlat, this%nlon)

    ! Return dorodm
    dorodm = this%dorodm

  end function get_dorodm


  !------------------------------------------------------------------
  ! get_lgdiss
  !
  ! Retun config lgdiss
  !------------------------------------------------------------------
  function get_lgdiss(this) result(lgdiss)

    class(qg_config_type), intent(in) :: this
    logical                           :: lgdiss

    ! Return lgdiss
    lgdiss = this%lgdiss

  end function get_lgdiss


  !------------------------------------------------------------------
  ! get_ddisdx
  !
  ! Retun config ddisdx
  !------------------------------------------------------------------
  function get_ddisdx(this) result(ddisdx)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: ddisdx(this%nlat, this%nlon)

    ! Return ddisdx
    ddisdx = this%ddisdx

  end function get_ddisdx


  !------------------------------------------------------------------
  ! get_ddisdy
  !
  ! Retun config ddisdy
  !------------------------------------------------------------------
  function get_ddisdy(this) result(ddisdy)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: ddisdy(this%nlat, this%nlon)

    ! Return ddisdy
    ddisdy = this%ddisdy

  end function get_ddisdy


  !------------------------------------------------------------------
  ! get_rdiss
  !
  ! Retun config rdiss
  !------------------------------------------------------------------
  function get_rdiss(this) result(rdiss)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: rdiss(this%nlat, this%nlon)

    ! Return rdiss
    rdiss = this%rdiss

  end function get_rdiss


  !------------------------------------------------------------------
  ! get_rinhel
  !
  ! Retun config rinhel
  !------------------------------------------------------------------
  function get_rinhel(this) result(rinhel)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: rinhel(this%nsh2,0:5)

    ! Return rinhel
    rinhel = this%rinhel

  end function get_rinhel


  !------------------------------------------------------------------
  ! get_rl1
  !
  ! Retun config rl1
  !------------------------------------------------------------------
  function get_rl1(this) result(rl1)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: rl1

    ! Return rl1
    rl1 = this%rl1

  end function get_rl1


  !------------------------------------------------------------------
  ! get_rl2
  !
  ! Retun config rl2
  !------------------------------------------------------------------
  function get_rl2(this) result(rl2)

    class(qg_config_type), intent(in) :: this
    real(r8kind)                      :: rl2

    ! Return rl2
    rl2 = this%rl2

  end function get_rl2


end module QG_Config
