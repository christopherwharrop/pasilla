module QG

  use kind

  implicit none

  ! Model configuration parameters
  character(len=3) ::  ft   ! Character string containing resolution
  integer :: nm             ! The truncation is of type T(riangular) nm
  integer :: nlon           ! Number of longitude points of the Gaussian grid
  integer :: nlat           ! Number of latitude  points of the Gaussian grid
  integer :: nvl            ! Number of vorticity levels in the vertical (should be set to 3)
  integer :: ntl            ! Number of temperature levels in the vertical (equal to nvl-1)
  integer :: nsh            ! Half of nsh2
  integer :: nsh2           ! Number of coefficients needed to define one level of the T nm model
  integer :: ngp            ! Number of grid points of the Gaussian grid
  real(r8kind)  :: dt       ! timestep in fraction of one day
  real(r8kind)  :: dtime    ! timestep in seconds
  real(r8kind)  :: dtt      ! timestep in non-dimensional units


  ! Model runtime parameters
  character(len=32)  :: obsfile  ! Name of observation file
  character(len=256) :: rootdir  ! Path of qgmodel directory
  integer :: nstepsperday
  integer :: nstepsbetweenoutput
  integer :: ndayskip
  integer :: nday

  ! Mathematical and physical constants
  real(r8kind) :: pi     ! value of pi
  real(r8kind) :: fzero  ! value of f at 45 degrees
  real(r8kind) :: dp     ! layer thicknesses [Pa]
  real(r8kind) :: om     ! angular velocity of earth [rad/s]
  real(r8kind) :: rgas   ! gas constant
  real(r8kind) :: grav   ! gravity acceleration [m/s^2]
  real(r8kind) :: radius ! radius of earth [m]

  real(r8kind), allocatable  :: phi(:)   ! Gauss points in radians
  real(r8kind), allocatable  :: cosfi(:) ! cosine of phi
  real(r8kind), allocatable  :: sinfi(:) ! sine of phi

  integer, allocatable :: nshm(:) ! contains numbers 22 down to 1 for index 0 to 21
  integer, allocatable :: ll(:)   ! contains total wavenumber n of each spherical harmonic of the corresponding index

  real(r8kind), allocatable  :: rm(:)       ! contains zonal wavenumber m of each spherical harmonic of the corresponding index for zonal derivative operator
  real(r8kind), allocatable  :: rinhel(:,:) ! Laplace and Helmholtz operator for Q-PSI inversion
  real(r8kind), allocatable  :: diss(:,:)   ! dissipation coefficients for each spherical harmonic
                                            !   diss(k,1) : hyperviscosity at the three levels
                                            !               (tdif sets timescale)
                                            !   diss(k,2) : Ekman friction at lower level
                                            !               (tdis sets timescale)

  logical :: lgdiss     ! if .true. then orography and land-sea mask dependent friction at the lower level plus Ekman friction, else only Ekman friction
  logical :: inf        ! if .true. then artificial PV forcing read from file
  logical :: obsf       ! if .true. PV forcing is calculated from observations in routine artiforc
  logical :: readstart  ! if .true. initial state is read from inputfile

  real(r8kind), allocatable  :: pp(:,:)  ! Legendre polynomials defined at Gausian latitudes
  real(r8kind), allocatable  :: pd(:,:)  ! mu derivative of Legendre polynomials
  real(r8kind), allocatable  :: pw(:,:)  ! weights for Legendre integrals

  real(r8kind), allocatable  :: rdiss(:,:)   ! landsea-mask/orography dependent friction
  real(r8kind), allocatable  :: ddisdx(:,:)  ! landsea-mask/orography dependent friction
  real(r8kind), allocatable  :: ddisdy(:,:)  ! landsea-mask/orography dependent friction

  real(r8kind)  :: rrdef1  ! Rossby radius of deformation of 200-500 thickness
  real(r8kind)  :: rrdef2  ! Rossby radius of deformation of 500-800 thickness
  real(r8kind)  :: rl1     ! one over Rossby rad. of def. squared of 200-500 thickness
  real(r8kind)  :: rl2     ! one over Rossby rad. of def. squared of 500-800 thickness
  real(r8kind)  :: relt1   ! nondimensional relaxation coefficient of 200-500 thickness
  real(r8kind)  :: relt2   ! nondimensional relaxation coefficient of 500-800 thickness
  real(r8kind)  :: tdis    ! Ekman dissipation timescale in days at lower level
  real(r8kind)  :: trel    ! relaxation time scale in days of the temperature
  real(r8kind)  :: tdif    ! dissipation timescale of scale-selective diffusion in days for wavenumber nm
  real(r8kind)  :: addisl  ! parameter used in the computation of the dissipation timescale at the lower level over land
  real(r8kind)  :: addish  ! parameter used in the computation of the dissipation timescale at the lower level as a function of topography
  real(r8kind)  :: h0      ! scale factor for the topographically induced upward motion at the lower level
  integer       :: idif    ! determines scale-selectivity of hyperviscosity; power of laplace operator

  real(r8kind), allocatable  :: psi(:,:)    ! stream function at the nvl levels
  real(r8kind), allocatable  :: psit(:,:)   ! thickness at the ntl levels
  real(r8kind), allocatable  :: qprime(:,:) ! potential vorticity
  real(r8kind), allocatable  :: dqprdt(:,:) ! time derivative of qprime
  real(r8kind), allocatable  :: for(:,:)    ! constant potential vorticity forcing at the nvl levels
  real(r8kind), allocatable  :: ws(:)       ! only used as portable workspace

  real(r8kind), allocatable  :: orog(:)     ! orography in m. divided by h0
  real(r8kind), allocatable  :: dorodl(:,:) ! derivative of orog wrt lambda
  real(r8kind), allocatable  :: dorodm(:,:) ! derivative of orag wrt sin(fi)

  real(r8kind), allocatable  :: trigd(:,:)  ! array used by the nag version of the fft
  real(r8kind), allocatable  :: trigi(:,:)  ! array used by the nag version of the fft
  real(r8kind), allocatable  :: wgg(:,:)    ! array used by the nag version of the fft

  real(r8kind), allocatable  :: psig(:,:,:)   ! grid values of dimensional streamfunction at the three levels
  real(r8kind), allocatable  :: qgpv(:,:,:)   ! grid values of dimensional pv at the three levels
  real(r8kind), allocatable  :: ug(:,:,:)     ! grid values of zonal velocity at the three levels in m/s
  real(r8kind), allocatable  :: vg(:,:,:)     ! grid values of meridional velocity at the three levels in m/s
  real(r8kind), allocatable  :: geopg(:,:,:)  ! geopotential on the grid

contains

  !-------------------------------------------------------------------------------
  ! allocate_comqg
  !-------------------------------------------------------------------------------
  subroutine allocate_comqg(resolution)

    integer, intent(in) :: resolution

    select case (resolution)

      case(106)
        nm = 106
        nlon = 320
        nlat = 160
        ft = "106"

      case(63)
        nm = 63
        nlon = 192
        nlat = 96
        ft = "63"

      case(42)
        nm = 42
        nlon = 128
        nlat = 64
        ft = "42"

      case(21)
        nm = 21
        nlon = 64
        nlat = 32
        ft = "21"

      case DEFAULT
        stop 'ERROR: Unsupported resolution'

    end select

    nvl = 3
    ntl = nvl - 1
    nsh = ((nm + 1) * (nm + 2)) / 2
    nsh2 = 2 * nsh
    ngp = nlon * nlat

    allocate(phi(nlat), cosfi(nlat), sinfi(nlat))
    allocate(nshm(0:nm), ll(nsh))
    allocate(rm(nsh), rinhel(nsh2,0:5), diss(nsh2,2))
    allocate(pp(nlat,nsh), pd(nlat,nsh), pw(nlat,nsh))
    allocate(rdiss(nlat,nlon), ddisdx(nlat,nlon), ddisdy(nlat,nlon))

    allocate(psi(nsh2,nvl), psit(nsh2,ntl))
    allocate(qprime(nsh2,nvl), dqprdt(nsh2,nvl), for(nsh2,nvl), ws(nsh2))
    allocate(orog(nsh2), dorodl(nlat,nlon), dorodm(nlat,nlon))
    allocate(trigd(nlon,2), trigi(nlon,2), wgg(nlat,nlon))
    allocate(psig(nlat,nlon,nvl), qgpv(nlat,nlon,nvl))
    allocate(ug(nlat,nlon,nvl), vg(nlat,nlon,nvl), geopg(nlat,nlon,nvl))

  end subroutine allocate_comqg


  !-------------------------------------------------------------------------------
  ! initqg
  !
  ! initialise parameters and operators and read initial state
  !-------------------------------------------------------------------------------
  subroutine initqg

    implicit none

    integer      :: resolution
    integer      :: i, j, k1, k2, k, l, m, n, ifail, ii, jj, i1, j1, nn
    real(r8kind) :: pigr4, dis, dif, rll
    real(r8kind) :: r1, a, b, c, d, e, sqn, rsqn
    real(r8kind) :: rnorm, rh0, dd, dlon
    real(r8kind), allocatable :: ininag(:,:)
    real(r8kind), allocatable :: agg(:,:), agg1(:,:), agg2(:,:) 
    real(r8kind), allocatable :: fmu(:,:), wsx(:)

    namelist /param/ tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2
    namelist /control/ resolution, nstepsperday, nstepsbetweenoutput, &
   &                   ndayskip, nday, obsfile, inf, obsf, readstart


    ! Default control namelist values
    resolution = 21
    nstepsperday = 36
    nstepsbetweenoutput = 36
    ndayskip = 0
    nday = 10
    obsfile = 'sf7910T106.shfs'
    inf = .false.
    obsf = .false.
    readstart = .false.

    ! Default param namelist values
    tdis = 3.0
    addisl = 0.5
    addish = 0.5
    trel = 25.0
    tdif = 3.0
    idif = 4
    h0 = 3.0
    rrdef1 = 0.110
    rrdef2 = 0.070

    ! Read innput namelist parameters
    open(15, file = 'namelist.input', status = 'old', form = 'formatted')
    read(15, nml = control)
    read(15, nml = param)
    close(15)

    ! Write namelist parameters
    open(16, file = 'namelist.output')
    write(16, nml = control)
    write(16, nml = param)
    close(16)

    ! Set model resolution dependent parameters and allocate model state variables
    call allocate_comqg(resolution)

    ! Allocate local data
    allocate(ininag(nlat, nlon))
    allocate(agg(nlat, nlon), agg1(nlat, nlon), agg2(nlat, nlon))
    allocate(fmu(nlat, 2), wsx(nsh2))

    ! Read model input from qgcoefT*
    open(11, file = './qgcoefT' // trim(ft) // '.dat', form = 'formatted')
    do i = 0, nm
      read(11, *) nshm(i)
    enddo
    do i = 1, nsh
      read(11, *) ll(i)
    enddo
    do k = 1, nsh
      do j = 1, nlat
        read(11, *) pp(j, k)
      enddo
    enddo
    do k = 1, nsh
      do j = 1, nlat
        read(11, *) pd(j, k)
      enddo
    enddo
    do k = 1, nsh
      do j = 1, nlat
        read(11, *) pw(j, k)
      enddo
    enddo
    close(11)

    ! Read initial streamfunction
    open(12, file = './qgstartT' // trim(ft) // '.dat', form = 'formatted')
    open(13, file = './qgbergT' // trim(ft) // '.dat', form = 'formatted')
    if (inf) then
      open(14, file = './qgpvforT' // trim(ft) // '.dat', form = 'formatted')
    endif


    pi = 4d0 * atan(1d0)
    radius = 6.37e+6 
    om = 4d0 * pi / (24d0 * 3600d0)

    pigr4 = 4.d0 * pi
    rl1 = 1.0d0 / rrdef1**2
    rl2 = 1.0d0 / rrdef2**2
    relt1 = max(0.0d0, rl1 / (trel * pigr4))
    relt2 = max(0.0d0, rl2 / (trel * pigr4))
    dis = max(0.0d0, 1.0d0 / (tdis * pigr4))
    rll = dble(ll(nsh))
    dif = max(0.0d0, 1.0d0 / (tdif * pigr4 * (rll * (rll + 1))**idif))

    ! time step of the model: 
    ! dt    : fraction of one day
    ! dtime : in seconds
    ! dtt   : dimensionless
    dt     = 1d0 / real(nstepsperday)
    dtime  = dt * (24d0 * 3600d0)
    dtt    = dt * pi * 4d0

    ! zonal derivative operator
    k2 = 0
    do m = 0, nm
      k1 = k2 + 1
      k2 = k2 + nshm(m)
      do k = k1, k2
        rm(k) = dble(m)
      enddo
    enddo

    ! laplace/helmholtz direct and inverse operators
    do j = 0, 5
      rinhel(1, j) = 0.0d0
    enddo
    diss(1, 1) = 0.0d0
    diss(1, 2) = 0.0d0
    do k = 2, nsh
      r1 = dble(ll(k) * (ll(k) + 1))
      a = -r1 - 3.0d0 * rl1
      b = -r1 - 3.0d0 * rl2
      c = -r1 - rl1
      d = -r1 - rl2
      e = a * d + b * c
      rinhel(k, 0) = -r1
      rinhel(k, 1) = -1.0d0 / r1
      rinhel(k, 2) =  d / e
      rinhel(k, 3) =  b / e
      rinhel(k, 4) = -c / e
      rinhel(k, 5) =  a / e
      diss(k, 2) = dis * r1
      diss(k, 1) = -dif * r1**idif
    enddo
    do j = 0, 5
      do k = 1, nsh
        rinhel(k + nsh, j) = rinhel(k, j)
      enddo
    enddo
    do j = 1, 2
      do k = 1, nsh
        diss(k + nsh, j) = diss(k, j)
      enddo
    enddo


    ! compensation for normalization in nag fft routines
    sqn = sqrt(dble(nlon))
    rsqn = 1d0 / sqn
    do k = 1, nsh
      do i = 1, nlat
        pp(i, k) = pp(i, k) * sqn
        pd(i, k) = pd(i, k) * sqn
        pw(i, k) = pw(i, k) * rsqn
      enddo
    enddo

    ! initialization of coefficients for fft
    do j = 1, nlon
      do i = 1, nlat
        ininag(i, j) = 1.0d0
      enddo
    enddo 
    ifail = 0
    call c06fpf(nlat, nlon, ininag, 'i', trigd, wgg, ifail)
    ifail = 0
    call c06fqf(nlat, nlon, ininag, 'i', trigi, wgg, ifail)

    ! orography and dissipation terms
    ! fmu(i, 1): sin(phi(i))
    ! fmu(i, 2): 1 - sin**2(phi(i))      
    rnorm = 1.0d0 / sqrt(3.0d0 * nlon)
    do i = 1, nlat
      fmu(i, 1) = rnorm * pp(i, 2)
      fmu(i, 2) = 1.d0 - fmu(i, 1)**2
      sinfi(i) = fmu(i, 1)
      phi(i) = asin(sinfi(i))
      cosfi(i) = cos(phi(i))
      phi(i) = 180d0 * phi(i) / pi
    enddo
    dlon = 360d0 / real(nlon)

    ! height of orography in meters
    do i = 1, nlon
      do j = 1, nlat
        read(13, *) agg1(J, I)
      enddo
    enddo
    rh0 = max(0.0d0, 0.001d0 / h0)
    do j = 1, nlon
      do i = 1, nlat
        agg(i, j) = fmu(i, 1) * agg1(i, j) * rh0
    !      agg(i, j) = agg1(i, j) * rh0
      enddo
    enddo

    ! surface dependent friction
    lgdiss = ((addisl .gt. 0.0) .or. (addish .gt. 0.0))
    orog = reshape(ggtosp (agg), (/nsh2/))
    call ddl (orog, ws)
    dorodl = sptogg (ws, pp)
    dorodm = sptogg (orog, pd)
    if (lgdiss) then
      do i = 1, nlon
        do j = 1, nlat
          read(13, *) agg2(j, i)
        enddo
      enddo
      do j = 1, nlon
        do i = 1, nlat
          agg(i, j) = 1.0d0 + addisl * agg2(i, j) + addish * (1.0d0 - exp(-0.001d0 * agg1(i, j)))
        enddo
      enddo

      ws = reshape(ggtosp (agg), (/nsh2/))
      call ddl (ws, wsx)
      rdiss = sptogg (ws, pp)
      ddisdx = sptogg (wsx, pp)
      ddisdy = sptogg (ws, pd)
      dd = 0.5d0 * diss(2, 2)
      do j = 1, nlon
        do i = 1, nlat
          ddisdx(i, j) = dd * ddisdx(i, j) / fmu(i, 2)
          ddisdy(i, j) = dd * ddisdy(i, j) * fmu(i, 2)
        enddo
      enddo

    endif

    ! forcing term
    do l = 1, 3
      do k = 1, nsh2
        for(k, l) = 0d0
      enddo
    enddo
    if (inf) then
      read(14, '(1e12.5)') ((for(k, l), k = 1, nsh2), l = 1, 3)
    endif
    if (obsf) then
      call artiforc
    endif

    if (readstart) then
      do l = 1, 3
        do k = 1, nsh2
          read(12, *) psi(k, l)
        enddo
      enddo
    else
      do l = 1, 3
        do k = 1, nsh2
          psi(k, l) = 0d0
        enddo
      enddo
    endif
    close(12)

    ! Potential vorticity and streamfunction fields
    call psitoq

    close(13)
    close(14)

    open(13, file = 'qgbergT' // trim(ft) // '.grads', form = 'unformatted')
    write(13) ((real(agg1(j, i)), i = 1, nlon), j = 1, nlat)
    write(13) ((real(agg2(j, i)), i = 1, nlon), j = 1, nlat)
    close(13)
    open(50, file = 'qgbergT' // trim(ft) // '.ctl', form = 'formatted')
    write(50, '(A)') 'dset ^qgbergT' // trim(ft) // '.grads'
    write(50, '(A)') 'undef 9.99e+10'
    write(50, '(A)') 'options sequential big_endian'
    write(50, '(A)') 'title three level QG model'
    write(50, '(A)') '*'
    write(50, '(A, i4, A, F19.14)') 'xdef ', nlon, ' linear  0.000 ', dlon
    write(50, '(A)') '*'
    write(50, '(A, I4, A, 1F19.14)') 'ydef ', nlat, ' levels ', phi(1)
    write(50, '(F19.14)') (phi(j), j = 2, nlat)
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

    open(14, file = 'qgpvforT' // trim(ft) // '.grads', form = 'unformatted')
    do l = 1, nvl
      agg1 = sptogg(for(1, l), pp)
      write(14) ((real(agg1(j, i)), i = 1, nlon), j = 1, nlat)
    enddo
    close(14)

    open(50, file = 'qgpvforT' // trim(ft) // '.ctl', form = 'formatted')
    write(50, '(A)') 'dset ^qgpvforT' // trim(ft) // '.grads'
    write(50, '(A)') 'undef 9.99e+10'
    write(50, '(A)') 'options sequential big_endian'
    write(50, '(A)') 'title three level QG model'
    write(50, '(A)') '*'
    write(50, '(A, i4, A, F19.14)') 'xdef ', nlon, ' linear  0.000 ', dlon
    write(50, '(A)') '*'
    write(50, '(A, I4, A, 1F19.14)') 'ydef ', nlat, ' levels ', phi(1)
    write(50, '(F19.14)') (phi(j), j = 2, nlat)
    write(50, '(A)') '*'
    write(50, '(A)') 'zdef  3 levels 800 500 200'
    write(50, '(A)') '*'
    write(50, '(A)') 'tdef 1 linear 1jan0001 1dy'
    write(50, '(A)') '*'
    write(50, '(A)') 'vars  1'
    write(50, '(A)') 'pvfor    3  99 pv forcing field [nondim]'
    write(50, '(A)') 'endvars'

    close(50)

    return

  end subroutine initqg


  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  ! NOTE psit is destroyed
  !----------------------------------------------------------------------
  subroutine ddt

    implicit none

    integer :: k, l, i, j
    real(r8kind) :: dum1, dum2

    ! advection of potential vorticity at upper level
    call jacob (psi(1, 1), qprime(1, 1), dqprdt(1, 1))

    ! advection of potential vorticity at middle level
    call jacob (psi(1, 2), qprime(1, 2), dqprdt(1, 2))

    ! advection of potential vorticity and dissipation at lower level
    call jacobd (psi(1, 3), qprime(1, 3), dqprdt(1, 3))

    ! relaxation of temperature and forcing
    do k = 1, nsh2
      dum1 = relt1 * psit(k, 1)
      dum2 = relt2 * psit(k, 2)
      dqprdt(k, 1) = dqprdt(k, 1) + dum1        + for(k, 1)
      dqprdt(k, 2) = dqprdt(k, 2) - dum1 + dum2 + for(k, 2)
      dqprdt(k, 3) = dqprdt(k, 3)        - dum2 + for(k, 3)
    enddo

    ! explicit horizontal diffusion
    do l = 1, 3
      do k = 1, nsh2
        dqprdt(k, l) = dqprdt(k, l) + diss(k, 1) * qprime(k, l)
      enddo
    enddo

    return

  end subroutine ddt


  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  subroutine jacob (psiloc, pvor, sjacob)

    implicit none

    real(r8kind), intent( in) :: psiloc(nsh2)
    real(r8kind), intent( in) :: pvor(nsh2)
    real(r8kind), intent(out) :: sjacob(nsh2)

    integer      :: i, j, k
    real(r8kind) :: vv(nsh2)
    real(r8kind) :: dpsidl(nlat, nlon),  dpsidm(nlat, nlon),  dvordl(nlat, nlon)
    real(r8kind) :: dvordm(nlat, nlon),  gjacob(nlat, nlon),  dpsidls(nsh2)

    ! space derivatives of potential vorticity
    call ddl (pvor, vv)
    dvordl = sptogg (vv, pp)
    dvordm = sptogg (pvor, pd)

    ! space derivatives of streamfunction
    call ddl (psiloc, dpsidls)
    dpsidl = sptogg (dpsidls, pp)
    dpsidm = sptogg (psiloc, pd)

    ! jacobian term
    do j = 1, nlon
      do i = 1, nlat
        gjacob(i, j) = dpsidm(i, j) * dvordl(i, j) - dpsidl(i, j) * dvordm(i, j)
      enddo
    enddo

    sjacob = reshape(ggtosp (gjacob), (/nsh2/))

    ! planetary vorticity advection
    do k = 1, nsh2
      sjacob(k) = sjacob(k) - dpsidls(k)
    enddo

    return

  end subroutine jacob


  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  subroutine jacobd (psiloc, pvor, sjacob)

    implicit none

    real(r8kind), intent( in) :: psiloc(nsh2)
    real(r8kind), intent( in) :: pvor(nsh2)
    real(r8kind), intent(out) :: sjacob(nsh2)

    integer      :: i, j, k
    real(r8kind) :: dpsidl(nlat, nlon),  dpsidm(nlat, nlon),  dvordl(nlat, nlon)
    real(r8kind) :: dvordm(nlat, nlon),  gjacob(nlat, nlon),  vv(nsh2)
    real(r8kind) :: azeta(nlat, nlon), dpsidls(nsh2)

    ! space derivatives of potential vorticity 
    call ddl (pvor, vv)
    dvordl = sptogg (vv, pp)
    dvordm = sptogg (pvor, pd)

    ! space derivatives of streamfunction
    call ddl (psiloc, dpsidls)
    dpsidl = sptogg (dpsidls, pp)
    dpsidm = sptogg (psiloc, pd)

    ! jacobian term + orographic forcing
    do j = 1, nlon
      do i = 1, nlat
        gjacob(i, j) = dpsidm(i, j) * (dvordl(i, j) + sinfi(i) * dorodl(i, j)) -  &
   &                   dpsidl(i, j) * (dvordm(i, j) + sinfi(i) * dorodm(i, j))
      enddo
    enddo

    ! dissipation 
    if (lgdiss) then

      !   spatially varying dissipation 
      do k = 1, nsh2
        vv(k) = diss(k, 2) * psiloc(k)
      enddo

      azeta = sptogg (vv, pp)

      do j = 1, nlon
        do i = 1, nlat
          gjacob(i, j) = gjacob(i, j) - dpsidm(i, j) * ddisdy(i, j) &
   &                   - dpsidl(i, j) * ddisdx(i, j)                &
   &                   + rdiss(i, j)  * azeta(i, j)       
        enddo
      enddo

      sjacob = reshape(ggtosp (gjacob), (/nsh2/))

    else

      !   uniform dissipation
      sjacob = reshape(ggtosp (gjacob), (/nsh2/))

      do k = 1, nsh2
        sjacob(k) = sjacob(k) + diss(k, 2) * psi(k, 3)
      enddo

    endif

    ! planetary vorticity advection
    do k = 1, nsh2
      sjacob(k) = sjacob(k) - dpsidls(k)
    enddo

    return

  end subroutine jacobd


  !-----------------------------------------------------------------------
  ! zonal derivative in spectral space
  ! input spectral field as
  ! output spectral field dadl which is as differentiated wrt lambda
  !-----------------------------------------------------------------------
  subroutine ddl (as, dadl)
      
    implicit none

    real(r8kind), intent( in) :: as(nsh, 2)
    real(r8kind), intent(out) :: dadl(nsh, 2)

    integer :: k
 
    do k = 1, nsh
      dadl(k, 1) = -rm(k) * as(k, 2)
      dadl(k, 2) =  rm(k) * as(k, 1)
    enddo
 
    return

  end subroutine ddl


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid
  ! input  spectral field as,  legendre polynomials pploc (pp or pd) 
  !        where pp are legendre polynomials and pd derivatives with
  !        respect to sin(fi)
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg (as, pploc) result(agg)

    implicit none

    real(r8kind), intent( in) :: as(nsh, 2)
    real(r8kind), intent( in) :: pploc(nlat, nsh)

    real(r8kind) :: agg(nlat, nlon)

    integer :: i, ifail, j, k, k1, k2, m, mi, mr, nlon1

   ! inverse legendre transform
    do j = 1, nlon
      do i = 1, nlat
        agg(i, j) = 0.0d0
      enddo
    enddo

    nlon1 = nlon + 1
    k2 = nshm(0)

    do k = 1, k2
      do i = 1, nlat
        agg(i, 1) = agg(i, 1) + as(k, 1) * pploc(i, k)
      enddo
    enddo

    do m = 1, nm
      mr = m + 1
      mi = nlon1 - m
      k1 = k2 + 1
      k2 = k2 + nshm(m)
      do k = k1, k2
        do i = 1, nlat
          agg(i, mr) = agg(i, mr) + as(k, 1) * pploc(i, k)
        enddo
        do i = 1, nlat
          agg(i, mi) = agg(i, mi) - as(k, 2) * pploc(i, k)
        enddo
      enddo
    enddo

    ! inverse fourier transform
    ifail = 0
    call c06fqf (nlat, nlon, agg, 'r', trigi, wgg, ifail)

    return

  end function sptogg
 

  !-----------------------------------------------------------------------
  ! conversion from gaussian grid (agg) to spectral coefficients (as)
  ! input array agg is destroyed
  ! output as contains spectral coefficients
  !-----------------------------------------------------------------------
  function ggtosp (agg_in) result (as)

    implicit none

    real(r8kind), intent(in) :: agg_in(nlat, nlon)
    real(r8kind)             :: as(nsh, 2)

    real(r8kind) :: agg(nlat, nlon)
    integer :: ir, ifail, j, k, k1, k2, m, mi, mr, nlon1, i

    ! Make a local copy of agg_in so it is not destroyed by c06fpf
    agg(:,:) = agg_in(:,:)

    ! fourier transform
    ifail = 0
    call c06fpf (nlat, nlon, agg, 'r', trigd, wgg, ifail)

    ! legendre transform
    do ir = 1, 2
      do k = 1, nsh
        as(k, ir) = 0.0d0
      enddo
    enddo

    nlon1 = nlon + 1

    k2 = nshm(0)

    do k = 1, k2
      do i = 1, nlat
        as(k, 1) = as(k, 1) + agg(i, 1) * pw(i, k)
      enddo
    enddo

    do m = 1, nm
      mr = m + 1
      mi = nlon1 - m
      k1 = k2 + 1
      k2 = k2 + nshm(m)
      do k = k1, k2
        do i = 1, nlat
          as(k, 1) = as(k, 1) + agg(i, mr) * pw(i, k)
          as(k, 2) = as(k, 2) + agg(i, mi) * pw(i, k)
        enddo
      enddo
    enddo

    return

  end function ggtosp
 

  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  subroutine qtopsi

    implicit none

    integer :: k
    real(r8kind) :: r3

    do k = 1, nsh2
      ws(k) = qprime(k, 1) + qprime(k, 3)
      psi(k, 1) = rinhel(k, 1) * (ws(k) + qprime(k, 2))
      psi(k, 2) = ws(k) - 2.d0 * qprime(k, 2)
      psi(k, 3) = qprime(k, 1) - qprime(k, 3)
    enddo

    do k = 1, nsh2
      psit(k, 1) = rinhel(k, 2) * psi(k, 2) + rinhel(k, 3) * psi(k, 3)
      psit(k, 2) = rinhel(k, 4) * psi(k, 2) + rinhel(k, 5) * psi(k, 3)
    enddo

    r3 = 1. / 3.
    do k = 1, nsh2
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
  subroutine psitoq 
      
    implicit none

    integer :: k

    do k = 1, nsh2
      psit(k, 1) = psi(k, 1) - psi(k, 2)
      psit(k, 2) = psi(k, 2) - psi(k, 3)
      qprime(k, 1) = rinhel(k, 0) * psi(k, 1) - rl1 * psit(k, 1)
      qprime(k, 2) = rinhel(k, 0) * psi(k, 2) + rl1 * psit(k, 1) - rl2 * psit(k, 2)
      qprime(k, 3) = rinhel(k, 0) * psi(k, 3) + rl2 * psit(k, 2)
    enddo

    return

  end subroutine psitoq


  !-----------------------------------------------------------------------
  !  computation of potential vorticity qout from stream function sfin
  !-----------------------------------------------------------------------
  subroutine psiq(sfin, qout)

    implicit none

    real(r8kind), intent( in) :: sfin(nsh2, nvl)
    real(r8kind), intent(out) :: qout(nsh2, nvl)

    real(r8kind) :: tus(nsh2)

    integer :: k

    do k = 1, nsh2
      tus(k) = rl1 * sfin(k, 1) - rl1 * sfin(k, 2)
    enddo

    do k = 1, nsh2
      qout(k, 1) = rinhel(k, 0) * sfin(k, 1) - tus(k)
      qout(k, 2) = rinhel(k, 0) * sfin(k, 2) + tus(k)
    enddo

    do k = 1, nsh2
      tus(k) = rl2 * sfin(k, 2) - rl2 * sfin(k, 3)
    enddo

    do k = 1, nsh2
      qout(k, 2) = qout(k, 2) - tus(k)
      qout(k, 3) = rinhel(k, 0) * sfin(k, 3) + tus(k)
    enddo

    return

  end subroutine psiq


  !-----------------------------------------------------------------------
  !  computation of streamfunction bb from potential vorticity qin
  !-----------------------------------------------------------------------
  subroutine qpsi(qin, sfout)

    implicit none

    real(r8kind), intent( in) :: qin(nsh2, nvl)
    real(r8kind), intent(out) :: sfout(nsh2, nvl)

    real(r8kind) :: tus(nsh2, ntl),  r3
    integer :: k

    do k = 1, nsh2
      ws(k) = qin(k, 1) + qin(k, 3)
      sfout(k, 1) = rinhel(k, 1) * (ws(k) + qin(k, 2))
      sfout(k, 2) = ws(k) - 2. * qin(k, 2)
      sfout(k, 3) = qin(k, 1) - qin(k, 3)
    enddo

    do k = 1, nsh2
      tus(k, 1) = rinhel(k, 2) * sfout(k, 2) + rinhel(k, 3) * sfout(k, 3)
      tus(k, 2) = rinhel(k, 4) * sfout(k, 2) + rinhel(k, 5) * sfout(k, 3)
    enddo

    r3 = 1. / 3
    do k = 1, nsh2
      sfout(k, 2) = r3 * (sfout(k, 1) - tus(k, 1) + tus(k, 2))
      sfout(k, 1) = sfout(k, 2) + tus(k, 1)
      sfout(k, 3) = sfout(k, 2) - tus(k, 2)
    enddo

    return

  end subroutine qpsi


  !-----------------------------------------------------------------------
  ! computation of thickness tus from potential vorticity qin
  !-----------------------------------------------------------------------
  subroutine qpsit(qin, tus)

    implicit none

    real(r8kind), intent( in) :: qin(nsh2, nvl)
    real(r8kind), intent(out) :: tus(nsh2, ntl)

    real(r8kind) :: r3, sfout(nsh2, nvl)
    integer :: k

    do k = 1, nsh2
      ws(k) = qin(k, 1) + qin(k, 3)
      sfout(k, 1) = rinhel(k, 1) * (ws(k) + qin(k, 2))
      sfout(k, 2) = ws(k) - 2. * qin(k, 2)
      sfout(k, 3) = qin(k, 1) - qin(k, 3)
    enddo

    do k = 1, nsh2
      tus(k, 1) = rinhel(k, 2) * sfout(k, 2) + rinhel(k, 3) * sfout(k, 3)
      tus(k, 2) = rinhel(k, 4) * sfout(k, 2) + rinhel(k, 5) * sfout(k, 3)
    enddo

    return

  end subroutine qpsit

 
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
  pure function fmtofs (y) result(z)

    implicit none

    real(r8kind), intent( in) :: y(nsh2, nvl)
    real(r8kind)              :: z(nsh2, nvl)

    integer ::  m, n, k, indx, l

    do l = 1, nvl
      k = 1
      do m = 0, nm
        do n = max(m, 1), nm
          k = k + 1
          if (m .eq. 0) then
            indx = n**2
          else
            indx = n**2 + 2 * m - 1
          end if
          z(indx, l) = y(k, l)
          if (m .ne. 0) z(indx + 1, l) = y(k + nsh, l)
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
  pure function fstofm (y, ntr) result(z)

    implicit none

    real(r8kind),intent(in)  :: y(nsh2, nvl)
    integer,     intent(in) :: ntr

    real(r8kind) :: z(nsh2, nvl)

    integer :: m, n, k, indx, i, l

    do l = 1, nvl
      do i = 1, nsh2
        z(i, l) = 0d0
      enddo
      k = 1
      do m = 0, nm
        do n = max(m, 1), nm
          k = k + 1
          if ((m .le. ntr).and.(n .le. ntr)) then
            if (m .eq. 0) then
              indx = n**2
            else
              indx = n**2 + 2 * m - 1
            end if
            z(k, l) = y(indx, l)
            if (m .ne. 0) z(k + nsh, l) = y(indx + 1, l)
          endif
        enddo
      enddo
    enddo

    return

  end function fstofm


  !-----------------------------------------------------------------------
  ! truncates y to ntr and writes to z which is formatted for
  ! lower resolution model Tntr.
  !-----------------------------------------------------------------------
  pure function truncate(y, ntr) result(yt)

    implicit none

    real(r8kind), intent( in) :: y(nsh2, nvl)
    integer,      intent( in) :: ntr
    real(r8kind)              :: yt(nsh2, nvl)

    integer :: m, n, k, indx, l, nshntr, i
    real(r8kind) :: z(nsh2, nvl)

    nshntr = (ntr + 1) * (ntr + 2) * 0.5

    do l = 1, nvl
      k = 1
      do m = 0, nm
        do n = max(m, 1), nm
          k = k + 1
          if (m .eq. 0) then
            indx = n**2
          else
            indx = n**2 + 2 * m - 1
          end if
          z(indx, l) = y(k, l)
          if (m .ne. 0) z(indx + 1, l) = y(k + nsh, l)
        enddo
      enddo
    enddo

    do l = 1, nvl
      do i = 1, nsh2
        yt(i, l) = 0d0
      enddo
      k = 1
      do m = 0, ntr
        do n = max(m, 1), ntr
          k = k + 1
            if (m .eq. 0) then
              indx = n**2
            else
              indx = n**2 + 2 * m - 1
            end if
            yt(k, l) = z(indx, l)
            if (m .ne. 0) yt(k + nshntr, l) = z(indx + 1, l)
        enddo
      enddo
    enddo

    return

  end function truncate


  !-----------------------------------------------------------------------
  ! performs a fourth order runge kutta time step at truncation nm
  ! with time step dt
  ! dqdt calculates the time derivative
  ! input  qprime at current time
  ! output qprime at current time plus dt
  !-----------------------------------------------------------------------
  subroutine forward

    implicit none

    integer :: k, l, nvar
    real(r8kind) :: dt2, dt6
    real(r8kind) :: y(nsh2, nvl), dydt(nsh2, nvl), yt(nsh2, nvl)
    real(r8kind) :: dyt(nsh2, nvl), dym(nsh2, nvl)

    nvar = (nm + 2) * nm
    dt2 = dtt * 0.5d0
    dt6 = dtt / 6d0
    y = fmtofs(qprime)
    call dqdt(y, dydt)
    do l = 1, nvl
      do k = 1, nvar
        yt(k, l) = y(k, l) + dt2 * dydt(k, l)
      enddo
    enddo
    call dqdt(yt, dyt)
    do l = 1, nvl
      do k = 1, nvar
        yt(k, l) = y(k, l) + dt2 * dyt(k, l)
      enddo
    enddo
    call dqdt(yt, dym)
    do l = 1, nvl
      do k = 1, nvar
        yt(k, l) = y(k, l) + dtt * dym(k, l)
        dym(k, l) = dyt(k, l) + dym(k, l)
      enddo
    enddo
    call dqdt(yt, dyt)
    do l = 1, nvl
      do k = 1, nvar
        y(k, l) = y(k, l) + dt6 * (dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
      enddo
    enddo
    qprime = fstofm(y, nm)

    return

  end subroutine forward


  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  subroutine dqdt(y, dydt)

    implicit none

    real(r8kind), intent( in) :: y(nsh2, nvl)
    real(r8kind), intent(out) :: dydt(nsh2, nvl)

    qprime = fstofm(y, nm)
    call qtopsi
    call ddt
    dydt = fmtofs(dqprdt)

    return

  end subroutine dqdt


  !-----------------------------------------------------------------------
  ! computation of geostrophic winds at all levels
  ! computes geopotential height in [m2 / s2[ = [gz] from streamfunction 
  ! by solving the linear balance equation: 
  ! del phi = (1 - mu**2 ) d psi / dmu + mu del psi
  ! the global mean value is not determined and set to zero
  !  
  !-----------------------------------------------------------------------
  subroutine gridfields

    implicit none

    integer :: i, j, k, l
    real(r8kind) :: facwind, facsf, facgp, facpv
    real(r8kind) :: dpsdl(nlat, nlon), dpsdm(nlat, nlon), psik(nsh2), vv(nsh2)
    real(r8kind) :: fmu(nlat)
    real(r8kind) :: delpsis(nsh2), delpsig(nlat, nlon)
    real(r8kind) :: dmupsig(nlat, nlon), delgeog(nlat, nlon)
    real(r8kind) :: delgeos(nsh2), geos(nsh2)

    call qtopsi      

    ! space derivatives of streamfunction
    facwind = radius * om
    facsf = om * (radius)**2
    facgp = (om**2) * (radius**2)
    facpv = om

    do i = 1, nlat
      fmu(i) = 1 - sin(pi * phi(i) / 180d0)**2
    enddo

    do l = 1, nvl
      psig(:, :, l) = sptogg(psi(1, l), pp)
      do j = 1, nlon
        do i = 1, nlat
          psig(i, j, l) = facsf * psig(i, j, l)
        enddo
      enddo

      qgpv(:, :, l) = sptogg(qprime(1, l), pp)
      do j = 1, nlon
        do i = 1, nlat
          qgpv(i, j, l) = facpv * qgpv(i, j, l)
        enddo
      enddo

      do k = 1, nsh2
        psik(k) = psi(k, l)
      enddo


      call ddl (psik, vv)
      dpsdl = sptogg (vv, pp)
      dpsdm = sptogg (psik, pd)

      do j = 1, nlon
        do i = 1, nlat
          ug(i, j, l) = -facwind * dpsdm(i, j) * cosfi(i)
          vg(i, j, l) = +facwind * dpsdl(i, j) / cosfi(i)
        enddo
      enddo

      ! solve linear balance equation
      call lap(psi(1, l), delpsis)
      delpsig = sptogg(delpsis, pp)
      dmupsig = sptogg(psi(1, l), pd)

      do j = 1, nlon
        do i = 1, nlat
          delgeog(i, j) = fmu(i) * dmupsig(i, j) + sinfi(i) * delpsig(i, j)
        enddo
      enddo

      delgeos = reshape(ggtosp(delgeog), (/nsh2/))
      call lapinv(delgeos, geos)
      geos(1) = 0.d0
      geopg(:, :, l) = sptogg(geos, pp)


      do j = 1, nlon
        do i = 1, nlat
          geopg(i, j, l) = facgp * geopg(i, j, l)
        enddo
      enddo

    enddo

    return

  end subroutine gridfields

      
  !-----------------------------------------------------------------------
  ! computation of laplace operator in spectral domain
  ! input  xs  field in spectral form
  ! output xsl laplace of xs in spectral form
  !-----------------------------------------------------------------------
  subroutine lap(xs, xsl)

    implicit none

    real(r8kind), intent( in) :: xs(nsh2)
    real(r8kind), intent(out) :: xsl(nsh2)

    integer :: k

    do k = 1, nsh2
      xsl(k) = xs(k) * rinhel(k, 0)
    enddo

    return

  end subroutine lap

      
  !-----------------------------------------------------------------------
  ! computation of laplace operator in spectral domain
  ! input  xsl field in spectral form
  ! output xs  inverse laplace of xs in spectral form
  !-----------------------------------------------------------------------
  subroutine lapinv(xsl, xs)

    implicit none

    real(r8kind), intent( in) :: xsl(nsh2)
    real(r8kind), intent(out) :: xs(nsh2)

    integer :: k

    do k = 1, nsh2
      xs(k) = xsl(k) * rinhel(k, 1)
    enddo

    return

  end subroutine lapinv


  !-----------------------------------------------------------------------
  ! computation of artifical forcing according to roads(1987)
  ! the forcing is computed from daily ecmwf data for the winter season
  ! the forcing is composed of a climatological forcing and a contribution
  ! of the eddy terms.
  !------------------------------------------------------------------------
  subroutine artiforc_iter

    implicit none

    integer :: i, j, k, l, iy, id, iday, nyb, nye, nd
    real(r4kind) :: psi4(nsh2, 3)
    real(r8kind) :: psim(nsh2, 3), sum(nsh2, 3)
    real(r8kind) :: eddf(nsh2, 3), totf(nsh2, 3), climf(nsh2, 3), forg(nlat, nlon)

    open(unit = 46, file = './qgmodelT42.T21', status = 'old', form = 'formatted')
    open(unit = 32, file = './qgpvforT42world.grads', form = 'unformatted')

    nyb = 83
    nye = 92
    nday = 91

    do l = 1, 3
      do k = 1, nsh2
        sum(k, l) = 0.
      enddo
    enddo

    iday = 0

    do iy = nyb, nye

      nd = 90
      if (iy .eq. 84 .or. iy .eq. 88 .or. iy .eq. 92) nd = 91

      do id = 1, nd
        read(46, *)i
        read(46, '(5e12.5)')((psi4(k, l), k = 1, nsh2), l = 1, 3)
        iday = iday + 1
        do l = 1, 3
          do k = 1, nsh2
            sum(k, l) = sum(k, l) + psi4(k, l)
          enddo
        enddo
      enddo
    enddo


    do l = 1, 3
      do k = 1, nsh2
        psim(k, l) = sum(k, l) / iday
      enddo
    enddo

    ! calculate the climatological forcing
    do l = 1, 3
      do k = 1, nsh2
        psi(k, l) = psim(k, l)
      enddo
    enddo

    call psitoq

    call ddt 

    do l = 1, 3
      do k = 1, nsh2
        climf(k, l) = -dqprdt(k, l)
      enddo
    enddo

    ! calculate the eddy forcing
    rewind(46)

    do l = 1, 3
      do k = 1, nsh2
        sum(k, l) = 0.
      enddo
    enddo
    iday = 0

    do iy = nyb, nye

      nd = 90
      if (iy .eq. 84 .or. iy .eq. 88 .or. iy .eq. 92) nd = 91

      do id = 1, nday
        read(46, *)i
        read(46, '(5e12.5)')((psi4(k, l), k = 1, nsh2), l = 1, 3)
        iday = iday + 1
        do l = 1, 3
          do k = 1, nsh2
            psi(k, l) = psi4(k, l) - psim(k, l)
          enddo
        enddo
        call psitoq
        call eddforc
        do l = 1, 3
          do k = 1, nsh2
            sum(k, l) = sum(k, l) + dqprdt(k, l)
          enddo
        enddo
      enddo
    enddo

    do l = 1, 3
      do k = 1, nsh2
        eddf(k, l) = -sum(k, l) / iday
      enddo
    enddo

    ! compute the total forcing
    do l = 1, 3
      do k = 1, nsh2
        for(k, l) = climf(k, l) + eddf(k, l)
      enddo
    enddo

    write(14, '(1E12.5)') ((for(k, l), k = 1, nsh2), l = 1, 3)

    do l = 1, 3
      forg = sptogg(for(1, l), pp)
      write(32) ((real(forg(i, j)), j = 1, nlon), i = 1, nlat)
    enddo

    close(46)
    close(32)

    return

  end subroutine artiforc_iter


  !-----------------------------------------------------------------------
  ! computation of the eddy forcing
  !-----------------------------------------------------------------------
  subroutine eddforc

    implicit none

    call jacobedd (psi(1, 1), qprime(1, 1), dqprdt(1, 1))

    call jacobedd (psi(1, 2), qprime(1, 2), dqprdt(1, 2))

    call jacobedd (psi(1, 3), qprime(1, 3), dqprdt(1, 3))

    return

  end subroutine eddforc
      

  !-----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  ! the only difference with the routine jacob 
  ! is that the planetary vorticity advection is omitted.
  !-----------------------------------------------------------------------
  subroutine jacobedd (psiloc, pvor, sjacob)

    implicit none

    real(r8kind), intent( in) :: psiloc(nsh2)
    real(r8kind), intent( in) :: pvor(nsh2)
    real(r8kind), intent(out) :: sjacob(nsh2)

    real(r8kind) :: vv(nsh2)
    integer :: i, j, k
    real(r8kind) :: dpsidls(nsh2)
    real(r8kind) :: dpsidl(ngp), dpsidm(ngp), dvordl(ngp), dvordm(ngp)
    real(r8kind) :: gjacob(ngp)

   ! space derivatives of potential vorticity
    call ddl (pvor, vv)
    dvordl = reshape(sptogg (vv, pp), (/ngp/))
    dvordm = reshape(sptogg (pvor, pd), (/ngp/))

    ! space derivatives of streamfunction
    call ddl (psiloc, dpsidls)
    dpsidl = reshape(sptogg (dpsidls, pp), (/ngp/))
    dpsidm = reshape(sptogg (psiloc, pd), (/ngp/))

    ! jacobian term
    do j = 1, ngp
        gjacob(j) = dpsidm(j) * dvordl(j) - dpsidl(j) * dvordm(j)
    enddo

    sjacob = reshape(ggtosp (gjacob), (/nsh2/))

    return

  end subroutine jacobedd


  !-----------------------------------------------------------------------
  ! computation of artifical forcing according to roads(1987)
  ! the forcing is computed from file obsfile
  !------------------------------------------------------------------------
  subroutine artiforc

    implicit none

    integer :: i, j, k, l, iday, fl, nvar
    real(r4kind) :: psi4(nsh2, 3)
    real(r8kind) :: sum(nsh2, 3), forg(nlat, nlon), dlon, psifs(nsh2, 3)

    nvar = (nm + 2) * nm

    dlon = 360d0 / real(nlon)

    open(14, file = 'qgpvforT' // trim(ft) // '.dat', form = 'formatted')

    open(unit = 46, file = './' // obsfile, status = 'old', form = 'unformatted')
    fl = index(obsfile, " ") - 1
    open(unit = 32, file = 'qgpvforT' // trim(ft) // '.grads', form = 'unformatted')

    open(unit = 99, file = './' // obsfile(1:fl) // trim(ft) // '.grads', form = 'unformatted')

    write(*, '(A, A)') "Calculating forcing from ", obsfile

    do l = 1, nvl
      do k = 1, nsh2
        sum(k, l) = 0.
        for(k, l) = 0d0
      enddo
    enddo

    ! calculate the mean tendency
    iday = 0

 10 continue
    do l = 1, nvl
      read(46, end = 20) (psi4(k, nvl - l + 1), k = 1, nvar)
    enddo

    iday = iday + 1
    do l = nvl, 1, -1
      do k = 1, nvar
        psifs(k, l) = psi4(k, l)
      enddo
    enddo
    psi = fstofm(psifs, nm)
    do l = nvl, 1, -1
      forg = sptogg(psi(1, l), pp)
      write(99) ((real(forg(j, i)), i = 1, nlon), j = 1, nlat)
    enddo
    call psitoq
    call ddt
    do l = 1, nvl
      do k = 1, nsh2
        sum(k, l) = sum(k, l) + dqprdt(k, l)
      enddo
    enddo
    goto 10
 20 continue
    close(99)

    do l = 1, nvl
      do k = 1, nsh2
        for(k, l) = -sum(k, l) / real(iday)
      enddo
    enddo

    write(14, '(1E12.5)') ((for(k, l), k = 1, nsh2), l = 1, nvl)

    do l = nvl, 1, -1
      forg = sptogg(for(1, l), pp)
      write(32) ((real(forg(i, j)), j = 1, nlon), i = 1, nlat)
    enddo

    close(46)
    close(32)

    open(50, file = './' // obsfile(1:fl) // trim(ft) // '.ctl', form = 'formatted')
    write(50, '(A)') 'dset ^' // obsfile(1:fl) // trim(ft) // '.grads'
    write(50, '(A)') 'undef 9.99e+10'
    write(50, '(A)') 'options sequential big_endian'
    write(50, '(A)') 'title three level QG model'
    write(50, '(A)') '*'
    write(50, '(A, i4, A, F19.14)') 'xdef ', nlon, ' linear  0.000 ', dlon
    write(50, '(A)') '*'
    write(50, '(A, I4, A, 1F19.14)') 'ydef ', nlat, ' levels ', phi(1)
    write(50, '(F19.14)') (phi(j), j = 2, nlat)
    write(50, '(A)') '*'
    write(50, '(A)') 'zdef  3 levels 800 500 200'
    write(50, '(A)') '*'
    write(50, '(A, I6, A)') 'tdef ', iday, ' linear 1jan0001 1dy'
    write(50, '(A)') '*'
    write(50, '(A)') 'vars  1'
    write(50, '(A)') 'sf    3  99 streamfunction (nondim)'
    write(50, '(A)') 'endvars'

    close(50)

    write(*, '(A, I6)') "Number of states used to calculate forcing: ", iday
    write(*, '(A)') "Forcing saved in files: "
    write(*, '(A)') 'qgpvforT' // trim(ft) // '.grads'
    write(*, '(A)') 'qgpvforT' // trim(ft) // '.dat'

    return

  end subroutine artiforc


  !-----------------------------------------------------------------------
  ! convert ecmwf winter analyses files to asc file to be read by
  ! artiforc
  !-----------------------------------------------------------------------
  subroutine analyses

    implicit none

    integer :: i, j, k, nyb, nye, npw, iy, id, l, iday, irec, nsh2ntr

    real(r8kind) :: psiloc(nsh2, 3), psiT21(nsh2, 3)
    real(r4kind) :: psi4(nsh2, 3), psig4(nlat, nlon, nvl), scalesf
    character(len=4) :: fy

    scalesf = 1d0 / (radius * radius * om)
    nsh2ntr = 22 * 23

    open(98, file = './anwin79_10T42.dat', form = 'formatted')
    open(96, file = './anwin79_10T21.dat', form = 'formatted')
  ! open(unit = 97, file = './eracheck.grads', form = 'unformatted')

    nyb = 1979
    nye = 2010
    npw = 90

    write(fy, '(I4.4)') iy

    do iy = nyb, nye
      write(fy, '(I4.4)') iy
      open(unit = 99, file = './era' // fy // 'djf.grads',  &
   &       form = 'unformatted', access = 'direct', recl = nlat * nlon * 4)

      do id = 1, npw
        do l = nvl, 1, -1
          irec = (id - 1) * 3 + nvl - l + 1
          read(99, rec = irec) ((psig4(j, i, l), i = 1, nlon), j = 1, nlat)
        enddo
        do l = nvl, 1, -1
          do i = 1, nlon
            do j = 1, nlat
              psig(j, i, l) = psig4(j, i, l) * scalesf
            enddo
          enddo
  !       write(97) ((real(psig(j, i, l)), i = 1, nlon), j = 1, nlat)
        enddo
        do l = 1, nvl
          psiloc(:, l) = reshape(ggtosp(psig(:, :, l)), (/nsh2/))
        enddo
        write(98, *) id
        write(98, '(5e12.5)')((real(psiloc(k, l)), k = 1, nsh2), l = 1, nvl)
        psiT21 = truncate(psiloc, 21)
        write(96, *) id
        write(96, '(5e12.5)')((real(psiT21(k, l)), k = 1, nsh2ntr), l = 1, nvl)
      enddo        
      close(99)
    enddo
    close(98)
  ! close(97)

  end subroutine analyses


  !-----------------------------------------------------------------------
  ! output streamfunction data to outputfile
  !-----------------------------------------------------------------------
  subroutine diagsf(istep)

    implicit none

    integer, intent(in) :: istep

    integer :: i, j, k, nout

    real(r8kind) :: psiloc(nsh2),  pvor(nsh2),  sjacob(nsh2), dlon
    real(r4kind) :: psi4(nsh2, 3)

    nout = nday * nstepsperday / nstepsbetweenoutput + 1

    dlon = 360d0 / real(nlon)

    if (istep .eq. 0) then
      open(52, file = 'qgmodelsfT' // trim(ft) // '.ctl', form = 'formatted')
      write(52, '(A)') 'dset ^qgmodelsfT' // trim(ft) // '.grads'
      write(52, '(A)') 'undef 9.99e+10'
      write(52, '(A)') 'options sequential big_endian'
      write(52, '(A)') 'title T' // trim(ft)
      write(52, '(A)') '*'
      write(52, '(A, I6, A, F19.14)') 'xdef ', nlon, ' linear  0.000 ', dlon
      write(52, '(A)') '*'
      write(52, '(A, I4, A, 1F19.14)') 'ydef ', nlat, ' levels ', phi(1)
      write(52, '(F19.14)') (phi(j), j = 2, nlat)
      write(52, '(A)') '*'
      write(52, '(A)') 'zdef  3 levels 800 500 200'
      write(52, '(A)') '*'
      write(52, '(A, I6, A)') 'tdef ', nout, ' linear 1jan0001 1dy'
      write(52, '(A)') '*'
      write(52, '(A)') 'vars  1'
      write(52, '(A)') 'psi    3  99 streamfunction [m2 / s]'
      write(52, '(A)') 'endvars'

      close(52)
      open(52, file = 'qgmodelsfT' // trim(ft) // '.grads', form = 'unformatted')
    endif

    if (mod(istep, nstepsbetweenoutput) .eq. 0) then
      call gridfields
      do k = nvl, 1, -1
        write(52) ((real(psig(j, i, k)), i = 1, nlon), j = 1, nlat)
      enddo
    endif

  end subroutine diagsf


  !-----------------------------------------------------------------------
  ! output model data to outputfile
  !-----------------------------------------------------------------------
  subroutine diag(istep)

    implicit none

    integer, intent(in) :: istep

    integer :: i, j, k, nout

    real(r8kind) :: psiloc(nsh2),  pvor(nsh2),  sjacob(nsh2),  dlon
    real(r4kind) :: psi4(nsh2, 3)

    nout = nday * nstepsperday / nstepsbetweenoutput + 1
    dlon = 360d0 / real(nlon)

    if (istep .eq. 0) then
      open(50, file = 'qgmodelT' // trim(ft) // '.ctl', form = 'formatted')
      write(50, '(A)') 'dset ^qgmodelT' // trim(ft) // '.grads'
      write(50, '(A)') 'undef 9.99e+10'
      write(50, '(A)') 'options sequential big_endian'
      write(50, '(A)') 'title T' // trim(ft)
      write(50, '(A)') '*'
      write(50, '(A, i4, A, F19.14)') 'xdef ', nlon, ' linear  0.000 ', dlon
      write(50, '(A)') '*'
      write(50, '(A, I4, A, 1F19.14)') 'ydef ', nlat, ' levels ', phi(1)
      write(50, '(F19.14)') (phi(j), j = 2, nlat)
      write(50, '(A)') '*'
      write(50, '(A)') 'zdef  3 levels 800 500 200'
      write(50, '(A)') '*'
      write(50, '(A, I6, A)') 'tdef ', nout, ' linear 1jan0001 1dy'
      write(50, '(A)') '*'
      write(50, '(A)') 'vars   5'
      write(50, '(A)') 'geopg  3  99 geopotential [m2 / s2]'
      write(50, '(A)') 'psi    3  99 streamfunction [m2 / s]'
      write(50, '(A)') 'pv     3  99 pv [1 / s]'
      write(50, '(A)') 'u      3  99 zonal velocity [m / s]'
      write(50, '(A)') 'v      3  99 meridional velocity [m / s]'
      write(50, '(A)') 'endvars'

      close(50)
      open(50, file = 'qgmodelT' // trim(ft) // '.grads', form = 'unformatted')
    endif

    if (mod(istep, nstepsbetweenoutput) .eq. 0) then
      call gridfields
      do k = nvl, 1, -1
        write(50) ((real(geopg(j, i, k)), i = 1, nlon), j = 1, nlat)
      enddo
      do k = nvl, 1, -1
        write(50) ((real(psig(j, i, k)), i = 1, nlon), j = 1, nlat)
      enddo
      do k = nvl, 1, -1
        write(50) ((real(qgpv(j, i, k)), i = 1, nlon), j = 1, nlat)
      enddo
      do k = nvl, 1, -1
        write(50) ((real(ug(j, i, k)), i = 1, nlon), j = 1, nlat)
      enddo
      do k = nvl, 1, -1
        write(50) ((real(vg(j, i, k)), i = 1, nlon), j = 1, nlat)
      enddo
    endif

  end subroutine diag
      
      
  !-----------------------------------------------------------------------
  ! output streamfunction state that can be read as initial state
  !-----------------------------------------------------------------------
  subroutine writestate

    implicit none

    integer :: i, j, k, l

    open(12, file = 'qgendT' // trim(ft) // '.dat', form = 'formatted')
    do l = 1, 3
      do k = 1, nsh2
         write(12, *) psi(k, l)
      enddo
    enddo
    close(12)

  end subroutine writestate

end module QG
