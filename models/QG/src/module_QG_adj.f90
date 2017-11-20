module QG_Model_ADJ

  use kind
  use Abstract_Model, only : abstract_model_type
  use QG_Model_TL
  use QG_Config
  use QG_GGSP
  use QG_Util

  implicit none

  private

  public :: qg_adj_type

  type, extends(abstract_model_type) :: qg_adj_type
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
    integer :: nmat           ! ? Needed for adjoint

    ! Current model step and simulaiton clock
    integer      :: step
    real(r8kind) :: clock

    ! Model time step
    real(r8kind) :: dtt                      ! dimensionless time step

    ! Model state
    real(r8kind), allocatable :: psi(:,:)        ! Stream function at the nvl levels
    real(r8kind), allocatable :: trajectory(:,:) ! Adjoint trajectory
    real(r8kind), allocatable :: psit(:,:)       ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprime(:,:)     ! Potential vorticity
    real(r8kind), allocatable :: psitb(:,:)       ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprimeb(:,:)     ! Potential vorticity

    ! Model Forcing
    real(r8kind), allocatable :: for(:,:)    ! Constant potential vorticity forcing at the nvl levels

    ! Spectral Coefficients
    integer,      allocatable :: nshm(:)     ! Contains numbers 22 down to 1 for index 0 to 21
    integer,      allocatable :: ll(:)       ! Contains total wavenumber n of each spherical harmonic of the corresponding index
    real(r8kind), allocatable :: pp(:,:)     ! Legendre polynomials defined at Gausian latitudes
    real(r8kind), allocatable :: pd(:,:)     ! Mu derivative of Legendre polynomials
    real(r8kind), allocatable :: pw(:,:)     ! Weights for Legendre integrals
    real(r8kind), allocatable :: derimu(:,:) ! Mu operator needed for adjoint

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
    final :: destructor_qg_adj
    procedure :: adv_nsteps
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
    procedure :: get_psig
    procedure :: get_lat_lon_grid
    procedure :: get_for
    procedure :: get_state_vector
    procedure :: get_location_vector
    procedure :: get_interpolation_weights
    procedure, private :: init_state
    procedure, private :: dqdt
    procedure, private :: dqdt_b
    procedure, private :: ddt
    procedure, private :: ddt_b
    procedure, private :: jacob
    procedure, private :: jacob_b
    procedure, private :: jacobd
    procedure, private :: jacobd_b
    procedure, private :: psitoq
    procedure, private :: qtopsi
    procedure, private :: qtopsi_b
    procedure, private :: lap
    procedure, private :: lapinv
    procedure, private :: fmtofs
    procedure, private :: fmtofs_b
    procedure, private :: fstofm
    procedure, private :: fstofm_b
  end type qg_adj_type

  interface qg_adj_type
    procedure :: constructor_qg_adj
  end interface

  ! Mathematical and physical constants
  real(r8kind), parameter :: pi = 4d0 * atan(1d0)     ! value of pi


contains


  !-------------------------------------------------------------------------------
  ! constructor_qg_adj
  !-------------------------------------------------------------------------------
  function constructor_qg_adj(config, state, state_vector, trajectory, trajectory_vector, for, step) result (qg_adj)

    type(qg_config_type),   intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:,:)
    real(r8kind), optional, intent(in) :: state_vector(:)
    real(r8kind), optional, intent(in) :: trajectory(:,:)
    real(r8kind), optional, intent(in) :: trajectory_vector(:)
    real(r8kind), optional, intent(in) :: for(:,:)
    integer,      optional, intent(in) :: step

    type(qg_adj_type)                :: qg_adj
    real(r8kind), allocatable          :: psig3d(:,:,:), psisp(:,:)
    real(r8kind), allocatable          :: trajg3d(:,:,:)

    real(r8kind) :: nsteps_per_day

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    ! Set the model config
    qg_adj%config = config

    ! Get model grid dimensions
    qg_adj%ft = config%get_ft()
    qg_adj%nm = config%get_nm()
    qg_adj%nlon = config%get_nlon()
    qg_adj%nlat = config%get_nlat()
    qg_adj%nvl = config%get_nvl()
    qg_adj%ntl = config%get_ntl()
    qg_adj%nsh = config%get_nsh()
    qg_adj%nsh2 = config%get_nsh2()

    ! Get model grid converter
    qg_adj%ggsp = config%get_ggsp()

    ! Get model forcing from configuration
    allocate(qg_adj%for, source = config%get_for())

    ! Get Nondimensional relaxation coefficients
    qg_adj%relt1 = config%get_relt1()
    qg_adj%relt2 = config%get_relt2()

    ! Get dissipation coefficients for each spherical harmonic
    allocate(qg_adj%diss, source = config%get_diss())

    ! Get gauss points in radians and sine and cosine of it
    allocate(qg_adj%phi, source = config%get_phi())
    allocate(qg_adj%sinfi, source = config%get_sinfi())
    allocate(qg_adj%cosfi, source = config%get_cosfi())

    ! Get derivatives of orog
    allocate(qg_adj%dorodl, source = config%get_dorodl())
    allocate(qg_adj%dorodm, source = config%get_dorodm())

    ! Get orography and land-sea mask dependent friction option
    qg_adj%lgdiss = config%get_lgdiss()

    ! Get landsea-mask/orography dependent friction
    allocate(qg_adj%rdiss, source = config%get_rdiss())
    allocate(qg_adj%ddisdx, source = config%get_ddisdx())
    allocate(qg_adj%ddisdy, source = config%get_ddisdy())

    ! Get Laplace and Helmholtz operator for Q-PSI inversion
    allocate(qg_adj%rinhel(qg_adj%nsh2,0:5))
    qg_adj%rinhel = config%get_rinhel()

    ! Get one over Rossby rad. of def. squared of 200-500, and 500-800 thicknesses
    qg_adj%rl1 = config%get_rl1()
    qg_adj%rl2 = config%get_rl2()

    ! Initialize model step and clock
    if (present(step)) then
      qg_adj%step = step
    else
      qg_adj%step = 0
    end if
    qg_adj%clock = qg_adj%step * config%get_time_step()

    ! Initialize time step of the model:
    nsteps_per_day = 24.0d0 * 3600.0d0 / real(config%get_time_step())
    qg_adj%dtt = (1d0 / nsteps_per_day) * pi * 4d0

    ! Initialize streamfunction
    if (present(state)) then
      call qg_adj%init_state(psi=state)
    else if (present(state_vector)) then
      allocate(psig3d(qg_adj%nlon, qg_adj%nlat, qg_adj%nvl))
      allocate(psisp(qg_adj%nsh2, qg_adj%nvl))
      psig3d = reshape(state_vector,(/qg_adj%nlon, qg_adj%nlat, qg_adj%nvl/))
      psig3d(:,:,:) = psig3d(:,:,:) / facsf
      psisp(:,1) = reshape(qg_adj%ggsp%ggtosp(transpose(psig3d(:,:,1))), (/qg_adj%nsh2/))
      psisp(:,2) = reshape(qg_adj%ggsp%ggtosp(transpose(psig3d(:,:,2))), (/qg_adj%nsh2/))
      psisp(:,3) = reshape(qg_adj%ggsp%ggtosp(transpose(psig3d(:,:,3))), (/qg_adj%nsh2/))
      call qg_adj%init_state(psi=psisp)
    else
      call qg_adj%init_state()
    end if

    ! Initialize trajectory
    if (present(trajectory)) then
      allocate(qg_adj%trajectory, source = trajectory)
    else if (present(trajectory_vector)) then
      allocate(trajg3d(qg_adj%nlon, qg_adj%nlat, qg_adj%nvl))
      allocate(qg_adj%trajectory(qg_adj%nsh2,qg_adj%nvl))
      trajg3d = reshape(trajectory_vector,(/qg_adj%nlon, qg_adj%nlat, qg_adj%nvl/))
      trajg3d(:,:,:) = trajg3d(:,:,:) / facsf
      qg_adj%trajectory(:,1) = reshape(qg_adj%ggsp%ggtosp(transpose(trajg3d(:,:,1))), (/qg_adj%nsh2/))
      qg_adj%trajectory(:,2) = reshape(qg_adj%ggsp%ggtosp(transpose(trajg3d(:,:,2))), (/qg_adj%nsh2/))
      qg_adj%trajectory(:,3) = reshape(qg_adj%ggsp%ggtosp(transpose(trajg3d(:,:,3))), (/qg_adj%nsh2/))
    else
      allocate(qg_adj%trajectory, source = qg_adj%psi)
    end if
    allocate(qg_adj%psitb(qg_adj%nsh2,qg_adj%ntl))
    allocate(qg_adj%qprimeb(qg_adj%nsh2,qg_adj%nvl))
    call qg_adj%psitoq(qg_adj%trajectory, qg_adj%psitb, qg_adj%qprimeb)

  end function constructor_qg_adj


  !------------------------------------------------------------------
  ! destructor_qg_adj
  !
  ! Deallocates pointers used by a qg_adj_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_qg_adj(this)


    type(qg_adj_type), intent(inout) :: this

    ! No pointers in qg_adj_type object so we do nothing

  end subroutine destructor_qg_adj


  !-------------------------------------------------------------------------------
  ! init_state
  !-------------------------------------------------------------------------------
  subroutine init_state(this, psi)

    use netcdf
    use netcdf_utilities

    class(qg_adj_type),   intent(inout) :: this
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

    class(qg_adj_type) :: this
    integer           :: nsteps

    integer :: step, k, l, nvar, i
    real(r8kind) :: dt2, dt6
    real(r8kind) :: y(this%nsh2, this%nvl), dydt(this%nsh2, this%nvl), yt(this%nsh2, this%nvl)
    real(r8kind) :: yt1(this%nsh2, this%nvl), yt2(this%nsh2, this%nvl), y1(this%nsh2, this%nvl)
    real(r8kind) :: dyt(this%nsh2, this%nvl), dym(this%nsh2, this%nvl)
    real(r8kind) :: yb(this%nsh2, this%nvl), dydtb(this%nsh2, this%nvl), ytb(this%nsh2, this%nvl)
    real(r8kind) :: dytb(this%nsh2, this%nvl), dymb(this%nsh2, this%nvl)
    real(r8kind) :: qprimeb_save(this%nsh2, this%nvl)
    real(r8kind) :: tempb
    type(qg_config_type) :: config
    type(qg_tl_type) :: qg_tl

    config = this%config

    if (nsteps > 0) then

      nvar = (this%nm + 2) * this%nm
      dt2 = this%dtt * 0.5d0
      dt6 = this%dtt / 6d0

      ! Advance the model forward in time n steps
      do step = 1, nsteps

        ! Initialize a QG TL model with the current state and trajectory
        qg_tl = qg_tl_type(this%config, state=this%psi, trajectory=this%trajectory)

        ! Advance the QG TL model one step
        call qg_tl%adv_nsteps(1)

        ! Save the current trajectory
        qprimeb_save = this%qprimeb

        ! Set the adjoint trajectory to the initial perturbation
        this%qprimeb = qg_tl%get_qprimed() - this%qprimeb

        ! Move the perturbation backward
        y = this%fmtofs(this%qprime)
        call this%DQDT(y, dydt)
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + dt2 * dydt(k, l)
          end do
        end do
        call this%DQDT(yt, dyt)
        yt1 = yt
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + dt2 * dyt(k, l)
          end do
        end do
        call this%DQDT(yt, dym)
        yt2 = yt
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + this%dtt * dym(k, l)
            dym(k, l) = dyt(k, l) + dym(k, l)
          end do
        end do
        call this%DQDT(yt, dyt)
        y1 = y
        do l = 1, this%nvl
          do k = 1, nvar
            y(k, l) = y(k, l) + dt6 * (dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
          end do
        end do
        yb = 0.0_8

       call this%fstofm_B(y, yb, this%nm, this%qprimeb)
!        yb = this%fmtofs(this%qprimeb)
!        this%qprimeb = 0.0_r8kind

        dymb = 0.0_8
        dytb = 0.0_8
        dydtb = 0.0_8
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            tempb = dt6 * yb(k, l)
            dydtb(k, l) = dydtb(k, l) + tempb
            dytb(k, l) = dytb(k, l) + tempb
            dymb(k, l) = dymb(k, l) + 2. * tempb
          end do
        end do
        ytb = 0.0_8

        call this%DQDT_B(yt, ytb, dyt, dytb)
        dytb = 0.0_8
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            dytb(k, l) = dytb(k, l) + dymb(k, l)
            yb(k, l) = yb(k, l) + ytb(k, l)
            dymb(k, l) = dymb(k, l) + this%dtt * ytb(k, l)
            ytb(k, l) = 0.0_8
          end do
        end do
        call this%DQDT_B(yt2, ytb, dym, dymb)
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            yb(k, l) = yb(k, l) + ytb(k, l)
            dytb(k, l) = dytb(k, l) + dt2 * ytb(k, l)
            ytb(k, l) = 0.0_8
          end do
        end do
        call this%DQDT_B(yt1, ytb, dyt, dytb)
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            yb(k, l) = yb(k, l) + ytb(k, l)
            dydtb(k, l) = dydtb(k, l) + dt2 * ytb(k, l)
            ytb(k, l) = 0.0_8
          end do
        end do
        call this%DQDT_B(y1, yb, dydt, dydtb)

        call this%fmtofs_B(this%qprime, this%qprimeb, yb)
!       this%qprimeb = this%fstofm(yb, this%nm)
!       yb = 0.0_r8kind

        ! Update the trajectory
        this%qprimeb = qprimeb_save - this%qprimeb

        ! Update model state with original trajectory plus its increment
        this%qprime = this%qprime - this%qprimeb

        ! Make stream function consistent with potential vorticity
        call this%qtopsi(this%qprimeb, this%trajectory, this%psitb)
        call this%qtopsi(this%qprime, this%psi, this%psit)

        ! Inrement the step count
        this%step = this%step - 1

      end do

    end if

  end subroutine adv_nsteps


  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  subroutine dqdt(this, y, dydt)

    class(qg_adj_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: y(:,:)
    real(r8kind),         intent(  out) :: dydt(:,:)

! qprime
    real(r8kind) :: local_qprime(this%nsh2, this%nvl)
! psi
    real(r8kind) :: local_psi(this%nsh2, this%nvl)
! psit
    real(r8kind) :: local_psit(this%nsh2, this%ntl)

    real(r8kind) :: dqprdt(this%nsh2,this%nvl) ! time derivative of qprime

    local_qprime = this%fstofm(y, this%nm)
    call this%qtopsi(local_qprime, local_psi, local_psit)            ! qprime --> psi and psit
    dqprdt = this%ddt(local_psi, local_psit, local_qprime, this%for) ! psi, psit, qprime, for, diss --> dqprdt
    dydt = this%fmtofs(dqprdt)

    return

  end subroutine dqdt


  !-----------------------------------------------------------------------
  !  Differentiation of dqdt in reverse (adjoint) mode:
  !   gradient     of useful results: tmp y dydt
  !   with respect to varying inputs: tmp y
  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  subroutine DQDT_B(this, y, yb, dydt, dydtb)

    class(qg_adj_type), intent(inout) :: this
    real(r8kind), INTENT(IN) :: y(:, :)
    real(r8kind) :: yb(:, :)
    real(r8kind) :: dydt(:, :)
    real(r8kind) :: dydtb(:, :)

    ! qprime
    real(r8kind) :: local_qprime(this%nsh2, this%nvl)
    real(r8kind) :: local_qprimeb(this%nsh2, this%nvl)
    ! psi
    real(r8kind) :: local_psi(this%nsh2, this%nvl)
    real(r8kind) :: local_psib(this%nsh2, this%nvl)
    ! psit
    real(r8kind) :: local_psit(this%nsh2, this%ntl)
    real(r8kind) :: local_psitb(this%nsh2, this%ntl)
    ! time derivative of qprime
    real(r8kind) :: dqprdt(this%nsh2, this%nvl)
    real(r8kind) :: dqprdtb(this%nsh2, this%nvl)

    local_psi = 0.0_r8kind
    local_psit = 0.0_r8kind
    local_qprime = this%fstofm(y, this%nm)

    call this%qtopsi(local_qprime, local_psi, local_psit)

    dqprdt = this%ddt(local_psi, local_psit, local_qprime, this%for)

    call this%fmtofs_b(dqprdt, dqprdtb, dydtb)

    call this%ddt_b(local_psi, local_psib, local_psit, local_psitb, local_qprime, local_qprimeb, this%for, dqprdtb)

    local_psi = 0.0_r8kind
    local_psit = 0.0_r8kind
    call this%QTOPSI_B(local_qprime, local_qprimeb, local_psi, local_psib, local_psit, local_psitb)

    call this%fstofm_B(y, yb, this%nm, local_qprimeb)

  end subroutine DQDT_B

  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  function ddt(this, psi, psit, qprime, for) result(dqprdt)

    class(qg_adj_type), intent(in) :: this
    real(r8kind),         intent(in) :: psi(this%nsh2,this%nvl)    ! stream function at the nvl levels
    real(r8kind),         intent(in) :: psit(this%nsh2,this%ntl)   ! thickness at the ntl levels
    real(r8kind),         intent(in) :: qprime(this%nsh2,this%nvl) ! potential vorticity
    real(r8kind),         intent(in) :: for(this%nsh2,this%nvl)    ! constant potential vorticity forcing at the nvl levels
    real(r8kind)                     :: dqprdt(this%nsh2,this%nvl)

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
  !  Differentiation of ddt in reverse (adjoint) mode:
  !   gradient     of useful results: tmp dqprdt
  !   with respect to varying inputs: tmp psi qprime psit
  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  subroutine DDT_B(this, psi, psib, psit, psitb, qprime, qprimeb, for, dqprdtb)

    class(qg_adj_type), intent(inout) :: this
    ! stream function at the nvl levels
    real(r8kind), INTENT(IN) :: psi(this%nsh2, this%nvl)
    real(r8kind) :: psib(this%nsh2, this%nvl)
    ! thickness at the ntl levels
    real(r8kind), INTENT(IN) :: psit(this%nsh2, this%ntl)
    real(r8kind) :: psitb(this%nsh2, this%ntl)
    ! potential vorticity
    real(r8kind), INTENT(IN) :: qprime(this%nsh2, this%nvl)
    real(r8kind) :: qprimeb(this%nsh2, this%nvl)
    ! constant potential vorticity forcing at the nvl levels
    real(r8kind), INTENT(IN) :: for(this%nsh2, this%nvl)
    real(r8kind) :: dqprdt(this%nsh2, this%nvl)
    real(r8kind) :: dqprdtb(this%nsh2, this%nvl)

    integer :: k, l, i, j
    real(r8kind) :: dum1, dum2
    real(r8kind) :: dum1b, dum2b
    real(r8kind), DIMENSION(this%nsh2) :: res
    real(r8kind), DIMENSION(this%nsh2) :: resb
    real(r8kind), DIMENSION(this%nsh2) :: res0
    real(r8kind), DIMENSION(this%nsh2) :: resb0

    ! advection of potential vorticity at upper level
    res = this%JACOB(psi(:, 1), qprime(:, 1))

    ! advection of potential vorticity at middle level
    res0 = this%JACOB(psi(:, 2), qprime(:, 2))

    ! advection of potential vorticity and dissipation at lower level
    qprimeb = 0.0_8
    do l=3,1,-1
      do k=this%nsh2,1,-1
        qprimeb(k, l) = qprimeb(k, l) + this%diss(k, 1) * dqprdtb(k, l)
      end do
    end do
    psitb = 0.0_8
    do k=this%nsh2,1,-1
      dum2b = dqprdtb(k, 2) - dqprdtb(k, 3)
      dum1b = dqprdtb(k, 1) - dqprdtb(k, 2)
      psitb(k, 2) = psitb(k, 2) + this%relt2 * dum2b
      psitb(k, 1) = psitb(k, 1) + this%relt1 * dum1b
    end do
    call this%JACOBD_B(psi(:, 3), psib(:, 3), qprime(:, 3), qprimeb(:, 3), dqprdtb(:, 3))
    dqprdtb(:, 3) = 0.0_8
    resb0 = dqprdtb(:, 2)
    dqprdtb(:, 2) = 0.0_8

    call this%JACOB_B(psi(:, 2), psib(:, 2), qprime(:, 2), qprimeb(:, 2), resb0)
    resb = dqprdtb(:, 1)

    call this%JACOB_B(psi(:, 1), psib(:, 1), qprime(:, 1), qprimeb(:, 1), resb)
  end subroutine DDT_B


  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacob (this, psiloc, pvor) result(sjacob)

    implicit none

    class(qg_adj_type), intent(in) :: this
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
  !  Differentiation of jacob in reverse (adjoint) mode:
  !   gradient     of useful results: tmp sjacob psiloc pvor
  !   with respect to varying inputs: tmp psiloc pvor
  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  subroutine jacob_b(this, psiloc, psilocb, pvor, pvorb, sjacobb)

    class(qg_adj_type), intent(in) :: this
    real(r8kind),       intent(in) :: psiloc(this%nsh2)
    real(r8kind)                   :: psilocb(this%nsh2)
    real(r8kind),       intent(in) :: pvor(this%nsh2)
    real(r8kind)                   :: pvorb(this%nsh2)
    real(r8kind)                   :: sjacobb(this%nsh2)

    integer :: i, j, k
    real(r8kind) :: vv(this%nsh2)
    real(r8kind) :: vvb(this%nsh2)
    real(r8kind) :: dpsidl(this%nlat, this%nlon), dpsidm(this%nlat, this%nlon), dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidlb(this%nlat, this%nlon), dpsidmb(this%nlat, this%nlon), dvordlb(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon), gjacob(this%nlat, this%nlon), dpsidls(this%nsh2)
    real(r8kind) :: dvordmb(this%nlat, this%nlon), gjacobb(this%nlat, this%nlon), dpsidlsb(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vv = reshape(ggsp%ddl(pvor), (/this%nsh2/))
    dvordl = ggsp%sptogg_pp(vv)
    dvordm = ggsp%sptogg_pd(pvor)

    ! space derivatives of streamfunction
    dpsidls = reshape(ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidl = ggsp%sptogg_pp(dpsidls)
    dpsidm = ggsp%sptogg_pd(psiloc)
    dpsidlsb = 0.0_8
    do k = this%nsh2, 1, -1
      dpsidlsb(k) = dpsidlsb(k) - sjacobb(k)
    end do

    gjacobb = ggsp%sptogg_pw(sjacobb)

    dpsidlb = 0.0_8
    dpsidmb = 0.0_8
    dvordlb = 0.0_8
    dvordmb = 0.0_8
    do j = this%nlon, 1, -1
      do i = this%nlat, 1, -1
        dpsidmb(i, j) = dpsidmb(i, j) + dvordl(i, j) * gjacobb(i, j)
        dvordlb(i, j) = dvordlb(i, j) + dpsidm(i, j) * gjacobb(i, j)
        dpsidlb(i, j) = dpsidlb(i, j) - dvordm(i, j) * gjacobb(i, j)
        dvordmb(i, j) = dvordmb(i, j) - dpsidl(i, j) * gjacobb(i, j)
        gjacobb(i, j) = 0.0_8
      end do
    end do

    psilocb = reshape(ggsp%ggtosp(dpsidmb), (/this%nsh2/))
    dpsidlsb = reshape(ggsp%ggtosp(dpsidlb), (/this%nsh2/))

    call ggsp%ddl_b(psiloc, psilocb, dpsidlsb)

    pvorb = reshape(ggsp%ggtosp(dvordmb), (/this%nsh2/))

    vvb = 0.0_8
    vvb = reshape(ggsp%ggtosp(dvordlb), (/this%nsh2/))

    call ggsp%ddl_b(pvor, pvorb, vvb)

  end subroutine jacob_b


  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacobd (this, psiloc, pvor) result(sjacob)

    class(qg_adj_type), intent( in) :: this
    real(r8kind),       intent( in) :: psiloc(this%nsh2)
    real(r8kind),       intent( in) :: pvor(this%nsh2)
    real(r8kind)                    :: sjacob(this%nsh2)

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
  !  Differentiation of jacobd in reverse (adjoint) mode:
  !   gradient     of useful results: tmp sjacob pvor
  !   with respect to varying inputs: tmp psiloc pvor
  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  subroutine jacobd_b(this, psiloc, psilocb, pvor, pvorb, sjacobb)

    class(qg_adj_type), intent(inout) :: this
    real(r8kind),       intent( in) :: psiloc(this%nsh2)
    real(r8kind),       intent(out) :: psilocb(this%nsh2)
    real(r8kind),       intent( in) :: pvor(this%nsh2)
    real(r8kind),       intent(out) :: pvorb(this%nsh2)
    real(r8kind),       intent( in) :: sjacobb(this%nsh2)

    integer :: i, j, k
    real(r8kind) :: dpsidl(this%nlat, this%nlon), dpsidm(this%nlat, this%nlon), dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidlb(this%nlat, this%nlon), dpsidmb(this%nlat, this%nlon), dvordlb(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon), gjacob(this%nlat, this%nlon), vv(this%nsh2)
    real(r8kind) :: dvordmb(this%nlat, this%nlon), gjacobb(this%nlat, this%nlon), vvb(this%nsh2)
    real(r8kind) :: azeta(this%nlat, this%nlon), dpsidls(this%nsh2)
    real(r8kind) :: azetab(this%nlat, this%nlon), dpsidlsb(this%nsh2)
    real(r8kind) :: tempb, tempb0

    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity 
    vv = reshape(ggsp%ddl(pvor), (/this%nsh2/))
    dvordl = ggsp%sptogg_pp(vv)
    dvordm = ggsp%sptogg_pd(pvor)

    ! space derivatives of streamfunction
    dpsidls = reshape(ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidl = ggsp%sptogg_pp(dpsidls)
    dpsidm = ggsp%sptogg_pd(psiloc)

    ! dissipation 
    dpsidlsb = 0.0_r8kind
    do k = this%nsh2, 1, -1
      dpsidlsb(k) = dpsidlsb(k) - sjacobb(k)
    end do

    if (this%lgdiss) then

      gjacobb = ggsp%sptogg_pw(sjacobb)
      dpsidlb = 0.0_r8kind
      dpsidmb = 0.0_r8kind
      azetab = 0.0_r8kind
      do j = this%nlon, 1, -1
        do i = this%nlat, 1, -1
          dpsidmb(i, j) = dpsidmb(i, j) - this%ddisdy(i, j)*gjacobb(i, j)
          azetab(i, j) = azetab(i, j) + this%rdiss(i, j)*gjacobb(i, j)
          dpsidlb(i, j) = dpsidlb(i, j) - this%ddisdx(i, j)*gjacobb(i, j)
        end do
      end do
      vvb = 0.0_r8kind
      vvb = reshape(ggsp%ggtosp(azetab), (/this%nsh2/))
      psilocb = 0.0_r8kind
      do k = this%nsh2, 1, -1
        psilocb(k) = psilocb(k) + this%diss(k, 2) * vvb(k)
        vvb(k) = 0.0_r8kind
      end do

    else

      psilocb = 0.0_r8kind
      do k = this%nsh2, 1, -1
        psilocb(k) = psilocb(k) + this%diss(k, 2) * sjacobb(k)
      end do
      gjacobb = ggsp%sptogg_pw(sjacobb)
      dpsidlb = 0.0_r8kind
      dpsidmb = 0.0_r8kind
      vvb = 0.0_r8kind

    end if

    dvordlb = 0.0_r8kind
    dvordmb = 0.0_r8kind
    do j = this%nlon, 1, -1
      do i = this%nlat, 1, -1
        dpsidmb(i, j) = dpsidmb(i, j) + (this%sinfi(i) * this%dorodl(i, j) + dvordl(i, j)) * gjacobb(i, j)
        dvordlb(i, j) = dvordlb(i, j) + dpsidm(i, j) * gjacobb(i, j)
        dpsidlb(i, j) = dpsidlb(i, j) - (this%sinfi(i) * this%dorodm(i, j) + dvordm(i, j)) * gjacobb(i, j)
        dvordmb(i, j) = dvordmb(i, j) - dpsidl(i, j) * gjacobb(i, j)
        gjacobb(i, j) = 0.0_8
      end do
    end do

    psilocb = reshape(ggsp%ggtosp(dpsidmb), (/this%nsh2/))
    dpsidlsb = reshape(ggsp%ggtosp(dpsidlb), (/this%nsh2/))

    call ggsp%ddl_b(psiloc, psilocb, dpsidlsb)

    pvorb = reshape(ggsp%ggtosp(dvordmb), (/this%nsh2/))
    vvb = reshape(ggsp%ggtosp(dvordlb), (/this%nsh2/))

    call ggsp%ddl_b(pvor, pvorb, vvb)

  end subroutine jacobd_b


  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  subroutine qtopsi(this, qprime, psi, psit)

    class(qg_adj_type), intent( in) :: this
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
  !  Differentiation of qtopsi in reverse (adjoint) mode:
  !   gradient     of useful results: psi psit qprime
  !   with respect to varying inputs: qprime
  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  subroutine QTOPSI_b(this, qprime, qprimeb, psi, psib, psit, psitb)

    ! potential vorticity
    class(qg_adj_type), intent( in) :: this
    real(r8kind), intent(IN) :: qprime(:, :)
    real(r8kind) :: qprimeb(:, :)
    ! stream function at the nvl levels
    real(r8kind) :: psi(:, :)
    real(r8kind) :: psib(:, :)
    ! thickness at the ntl levels
    real(r8kind) :: psit(:, :)
    real(r8kind) :: psitb(:, :)

    integer :: k
    real(r8kind) :: r3
    ! only used as portable workspace
    real(r8kind) :: ws(this%nsh2)
    real(r8kind) :: wsb(this%nsh2)
    INTRINSIC SIZE
    real(r8kind) :: tempb
    real(r8kind) :: tempb0
    integer :: ad_to
    integer :: ad_to0
    integer :: ad_to1

    do k=1,SIZE(psi, 1)
      ws(k) = qprime(k, 1) + qprime(k, 3)
      psi(k, 1) = this%rinhel(k, 1)*(ws(k)+qprime(k, 2))
      psi(k, 2) = ws(k) - 2.d0*qprime(k, 2)
      psi(k, 3) = qprime(k, 1) - qprime(k, 3)
    end do
    do k=1,SIZE(psit, 1)
      psit(k, 1) = this%rinhel(k, 2)*psi(k, 2) + this%rinhel(k, 3)*psi(k, 3)
      psit(k, 2) = this%rinhel(k, 4)*psi(k, 2) + this%rinhel(k, 5)*psi(k, 3)
    end do
    r3 = 1./3.
    do k=1,SIZE(psi, 1)
      psi(k, 2) = r3*(psi(k, 1)-psit(k, 1)+psit(k, 2))
      psi(k, 1) = psi(k, 2) + psit(k, 1)
      psi(k, 3) = psi(k, 2) - psit(k, 2)
    end do
    do k=SIZE(psi, 1),1,-1
      psib(k, 2) = psib(k, 2) + psib(k, 3)
      psitb(k, 2) = psitb(k, 2) - psib(k, 3)
      psib(k, 3) = 0.0_r8kind
      psib(k, 2) = psib(k, 2) + psib(k, 1)
      psitb(k, 1) = psitb(k, 1) + psib(k, 1)
      psib(k, 1) = 0.0_r8kind
      tempb0 = r3*psib(k, 2)
      psib(k, 1) = psib(k, 1) + tempb0
      psitb(k, 1) = psitb(k, 1) - tempb0
      psitb(k, 2) = psitb(k, 2) + tempb0
      psib(k, 2) = 0.0_r8kind
    end do
    do k=SIZE(psit, 1),1,-1
      psib(k, 2) = psib(k, 2) + this%rinhel(k, 4)*psitb(k, 2)
      psib(k, 3) = psib(k, 3) + this%rinhel(k, 5)*psitb(k, 2)
      psitb(k, 2) = 0.0_r8kind
      psib(k, 2) = psib(k, 2) + this%rinhel(k, 2)*psitb(k, 1)
      psib(k, 3) = psib(k, 3) + this%rinhel(k, 3)*psitb(k, 1)
      psitb(k, 1) = 0.0_r8kind
    end do
    wsb = 0.0_r8kind
    do k=SIZE(psi, 1),1,-1
      qprimeb(k, 1) = qprimeb(k, 1) + psib(k, 3)
      qprimeb(k, 3) = qprimeb(k, 3) - psib(k, 3)
      psib(k, 3) = 0.0_r8kind
      wsb(k) = wsb(k) + psib(k, 2)
      qprimeb(k, 2) = qprimeb(k, 2) - 2.d0*psib(k, 2)
      psib(k, 2) = 0.0_r8kind
      tempb = this%rinhel(k, 1)*psib(k, 1)
      wsb(k) = wsb(k) + tempb
      qprimeb(k, 2) = qprimeb(k, 2) + tempb
      psib(k, 1) = 0.0_r8kind
      qprimeb(k, 1) = qprimeb(k, 1) + wsb(k)
      qprimeb(k, 3) = qprimeb(k, 3) + wsb(k)
      wsb(k) = 0.0_r8kind
    end do

  end subroutine QTOPSI_b


  !-----------------------------------------------------------------------
  ! computation of potential vorticity from stream function
  ! input psi streamfunction
  ! output qprime,  the potential vorticity and psit,  the layer thick.
  !-----------------------------------------------------------------------
  subroutine psitoq(this, psi, psit, qprime)
      
    class(qg_adj_type), intent( in) :: this
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

    class(qg_adj_type), intent(in) :: this
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
  !  Differentiation of fmtofs in reverse (adjoint) mode:
  !   gradient     of useful results: z
  !   with respect to varying inputs: y
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
  subroutine fmtofs_b(this, y, yb, zb)

    class(qg_adj_type) :: this
    real(r8kind), intent(IN) :: y(:, :)
    real(r8kind) :: yb(:, :)
    real(r8kind), DIMENSION(SIZE(y, 1), SIZE(y, 2)) :: z
    real(r8kind), DIMENSION(SIZE(y, 1), SIZE(y, 2)) :: zb

    integer :: m, n, k, indx, l
    integer, DIMENSION(SIZE(y, 2)) :: k2

    integer :: max1
    integer :: branch
    integer :: ad_from
    integer :: ad_to

    k = 0
    do l=1,SIZE(y, 2)
      k2(l) = k
      k = 1
      do m=0,this%nm
        IF (m .LT. 1) THEN
          max1 = 1
        ELSE
          max1 = m
        end IF
        ad_from = max1
        do n=ad_from, this%nm
          k = k + 1
          IF (m .EQ. 0) THEN
            indx = n**2
          ELSE
            indx = n**2 + 2*m - 1
          end IF
        end do
      end do
    end do
    yb = 0.0_r8kind
    do l=SIZE(y, 2),1,-1
      do m= this%nm,0,-1
        IF (m .LT. 1) THEN
          max1 = 1
        ELSE
          max1 = m
        end IF
        ad_from = max1
        do n= this%nm,ad_from,-1
          IF (m .EQ. 0) THEN
            indx = n**2
          ELSE
            indx = n**2 + 2*m - 1
          end IF
          IF (m .NE. 0) THEN
            yb(k+this%nsh, l) = yb(k+this%nsh, l) + zb(indx+1, l)
            zb(indx+1, l) = 0.0_r8kind
          end IF
          yb(k, l) = yb(k, l) + zb(indx, l)
          zb(indx, l) = 0.0_r8kind
          k = k - 1
        end do
      end do
      k = k2(l)
    end do

  end subroutine fmtofs_b


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

    class(qg_adj_type), intent(in) :: this
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
  !  Differentiation of fstofm in reverse (adjoint) mode:
  !   gradient     of useful results: y z
  !   with respect to varying inputs: y
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
  subroutine fstofm_b(this, y, yb, ntr, zb)

    class(qg_adj_type), intent(in) :: this
    real(r8kind), intent(IN) :: y(:, :)
    real(r8kind) :: yb(:, :)
    integer, intent(IN) :: ntr
    real(r8kind), DIMENSION(SIZE(y, 1), SIZE(y, 2)) :: z
    real(r8kind), DIMENSION(SIZE(y, 1), SIZE(y, 2)) :: zb

    integer :: m, n, k, indx, i, l
    integer, DIMENSION(SIZE(y, 2)) :: k2

    integer :: max1
    integer :: ad_to
    integer :: branch
    integer :: ad_from
    integer :: ad_to0

    k = 0
    do l=1,SIZE(y, 2)

      k2(l) = k
      k = 1
      do m=0,this%nm
        IF (m .LT. 1) THEN
          max1 = 1
        ELSE
          max1 = m
        end IF
        ad_from = max1
        do n=ad_from,this%nm
          k = k + 1
          IF (m .LE. ntr .AND. n .LE. ntr) THEN
            IF (m .EQ. 0) THEN
              indx = n**2
            ELSE
              indx = n**2 + 2*m - 1
            end IF
          end IF
        end do
      end do
    end do
    ad_to0 = l - 1
    do l=ad_to0,1,-1
      do m=this%nm,0,-1
        IF (m .LT. 1) THEN
          max1 = 1
        ELSE
          max1 = m
        end IF
        ad_from = max1
        do n=this%nm,ad_from,-1
          IF (m .LE. ntr .AND. n .LE. ntr) THEN
            IF (m .EQ. 0) THEN
              indx = n**2
            ELSE
              indx = n**2 + 2*m - 1
            end IF
            IF (m .NE. 0) THEN
              yb(indx+1, l) = yb(indx+1, l) + zb(k+this%nsh, l)
              zb(k+this%nsh, l) = 0.0_r8kind
            end IF
            yb(indx, l) = yb(indx, l) + zb(k, l)
            zb(k, l) = 0.0_r8kind
          end if
          k = k - 1
        end do
      end do

      k = k2(l)

    end do

  end subroutine fstofm_b


  !-----------------------------------------------------------------------
  ! computation of geostrophic winds at all levels
  ! computes geopotential height in [m2 / s2[ = [gz] from streamfunction 
  ! by solving the linear balance equation: 
  ! del phi = (1 - mu**2 ) d psi / dmu + mu del psi
  ! the global mean value is not determined and set to zero
  !-----------------------------------------------------------------------
  subroutine gridfields(this, lat, lon, lvl, geopg, psig, forg, qgpv, ug, vg)

    class(qg_adj_type), intent( in) :: this
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

    class(qg_adj_type), intent( in) :: this
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

    class(qg_adj_type), intent( in) :: this
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

    class(qg_adj_type),  intent(in) :: this
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

    class(qg_adj_type),   intent(in) :: this
    real(r8kind),           intent(in) :: xsl(:)
    real(r8kind), dimension(size(xsl)) :: xs

    integer :: k

    do k = 1, size(xsl)
      xs(k) = xsl(k) * this%rinhel(k, 1)
    enddo

    return

  end function lapinv


  !-------------------------------------------------------------------------------
  ! get_config
  !-------------------------------------------------------------------------------
  function get_config(this) result(config)

    class(qg_adj_type), intent(in) :: this
    type(qg_config_type)             :: config

    config = this%config

  end function get_config


  !-------------------------------------------------------------------------------
  ! get_step
  !-------------------------------------------------------------------------------
  function get_step(this) result(step)

    class(qg_adj_type), intent(in) :: this
    integer                          :: step

    step = this%step

  end function get_step


  !-------------------------------------------------------------------------------
  ! get_clock
  !-------------------------------------------------------------------------------
  function get_clock(this) result(clock)

    class(qg_adj_type), intent(in) :: this
    real(r8kind)                     :: clock

    clock = this%clock

  end function get_clock


  !-------------------------------------------------------------------------------
  ! get_nlat
  !-------------------------------------------------------------------------------
  function get_nlat(this) result(nlat)

    class(qg_adj_type), intent(in) :: this
    integer                          :: nlat

    nlat = this%nlat

  end function get_nlat


  !-------------------------------------------------------------------------------
  ! get_nlon
  !-------------------------------------------------------------------------------
  function get_nlon(this) result(nlon)

    class(qg_adj_type), intent(in) :: this
    integer                          :: nlon

    nlon = this%nlon

  end function get_nlon


  !-------------------------------------------------------------------------------
  ! get_nsh2
  !-------------------------------------------------------------------------------
  function get_nsh2(this) result(nsh2)

    class(qg_adj_type), intent(in) :: this
    integer                          :: nsh2

    nsh2 = this%nsh2

  end function get_nsh2

  !-------------------------------------------------------------------------------
  ! get_nvl
  !-------------------------------------------------------------------------------
  function get_nvl(this) result(nvl)

    class(qg_adj_type), intent(in) :: this
    integer                          :: nvl

    nvl = this%nvl

  end function get_nvl


  !-------------------------------------------------------------------------------
  ! get_psi
  !-------------------------------------------------------------------------------
  function get_psi(this) result(psi)

    class(qg_adj_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: psi

    psi = this%psi

  end function get_psi


  !-------------------------------------------------------------------------------
  ! get_trajectory
  !-------------------------------------------------------------------------------
  function get_trajectory(this) result(trajectory)

    class(qg_adj_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: trajectory

    trajectory = this%trajectory

  end function get_trajectory


  !-------------------------------------------------------------------------------
  ! get_for
  !-------------------------------------------------------------------------------
  function get_for(this) result(for)

    class(qg_adj_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: for

    for = this%for

  end function get_for


  !-------------------------------------------------------------------------------
  ! get_state_vector
  !-------------------------------------------------------------------------------
  function get_state_vector(this) result(state_vector)

    class(qg_adj_type),            intent(in) :: this
    real(r8kind), dimension(this%nlat * this%nlon * this%nvl) :: state_vector

    state_vector = reshape(this%get_psig(),(/this%nlat * this%nlon * this%nvl/))

  end function get_state_vector


  !-------------------------------------------------------------------------------
  ! get_location_vector
  !-------------------------------------------------------------------------------
  function get_location_vector(this) result(location_vector)

    class(qg_adj_type),            intent(in) :: this
    real(r8kind), dimension(this%nlat * this%nlon * this%nvl, 3) :: location_vector

    location_vector = reshape(this%get_lat_lon_grid(),(/this%nlat * this%nlon * this%nvl, 3/))

  end function get_location_vector


  !------------------------------------------------------------------
  ! get_interpolation_weights
  !------------------------------------------------------------------
  subroutine get_interpolation_weights(this, latx, lonx, lvl, NW_index, NE_index, SW_index, SE_index, NW_weight, NE_weight, SW_weight, SE_weight)

    class(qg_adj_type), intent( in) :: this
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

end module QG_Model_ADJ
