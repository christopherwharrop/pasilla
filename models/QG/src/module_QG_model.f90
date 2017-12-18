module QG_Model

  use kind
  use Abstract_Model, only : abstract_model_type
  use QG_Config
  use QG_GGSP

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! qg_model_type
  !-------------------------------------------------------------------------------
  public :: qg_model_type

  type, extends(abstract_model_type) :: qg_model_type
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

    ! Current model step and simulaiton clock
    integer      :: step
    real(r8kind) :: clock

    ! Model time step
    real(r8kind) :: dtt                      ! dimensionless time step

    ! Model state
    real(r8kind), allocatable :: psi(:,:)    ! Stream function at the nvl levels
    real(r8kind), allocatable :: psit(:,:)   ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprime(:,:) ! Potential vorticity

    ! Model Forcing
    real(r8kind), allocatable :: for(:,:)    ! Constant potential vorticity forcing at the nvl levels

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

  contains
    final :: destructor_qg_model
    procedure :: adv_nsteps => adv_nsteps_model
    procedure :: gridfields
    procedure :: get_config
    procedure :: get_step
    procedure :: get_clock
    procedure :: get_nlat
    procedure :: get_nlon
    procedure :: get_nsh2
    procedure :: get_nvl
    procedure :: get_psi
    procedure :: get_psig
    procedure :: get_lat_lon_grid
    procedure :: get_for
    procedure :: get_state_vector
    procedure :: get_location_vector
    procedure :: get_interpolation_weights
    procedure :: ggvtoss
    procedure :: sstoggv
    procedure :: rebalance
    procedure, private :: init_state
    procedure, private :: dqdt
    procedure, private :: ddt
    procedure, private :: jacob
    procedure, private :: jacobd
    procedure, private :: psitoq
    procedure, private :: qtopsi
    procedure, private :: lap
    procedure, private :: lapinv
    procedure, private :: fmtofs
    procedure, private :: fstofm
  end type qg_model_type

  interface qg_model_type
    procedure :: constructor_qg_model
  end interface

  !-------------------------------------------------------------------------------
  ! qg_tl_type
  !-------------------------------------------------------------------------------
  public :: qg_tl_type

  type, extends(qg_model_type) :: qg_tl_type
    private

    ! Model state
    real(r8kind), allocatable :: trajectory(:,:) ! Tangent Linear trajectory
    real(r8kind), allocatable :: psitd(:,:)       ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprimed(:,:)     ! Potential vorticity
    
  contains
    final :: destructor_qg_tl
    procedure :: adv_nsteps => adv_nsteps_tl
    procedure :: get_trajectory => get_trajectory_tl
    procedure :: get_qprimed
    procedure, private :: dqdt_d
    procedure, private :: ddt_d
    procedure, private :: jacob_d
    procedure, private :: jacobd_d
  end type qg_tl_type

  interface qg_tl_type
    procedure :: constructor_qg_tl
  end interface

  interface
    module function constructor_qg_tl(config, state, state_vector, trajectory, trajectory_vector, for, step) result (qg_tl)
      type(qg_config_type),   intent(in) :: config
      real(r8kind), optional, intent(in) :: state(:,:)
      real(r8kind), optional, intent(in) :: state_vector(:)
      real(r8kind), optional, intent(in) :: trajectory(:,:)
      real(r8kind), optional, intent(in) :: trajectory_vector(:)
      real(r8kind), optional, intent(in) :: for(:,:)
      integer,      optional, intent(in) :: step
      type(qg_tl_type)                   :: qg_tl
    end function constructor_qg_tl
    module elemental subroutine destructor_qg_tl(this)
      type(qg_tl_type), intent(inout) :: this
    end subroutine destructor_qg_tl
    module subroutine adv_nsteps_tl(this, nsteps)
      class(qg_tl_type), intent(inout) :: this
      integer          , intent(   in) :: nsteps
    end subroutine adv_nsteps_tl
    module subroutine dqdt_d(this, y, yd, dydt, dydtd)
      class(qg_tl_type) :: this
      real(r8kind), intent( in) :: y(:,:)
      real(r8kind), intent( in) :: yd(:,:)
      real(r8kind), intent(out) :: dydt(:,:)
      real(r8kind), intent(out) :: dydtd(:,:)
    end subroutine dqdt_d
    module subroutine ddt_d(this, psi, psid, psit, psitd, qprime, qprimed, for, dqprdt, dqprdtd)
      class(qg_tl_type), intent( in) :: this
      real(r8kind),      intent( in) :: psi(:,:)    ! stream function at the nvl levels
      real(r8kind),      intent( in) :: psid(:,:)    ! stream function at the nvl levels
      real(r8kind),      intent( in) :: psit(:,:)   ! thickness at the ntl levels
      real(r8kind),      intent( in) :: psitd(:,:)   ! thickness at the ntl levels
      real(r8kind),      intent( in) :: qprime(:,:) ! potential vorticity
      real(r8kind),      intent( in) :: qprimed(:,:) ! potential vorticity
      real(r8kind),      intent( in) :: for(:,:)    ! constant potential vorticity forcing at the nvl levels
      real(r8kind),      intent(out) :: dqprdt(:,:)
      real(r8kind),      intent(out) :: dqprdtd(:,:)
    end subroutine ddt_d
    module subroutine jacob_d(this, psiloc, psilocd, pvor, pvord, sjacob, sjacobd)
      class(qg_tl_type), intent(in) :: this
      real(r8kind), intent( in) :: psiloc(:)
      real(r8kind), intent( in) :: psilocd(:)
      real(r8kind), intent( in) :: pvor(:)
      real(r8kind), intent( in) :: pvord(:)
      real(r8kind), intent(out) :: sjacob(:)
      real(r8kind), intent(out) :: sjacobd(:)
    end subroutine jacob_d
    module subroutine jacobd_d(this, psiloc, psilocd, pvor, pvord, sjacob, sjacobd)
      class(qg_tl_type), intent( in) :: this
      real(r8kind),      intent( in) :: psiloc(:)
      real(r8kind),      intent( in) :: psilocd(:)
      real(r8kind),      intent( in) :: pvor(:)
      real(r8kind),      intent( in) :: pvord(:)
      real(r8kind),      intent(out) :: sjacob(:)
      real(r8kind),      intent(out) :: sjacobd(:)
    end subroutine jacobd_d
    module function get_trajectory_tl(this) result(trajectory)
      class(qg_tl_type),            intent(in) :: this
      real(r8kind), allocatable :: trajectory(:,:)
    end function get_trajectory_tl
    module function get_qprimed(this) result(qprimed)
      class(qg_tl_type),            intent(in) :: this
      real(r8kind), allocatable :: qprimed(:,:)
    end function get_qprimed
  end interface

  !-------------------------------------------------------------------------------
  ! qg_adj_type
  !-------------------------------------------------------------------------------
  public :: qg_adj_type

  type, extends(qg_model_type) :: qg_adj_type
    private

    ! Model state
    real(r8kind), allocatable :: trajectory(:,:) ! Adjoint trajectory
    real(r8kind), allocatable :: psitb(:,:)      ! Thickness at the ntl levels
    real(r8kind), allocatable :: qprimeb(:,:)    ! Potential vorticity

  contains
    final :: destructor_qg_adj
    procedure :: adv_nsteps => adv_nsteps_adj
    procedure :: get_trajectory => get_trajectory_adj
    procedure, private :: dqdt_b
    procedure, private :: ddt_b
    procedure, private :: jacob_b
    procedure, private :: jacobd_b
    procedure, private :: qtopsi_b
    procedure, private :: fmtofs_b
    procedure, private :: fstofm_b
  end type qg_adj_type

  interface qg_adj_type
    procedure :: constructor_qg_adj
  end interface

  interface
    module function constructor_qg_adj(config, state, state_vector, trajectory, trajectory_vector, for, step) result (qg_adj)
      type(qg_config_type),   intent(in) :: config
      real(r8kind), optional, intent(in) :: state(:,:)
      real(r8kind), optional, intent(in) :: state_vector(:)
      real(r8kind), optional, intent(in) :: trajectory(:,:)
      real(r8kind), optional, intent(in) :: trajectory_vector(:)
      real(r8kind), optional, intent(in) :: for(:,:)
      integer,      optional, intent(in) :: step
      type(qg_adj_type)                  :: qg_adj
    end function constructor_qg_adj
    module elemental subroutine destructor_qg_adj(this)
      type(qg_adj_type), intent(inout) :: this
    end subroutine destructor_qg_adj
    module subroutine adv_nsteps_adj(this, nsteps)
      class(qg_adj_type), intent(inout) :: this
      integer           , intent(   in) :: nsteps
    end subroutine adv_nsteps_adj
    module subroutine dqdt_b(this, y, yb, dydt, dydtb)
      class(qg_adj_type), intent(inout) :: this
      real(r8kind),       intent(   in) :: y(:,:)
      real(r8kind)                      :: yb(:,:)
      real(r8kind)                      :: dydt(:,:)
      real(r8kind)                      :: dydtb(:,:)
    end subroutine dqdt_b
    module subroutine ddt_b(this, psi, psib, psit, psitb, qprime, qprimeb, for, dqprdtb)
      class(qg_adj_type), intent(inout) :: this    
      real(r8kind),       intent(   in) :: psi(:,:)
      real(r8kind)                      :: psib(:,:)
      real(r8kind),       intent(   in) :: psit(:,:)
      real(r8kind)                      :: psitb(:,:)
      real(r8kind),       intent(   in) :: qprime(:,:)
      real(r8kind)                      :: qprimeb(:,:)
      real(r8kind),       intent(   in) :: for(:,:)
      real(r8kind)                      :: dqprdtb(:,:)
    end subroutine ddt_b
    module subroutine jacob_b(this, psiloc, psilocb, pvor, pvorb, sjacobb)
      class(qg_adj_type), intent(in) :: this
      real(r8kind),       intent(in) :: psiloc(:)
      real(r8kind)                   :: psilocb(:)
      real(r8kind),       intent(in) :: pvor(:)
      real(r8kind)                   :: pvorb(:)
      real(r8kind)                   :: sjacobb(:)
    end subroutine jacob_b
    module subroutine jacobd_b(this, psiloc, psilocb, pvor, pvorb, sjacobb)
      class(qg_adj_type), intent(inout) :: this
      real(r8kind),       intent(   in) :: psiloc(:)
      real(r8kind),       intent(  out) :: psilocb(:)
      real(r8kind),       intent(   in) :: pvor(:)
      real(r8kind),       intent(  out) :: pvorb(:)
      real(r8kind),       intent(   in) :: sjacobb(:)
    end subroutine jacobd_b
    module subroutine qtopsi_b(this, qprime, qprimeb, psi, psib, psit, psitb)
      class(qg_adj_type), intent(in) :: this
      real(r8kind),       intent(in) :: qprime(:,:)   ! potential vorticity
      real(r8kind)                   :: qprimeb(:,:)
      real(r8kind)                   :: psi(:,:)      ! stream function at the nvl levels
      real(r8kind)                   :: psib(:,:)
      real(r8kind)                   :: psit(:,:)     ! thickness at the ntl levels
      real(r8kind)                   :: psitb(:,:)
    end subroutine qtopsi_b
    module subroutine fmtofs_b(this, y, yb, zb)
      class(qg_adj_type)       :: this
      real(r8kind), intent(in) :: y(:, :)
      real(r8kind)             :: yb(:, :)
      real(r8kind)             :: zb(:,:)
    end subroutine fmtofs_b
    module subroutine fstofm_b(this, y, yb, ntr, zb)
      class(qg_adj_type), intent(in) :: this
      real(r8kind),       intent(in) :: y(:,:)
      real(r8kind)                   :: yb(:,:)
      integer,            intent(in) :: ntr
      real(r8kind)                   :: zb(:,:)
    end subroutine fstofm_b
    module function get_trajectory_adj(this) result(trajectory)
      class(qg_adj_type), intent(in) :: this
      real(r8kind), allocatable      :: trajectory(:,:)
    end function get_trajectory_adj
  end interface

  !-------------------------------------------------------------------------------
  ! Constants
  !-------------------------------------------------------------------------------

  ! Mathematical and physical constants
  real(r8kind), parameter :: pi = 4d0 * atan(1d0)     ! value of pi


contains


  !-------------------------------------------------------------------------------
  ! constructor_qg_model
  !-------------------------------------------------------------------------------
  function constructor_qg_model(config, state, state_vector, step) result (qg_model)

    type(qg_config_type),   intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:,:)
    real(r8kind), optional, intent(in) :: state_vector(:)
    integer,      optional, intent(in) :: step

    type(qg_model_type)                :: qg_model
    real(r8kind), allocatable          :: psig3d(:,:,:), psisp(:,:)

    real(r8kind) :: nsteps_per_day

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    ! Set the model config
    qg_model%config = config

    ! Get model grid dimensions
    qg_model%ft = config%get_ft()
    qg_model%nm = config%get_nm()
    qg_model%nlon = config%get_nlon()
    qg_model%nlat = config%get_nlat()
    qg_model%nvl = config%get_nvl()
    qg_model%ntl = config%get_ntl()
    qg_model%nsh = config%get_nsh()
    qg_model%nsh2 = config%get_nsh2()

    ! Get model grid converter
    qg_model%ggsp = config%get_ggsp()

    ! Get model forcing from configuration
    allocate(qg_model%for, source = config%get_for())

    ! Get Nondimensional relaxation coefficients
    qg_model%relt1 = config%get_relt1()
    qg_model%relt2 = config%get_relt2()

    ! Get dissipation coefficients for each spherical harmonic
    allocate(qg_model%diss, source = config%get_diss())

    ! Get gauss points in radians and sine and cosine of it
    allocate(qg_model%phi, source = config%get_phi())
    allocate(qg_model%sinfi, source = config%get_sinfi())
    allocate(qg_model%cosfi, source = config%get_cosfi())

    ! Get derivatives of orog
    allocate(qg_model%dorodl, source = config%get_dorodl())
    allocate(qg_model%dorodm, source = config%get_dorodm())

    ! Get orography and land-sea mask dependent friction option
    qg_model%lgdiss = config%get_lgdiss()

    ! Get landsea-mask/orography dependent friction
    allocate(qg_model%rdiss, source = config%get_rdiss())
    allocate(qg_model%ddisdx, source = config%get_ddisdx())
    allocate(qg_model%ddisdy, source = config%get_ddisdy())

    ! Get Laplace and Helmholtz operator for Q-PSI inversion
    allocate(qg_model%rinhel(qg_model%nsh2,0:5))
    qg_model%rinhel = config%get_rinhel()

    ! Get one over Rossby rad. of def. squared of 200-500, and 500-800 thicknesses
    qg_model%rl1 = config%get_rl1()
    qg_model%rl2 = config%get_rl2()

    ! Initialize model step and clock
    if (present(step)) then
      qg_model%step = step
    else
      qg_model%step = 0
    end if
    qg_model%clock = qg_model%step * config%get_time_step()

    ! Initialize time step of the model:
    nsteps_per_day = 24.0d0 * 3600.0d0 / real(config%get_time_step())
    qg_model%dtt = (1d0 / nsteps_per_day) * pi * 4d0

    ! Initialize streamfunction
    if (present(state)) then
      call qg_model%init_state(psi=state)
    else if (present(state_vector)) then
!      call qg_model%init_state(psi=reshape(state_vector,(/qg_model%nsh2, qg_model%nvl/)))
      allocate(psig3d(qg_model%nlon, qg_model%nlat, qg_model%nvl))
      allocate(psisp(qg_model%nsh2, qg_model%nvl))
      psig3d = reshape(state_vector,(/qg_model%nlon, qg_model%nlat, qg_model%nvl/))
      psig3d(:,:,:) = psig3d(:,:,:) / facsf
      psisp(:,1) = reshape(qg_model%ggsp%ggtosp(transpose(psig3d(:,:,1))), (/qg_model%nsh2/))
      psisp(:,2) = reshape(qg_model%ggsp%ggtosp(transpose(psig3d(:,:,2))), (/qg_model%nsh2/))
      psisp(:,3) = reshape(qg_model%ggsp%ggtosp(transpose(psig3d(:,:,3))), (/qg_model%nsh2/))
      call qg_model%init_state(psi=psisp)
    else
      call qg_model%init_state()
    end if

  end function constructor_qg_model


  !------------------------------------------------------------------
  ! destructor_qg_model
  !
  ! Deallocates pointers used by a qg_model_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_qg_model(this)

    type(qg_model_type), intent(inout) :: this

    ! No pointers in qg_model_type object so we do nothing

  end subroutine destructor_qg_model


  !-------------------------------------------------------------------------------
  ! init_state
  !-------------------------------------------------------------------------------
  subroutine init_state(this, psi)

    use netcdf
    use netcdf_utilities

    class(qg_model_type),   intent(inout) :: this
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
  subroutine adv_nsteps_model(this, nsteps)

    class(qg_model_type), intent(inout) :: this
    integer             , intent(   in) :: nsteps

    integer :: step, k, l, nvar
    real(r8kind) :: dt2, dt6
    real(r8kind) :: y(this%nsh2, this%nvl), dydt(this%nsh2, this%nvl), yt(this%nsh2, this%nvl)
    real(r8kind) :: dyt(this%nsh2, this%nvl), dym(this%nsh2, this%nvl)

    if (nsteps > 0) then

      nvar = (this%nm + 2) * this%nm
      dt2 = this%dtt * 0.5d0
      dt6 = this%dtt / 6d0

      ! Advance the model forward in time n steps
      do step = 1, nsteps
        y = this%fmtofs(this%qprime)
        call this%dqdt(y, dydt)
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + dt2 * dydt(k, l)
          enddo
        enddo
        call this%dqdt(yt, dyt)
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + dt2 * dyt(k, l)
          enddo
        enddo
        call this%dqdt(yt, dym)
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + this%dtt * dym(k, l)
            dym(k, l) = dyt(k, l) + dym(k, l)
          enddo
        enddo
        call this%dqdt(yt, dyt)
        do l = 1, this%nvl
          do k = 1, nvar
            y(k, l) = y(k, l) + dt6 * (dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
          enddo
        enddo
        this%qprime = this%fstofm(y, this%nm)

        ! Inrement the step count
        this%step = this%step + 1

      end do

      ! Make stream function consistent with potential vorticity
      call this%qtopsi(this%qprime, this%psi, this%psit)

    end if

  end subroutine adv_nsteps_model


  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  subroutine dqdt(this, y, dydt)

    class(qg_model_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: y(:,:)
    real(r8kind),         intent(  out) :: dydt(:,:)

    real(r8kind) :: qprime(this%nsh2,this%nvl) ! qprime
    real(r8kind) :: psi(this%nsh2,this%nvl)    ! psi
    real(r8kind) :: psit(this%nsh2,this%ntl)   ! psit
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
  function ddt(this, psi, psit, qprime, for) result(dqprdt)

    class(qg_model_type), intent(in) :: this
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
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacob (this, psiloc, pvor) result(sjacob)

    implicit none

    class(qg_model_type), intent(in) :: this
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
  function jacobd (this, psiloc, pvor) result(sjacob)

    class(qg_model_type), intent(in) :: this
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


  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  subroutine qtopsi(this, qprime, psi, psit)

    class(qg_model_type), intent( in) :: this
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
      
    class(qg_model_type), intent( in) :: this
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

    class(qg_model_type), intent(in) :: this
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

    class(qg_model_type), intent(in) :: this
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

    class(qg_model_type), intent( in) :: this
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

    class(qg_model_type), intent( in) :: this
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

    class(qg_model_type), intent( in) :: this
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

    class(qg_model_type),  intent(in) :: this
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

    class(qg_model_type),   intent(in) :: this
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

    class(qg_model_type), intent(in) :: this
    type(qg_config_type)             :: config

    config = this%config

  end function get_config


  !-------------------------------------------------------------------------------
  ! get_step
  !-------------------------------------------------------------------------------
  function get_step(this) result(step)

    class(qg_model_type), intent(in) :: this
    integer                          :: step

    step = this%step

  end function get_step


  !-------------------------------------------------------------------------------
  ! get_clock
  !-------------------------------------------------------------------------------
  function get_clock(this) result(clock)

    class(qg_model_type), intent(in) :: this
    real(r8kind)                     :: clock

    clock = this%clock

  end function get_clock


  !-------------------------------------------------------------------------------
  ! get_nlat
  !-------------------------------------------------------------------------------
  function get_nlat(this) result(nlat)

    class(qg_model_type), intent(in) :: this
    integer                          :: nlat

    nlat = this%nlat

  end function get_nlat


  !-------------------------------------------------------------------------------
  ! get_nlon
  !-------------------------------------------------------------------------------
  function get_nlon(this) result(nlon)

    class(qg_model_type), intent(in) :: this
    integer                          :: nlon

    nlon = this%nlon

  end function get_nlon


  !-------------------------------------------------------------------------------
  ! get_nsh2
  !-------------------------------------------------------------------------------
  function get_nsh2(this) result(nsh2)

    class(qg_model_type), intent(in) :: this
    integer                          :: nsh2

    nsh2 = this%nsh2

  end function get_nsh2

  !-------------------------------------------------------------------------------
  ! get_nvl
  !-------------------------------------------------------------------------------
  function get_nvl(this) result(nvl)

    class(qg_model_type), intent(in) :: this
    integer                          :: nvl

    nvl = this%nvl

  end function get_nvl


  !-------------------------------------------------------------------------------
  ! get_psi
  !-------------------------------------------------------------------------------
  function get_psi(this) result(psi)

    class(qg_model_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: psi

    psi = this%psi

  end function get_psi


  !-------------------------------------------------------------------------------
  ! get_for
  !-------------------------------------------------------------------------------
  function get_for(this) result(for)

    class(qg_model_type),            intent(in) :: this
    real(r8kind), dimension(this%nsh2,this%nvl) :: for

    for = this%for

  end function get_for


  !-------------------------------------------------------------------------------
  ! get_state_vector
  !-------------------------------------------------------------------------------
  function get_state_vector(this) result(state_vector)

    class(qg_model_type),            intent(in) :: this
    real(r8kind), dimension(this%nlat * this%nlon * this%nvl) :: state_vector

    state_vector = reshape(this%get_psig(),(/this%nlat * this%nlon * this%nvl/))

  end function get_state_vector


  !-------------------------------------------------------------------------------
  ! get_location_vector
  !-------------------------------------------------------------------------------
  function get_location_vector(this) result(location_vector)

    class(qg_model_type),            intent(in) :: this
    real(r8kind), dimension(this%nlat * this%nlon * this%nvl, 3) :: location_vector

    location_vector = reshape(this%get_lat_lon_grid(),(/this%nlat * this%nlon * this%nvl, 3/))

  end function get_location_vector


  !------------------------------------------------------------------
  ! get_interpolation_weights
  !------------------------------------------------------------------
  subroutine get_interpolation_weights(this, latx, lonx, lvl, NW_index, NE_index, SW_index, SE_index, NW_weight, NE_weight, SW_weight, SE_weight)

    class(qg_model_type), intent( in) :: this
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


  !------------------------------------------------------------------
  ! distance
  !------------------------------------------------------------------
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


  !------------------------------------------------------------------
  ! ggvtoss
  !------------------------------------------------------------------
  function ggvtoss(this, ggv) result (ss)

    class(qg_model_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: ggv(:)

    real(r8kind), dimension(this%nsh2, this%nvl) :: ss

    real(r8kind), allocatable :: gg3d(:,:,:)

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    allocate(gg3d(this%nlon, this%nlat, this%nvl))
    gg3d = reshape(ggv,(/this%nlon, this%nlat, this%nvl/))
    gg3d(:,:,:) = gg3d(:,:,:) / facsf
    ss(:,1) = reshape(this%ggsp%ggtosp(transpose(gg3d(:,:,1))), (/this%nsh2/))
    ss(:,2) = reshape(this%ggsp%ggtosp(transpose(gg3d(:,:,2))), (/this%nsh2/))
    ss(:,3) = reshape(this%ggsp%ggtosp(transpose(gg3d(:,:,3))), (/this%nsh2/))
print *, "ggvtoss: ss(1,:) = ", ss(1,:)
!    ss(1,:) = 0.0_r8kind

  end function


  !------------------------------------------------------------------
  ! sstoggv
  !------------------------------------------------------------------
  function sstoggv(this, ss) result (ggv)

    class(qg_model_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: ss(:,:)

    real(r8kind), dimension(this%nlon * this%nlat * this%nvl) :: ggv

    real(r8kind), allocatable :: gg3d(:,:,:)

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    allocate(gg3d(this%nlon, this%nlat, this%nvl))
    gg3d(:,:,1) = transpose(this%ggsp%sptogg_pp(ss(:, 1)))
    gg3d(:,:,2) = transpose(this%ggsp%sptogg_pp(ss(:, 2)))
    gg3d(:,:,3) = transpose(this%ggsp%sptogg_pp(ss(:, 3)))
    gg3d(:,:,:) = gg3d(:,:,:) * facsf
    ggv = reshape(gg3d,(/this%nlon * this%nlat * this%nvl/))

  end function


  !------------------------------------------------------------------
  ! rebalance
  !------------------------------------------------------------------
  function rebalance(this, ggv) result (ggvnew)

    class(qg_model_type), intent(inout) :: this
    real(r8kind),         intent(   in) :: ggv(:)

    real(r8kind), dimension(this%nlat * this%nlon * this%nvl) :: ggvnew

    real(r8kind), allocatable :: ss(:,:)
    real(r8kind), allocatable :: gg3d(:,:,:)

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    allocate(ss(this%nsh2, this%nvl))

    ss = this%ggvtoss(ggv)
print *, "rebalance: BEFORE ss(1,:) = ", ss(1,:)
!    ss(:,1) = ss(:,1) - sum(ss(:,1)) / this%nsh2
!    ss(:,2) = ss(:,2) - sum(ss(:,2)) / this%nsh2
!    ss(:,3) = ss(:,3) - sum(ss(:,3)) / this%nsh2
    ss(2:this%nsh2,1) = ss(2:this%nsh2,1) - sum(ss(2:this%nsh2,1)) / this%nsh2 - 1
    ss(2:this%nsh2,2) = ss(2:this%nsh2,2) - sum(ss(2:this%nsh2,2)) / this%nsh2 - 1
    ss(2:this%nsh2,3) = ss(2:this%nsh2,3) - sum(ss(2:this%nsh2,3)) / this%nsh2 - 1
print *, "rebalance: AFTER ss(1,:) = ", ss(1,:)

    ggvnew = this%sstoggv(ss)

  end function

end module QG_Model
