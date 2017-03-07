module QG_Config

  use kind,   only : r8kind
!  use Config, only : config_type

  implicit none

  private

  public :: qg_config_type

  type :: qg_config_type
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

  contains
    final :: destructor_qg_config
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
  end type qg_config_type

  interface qg_config_type
    procedure constructor_arglist
    procedure constructor_namelist_file
    procedure constructor_namelist_unit
  end interface


contains


  !------------------------------------------------------------------
  ! constructor_arglist
  !
  ! Returns an initialized qg_config_type object
  !------------------------------------------------------------------
  function constructor_arglist(resolution, time_step, obsfile, inf, obsf, tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2) result(qg_config)

    integer,           intent(in) :: resolution ! Model resolution
    integer,           intent(in) :: time_step  ! Model time step (in seconds)
    character(len=32), intent(in) :: obsfile    ! Name of observation file
    logical,           intent(in) :: inf        ! If .true. then artificial PV forcing read from file
    logical,           intent(in) :: obsf       ! If .true. PV forcing is calculated from observations in routine artiforc
    real(r8kind),      intent(in) :: tdis       ! Ekman dissipation timescale in days at lower level
    real(r8kind),      intent(in) :: addisl     ! Parameter used in the computation of the dissipation timescale at the lower level over land
    real(r8kind),      intent(in) :: addish     ! Parameter used in the computation of the dissipation timescale at the lower level as a function of topography
    real(r8kind),      intent(in) :: trel       ! Relaxation time scale in days of the temperature
    real(r8kind),      intent(in) :: tdif       ! Dissipation timescale of scale-selective diffusion in days for wavenumber nm
    integer,           intent(in) :: idif       ! Determines scale-selectivity of hyperviscosity; power of laplace operator
    real(r8kind),      intent(in) :: h0         ! Scale factor for the topographically induced upward motion at the lower level
    real(r8kind),      intent(in) :: rrdef1     ! Rossby radius of deformation of 200-500 thickness
    real(r8kind),      intent(in) :: rrdef2     ! Rossby radius of deformation of 500-800 thickness
    type(qg_config_type)          :: qg_config

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


end module QG_Config
