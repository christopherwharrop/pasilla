module L96_Config

  use kind,            only : r8kind
  use Abstract_Config, only : abstract_config_type

  implicit none

  private

  public :: l96_config_type

  type, extends(abstract_config_type) :: l96_config_type
      private
      integer      :: nx
      real(r8kind) :: time_step
      real(r8kind) :: forcing
  contains
      final :: destructor
      procedure :: get_nx
      procedure :: get_time_step
      procedure :: get_forcing
      procedure :: print
  end type l96_config_type

  interface l96_config_type
    procedure constructor_arglist
    procedure constructor_namelist_file
    procedure constructor_namelist_unit
  end interface

contains

  !------------------------------------------------------------------
  ! constructor_arglist
  !
  ! Returns an initialized l96_config_type object
  !------------------------------------------------------------------
  type(l96_config_type) function constructor_arglist(nx, time_step, forcing)

    integer,      intent(in) :: nx
    real(r8kind), intent(in) :: time_step
    real(r8kind), intent(in) :: forcing

    ! Initialize model parameters    
    constructor_arglist%nx = nx
    constructor_arglist%time_step = time_step
    constructor_arglist%forcing = forcing

  end function constructor_arglist


  !------------------------------------------------------------------
  ! constructor_namelist_file
  !
  ! Returns an initialized l96_config_type object
  !------------------------------------------------------------------
  type(l96_config_type) function constructor_namelist_file(nl_filename)

    ! Namelist filename
    character(len=*) :: nl_filename

    ! Namelist file descriptor
    integer :: nl_unit

    ! Define namelists and default values
    integer      :: nx = 40
    real(r8kind) :: time_step = 0.05_r8kind
    real(r8kind) :: forcing = 8.00_r8kind
    namelist /params/ nx, forcing, time_step

    ! Open the namelist file
    open(newunit=nl_unit, file=trim(nl_filename), form='formatted', status='old')

    ! Read the configuration
    read(nl_unit, nml=params)

    ! Initialize model parameters    
    constructor_namelist_file%nx = nx
    constructor_namelist_file%time_step = time_step
    constructor_namelist_file%forcing = forcing

    ! Close the namelist
    close(nl_unit)

  end function constructor_namelist_file


  !------------------------------------------------------------------
  ! constructor_namelist_unit
  !
  ! Returns an initialized l96_config_type object
  !------------------------------------------------------------------
  type(l96_config_type) function constructor_namelist_unit(nl_unit)

    ! Namelist unit number (must be an open file)
    integer :: nl_unit

    ! Define namelists and default values
    integer      :: nx = 40
    real(r8kind) :: time_step = 0.05_r8kind
    real(r8kind) :: forcing = 8.00_r8kind
    namelist /params/ nx, forcing, time_step

    ! Read the configuration
    read(nl_unit, nml=params)

    ! Initialize model parameters    
    constructor_namelist_unit%nx = nx
    constructor_namelist_unit%time_step = time_step
    constructor_namelist_unit%forcing = forcing

    ! Close the namelist
    rewind(nl_unit)

  end function constructor_namelist_unit


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a l96_config_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(l96_config_type), intent(inout) :: this

    ! No pointers in l96_config_type object so we do nothing

  end subroutine destructor


  !------------------------------------------------------------------
  ! get_nx
  !
  ! Returns the length of the x dimension
  !------------------------------------------------------------------
  pure integer function get_nx(this)

    class(L96_Config_Type), intent(in) :: this

    get_nx = this%nx

  end function get_nx


  !------------------------------------------------------------------
  ! get_forcing
  !
  ! Returns the value of the forcing parameter
  !------------------------------------------------------------------
  pure real(r8kind) function get_forcing(this)

    class(L96_Config_Type), intent(in) :: this

    get_forcing = this%forcing

  end function get_forcing


  !------------------------------------------------------------------
  ! get_time_step
  !
  ! Returns the value of the delta_t parameter
  !------------------------------------------------------------------
  pure real(r8kind) function get_time_step(this)

    class(L96_Config_Type), intent(in) :: this

    get_time_step = this%time_step

  end function get_time_step


  !------------------------------------------------------------------
  ! print
  !
  ! Prints contents of config object
  !------------------------------------------------------------------
  subroutine print(this)

    class(L96_Config_Type), intent(in) :: this

    character(len=32) :: numstr

    write(numstr,'(I)') this%nx
    write(*,'(A20,A)') "nx = ", adjustl(numstr)
    write(numstr,'(F16.13)') this%time_step
    write(*,'(A20,A)') "time_step = ", adjustl(numstr)
    write(numstr,'(F16.13)') this%forcing
    write(*,'(A20,A)') "forcing = ", adjustl(numstr)

  end subroutine print


end module L96_Config
