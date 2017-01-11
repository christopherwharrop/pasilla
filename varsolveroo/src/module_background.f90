module background

  use kind, only             : r8kind
  use module_constants, only : PI
  use config, only           : Config_Type
  use L96_Config, only       : l96_config_type
  use L96_Model,  only       : l96_model_type
  use L96_Reader, only       : l96_reader_type

  implicit none

  private

  public :: Background_Type

  type Background_Type
      private
      ! instance variable
      type(l96_config_type)     :: config
      integer                   :: npoints
      integer                   :: ntimes
      integer                   :: step
      real(r8kind), allocatable :: state(:,:)
      integer, allocatable      :: position(:,:)
      integer, allocatable      :: time(:)
  contains
      ! methods
      final :: destructor
      procedure :: get_config
      procedure :: get_npoints
      procedure :: get_ntimes
      procedure :: get_step
      procedure :: get_state_element
      procedure :: get_state_vector
      procedure :: get_state_at_time
  end type Background_Type

  interface Background_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Background
  !------------------------------------------------------------------
  type(Background_Type) function constructor(cfg)

    class(Config_Type), intent(in) :: cfg

    type(l96_model_type)  :: model
    type(l96_config_type) :: model_config
    type(l96_reader_type) :: reader
    character(len=256)    :: filename
    integer               :: t, tt

    constructor%ntimes = cfg%get_ntimes()

    reader = l96_reader_type("NETCDF")

    ! Allocate time array
    allocate (constructor%time(cfg%get_ntimes()))

    ! Read in model background for each time
    do t = 1, cfg%get_ntimes()
       constructor%time(t) = t
       tt = t
       if (cfg%get_method() .le. 2) tt = 2
       write(filename,'(A,I0.7)') 'bkgin_', tt
       call reader%read(model, filename)
       model_config = model%get_config()
       constructor%config = model_config
       if (t==1) then
         constructor%npoints = model_config%get_nx()
         allocate (constructor%state(constructor%ntimes, constructor%npoints))
         allocate (constructor%position(constructor%ntimes, constructor%npoints))
       end if
       if (tt == 2) constructor%step = model%get_step()
!       constructor%position(t,:) = model%location(:)
       constructor%state(t,:) = model%get_state()
    end do

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Background object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Background_Type), intent(inout) :: this

    ! No pointers in Background object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! get_config
  !
  ! Returns the model configuration associated with the background
  !------------------------------------------------------------------
  type(l96_config_type) function get_config(this)

    class(Background_Type) :: this

    get_config = this%config

  end function get_config


  !------------------------------------------------------------------
  ! get_npoints
  !
  ! Returns the size of the state vector
  !------------------------------------------------------------------
  integer function get_npoints(this)

    class(Background_Type) :: this

    get_npoints = this%npoints

  end function get_npoints


  !------------------------------------------------------------------
  ! get_ntimes
  !
  ! Returns the size of the time dimension of the state vector
  !------------------------------------------------------------------
  integer function get_ntimes(this)

    class(Background_Type) :: this

    get_ntimes = this%ntimes

  end function get_ntimes


  !------------------------------------------------------------------
  ! get_step
  !
  ! Returns the model time step for this background
  !------------------------------------------------------------------
  integer function get_step(this)

    class(Background_Type) :: this

    get_step = this%step

  end function get_step


  !------------------------------------------------------------------
  ! get_state_element
  !
  ! Returns the state at element t,i
  !------------------------------------------------------------------
  real(r8kind) function get_state_element(this, t, i)

    class(Background_Type) :: this
    integer                :: t, i

    get_state_element = this%state(t,i)

  end function get_state_element


  !------------------------------------------------------------------
  ! get_state_vector
  !
  ! Returns the state vector
  !------------------------------------------------------------------
  function get_state_vector(this)

    class(Background_Type) :: this
    real(r8kind), dimension(this%ntimes, this%npoints) :: get_state_vector

    get_state_vector = this%state

  end function get_state_vector


  !------------------------------------------------------------------
  ! get_state_at_time
  !
  ! Returns the state vector at time t
  !------------------------------------------------------------------
  function get_state_at_time(this, t)

    class(Background_Type) :: this
    integer                :: t
    real(r8kind), dimension(this%npoints) :: get_state_at_time

    get_state_at_time = this%state(t,:)

  end function get_state_at_time


end module background
