module observations

  use kind, only   : r8kind
  use config, only : Config_Type

  implicit none

  private

  public :: Observations_Type

  type Observations_Type
      private
      ! instance variable
      integer                   :: nobs
      real(r8kind), allocatable :: value(:)
      real(r8kind), allocatable :: position(:)
      integer, allocatable      :: time(:)
  contains
      ! methods
      final     :: destructor
      procedure :: get_nobs
      procedure :: get_time_element
      procedure :: get_time_vector
      procedure :: get_position_element
      procedure :: get_position_vector
      procedure :: get_value_element
      procedure :: get_value_vector
  end type Observations_Type

  interface Observations_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Observations
  !------------------------------------------------------------------
  type(Observations_Type) function constructor(cfg)

    class(Config_Type), intent(in) :: cfg

    integer            :: i
    character(len=128) :: filename   ! name of output file
    integer            :: fileunit


    ! Construct name of obs input file
    write(filename, '(A,I1,A)') 'lorenz96obs_', cfg%get_method(), '.txt'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted', status='old')

    ! Read the number of obs
    read(fileunit, '(I)') constructor%nobs

    ! Allocate object arrays
    allocate(constructor%value(constructor%nobs))
    allocate(constructor%position(constructor%nobs))
    allocate(constructor%time(constructor%nobs))

    do i=1,constructor%nobs
      read(fileunit, '(I,2F12.1)') constructor%time(i), constructor%position(i), constructor%value(i)
    end do

    close(fileunit)

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Observations object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Observations_Type), intent(inout) :: this

    ! No pointers in Observations object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! get_nobs
  !
  ! Return the number of obs
  !------------------------------------------------------------------
  integer function get_nobs(this)

    class(Observations_Type) :: this

    get_nobs = this%nobs

  end function get_nobs


  !------------------------------------------------------------------
  ! get_time_element
  !
  ! Return the time of observation i
  !------------------------------------------------------------------
  function get_time_element(this, i)

    class(Observations_Type) :: this
    integer :: i

    integer :: get_time_element

    get_time_element = this%time(i)

  end function get_time_element


  !------------------------------------------------------------------
  ! get_time_vector
  !
  ! Return the time vector
  !------------------------------------------------------------------
  function get_time_vector(this)

    class(Observations_Type) :: this

    integer, dimension(this%nobs) :: get_time_vector

    get_time_vector = this%time

  end function get_time_vector


  !------------------------------------------------------------------
  ! get_position_vector
  !
  ! Return the observation position vector
  !------------------------------------------------------------------
  function get_position_vector(this)

    class(Observations_Type) :: this

    integer, dimension(this%nobs) :: get_position_vector

    get_position_vector = this%position

  end function get_position_vector


  !------------------------------------------------------------------
  ! get_position_element
  !
  ! Return the position of observation i
  !------------------------------------------------------------------
  function get_position_element(this, i)

    class(Observations_Type) :: this
    integer                  :: i

    integer :: get_position_element

    get_position_element = this%position(i)

  end function get_position_element


  !------------------------------------------------------------------
  ! get_value_element
  !
  ! Return the value of observation i
  !------------------------------------------------------------------
  function get_value_element(this, i)

    class(Observations_Type) :: this
    integer                  :: i

    real(r8kind) :: get_value_element

    get_value_element = this%value(i)

  end function get_value_element


  !------------------------------------------------------------------
  ! get_value_vector
  !
  ! Return the observation values
  !------------------------------------------------------------------
  function get_value_vector(this) result(obs_vector)

    class(Observations_Type)           :: this
    real(r8kind), dimension(this%nobs) :: obs_vector

    obs_vector = this%value

  end function get_value_vector


end module observations
