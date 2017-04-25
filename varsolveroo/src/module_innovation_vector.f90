module innovation_vector

  use kind, only         : r8kind
  use background, only   : Background_Type
  use observations, only : Observations_Type
  use observation_operator, only : Observation_Operator_type

  implicit none

  private

  public :: Innovation_Vector_Type

  type Innovation_Vector_Type
      ! instance variable
      private
      integer                   :: nobs
      real(r8kind), allocatable :: value(:)
      integer, allocatable      :: position(:)
      integer, allocatable      :: time(:)
  contains
      ! methods
      procedure :: get_value_vector
      final :: destructor
  end type Innovation_Vector_Type

  interface Innovation_Vector_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Innovation_vector
  !------------------------------------------------------------------
  type(Innovation_vector_Type) function constructor(background, observations, H)

    class(Background_Type), intent(in)           :: background
    class(Observations_Type), intent(in)         :: observations
    class(Observation_Operator_Type), intent(in) :: H

    integer :: i, t, nobs

    nobs = observations%get_nobs()

    ! Initialize the length of the innovation vector
    constructor%nobs = nobs

    ! Allocate space for the innovation vector values, positions, and times
    allocate(constructor%value(nobs))
    allocate(constructor%position(nobs))
    allocate(constructor%time(nobs))

    ! Initialize the innovation vector positions and times
    constructor%position = observations%get_position_vector()
    constructor%time = observations%get_time_vector()

    constructor%value = observations%get_value_vector()
    do t = 1, background%get_ntimes()
      constructor%value = constructor%value - matmul(H%get_operator_at_time(t), background%get_state_at_time(t))
    end do

    ! Calculate the innovation vector values
    do i = 1, nobs
       write(*,'(A8,3I5,2F10.4)') "INO ", i, constructor%time(i), constructor%position(i), constructor%value(i), background%get_state_element(constructor%time(i), constructor%position(i))
    end do

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Innovation_Vector object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Innovation_Vector_Type), intent(inout) :: this

    ! No pointers in Innovation_vector object so we do nothing

  end subroutine


  function get_value_vector(this)

    class(Innovation_Vector_Type) :: this

    real(r8kind), dimension(this%nobs) :: get_value_vector

    get_value_vector = this%value

  end function get_value_vector


end module innovation_vector
