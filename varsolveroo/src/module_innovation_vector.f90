module innovation_vector

  use kind, only         : r8kind
  use background, only   : Background_Type
  use observations, only : Observations_Type

  implicit none

  private

  public :: Innovation_Vector_Type

  type Innovation_Vector_Type
      private
      integer, public                   :: nobs
      real(r8kind), allocatable, public :: value(:)
      integer, allocatable, public      :: position(:)
      integer, allocatable, public      :: time(:)
! instance variable
  contains
      final :: destructor
! methods
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
  type(Innovation_vector_Type) function constructor(background, observations)

    class(Background_Type), intent(in)   :: background
    class(Observations_Type), intent(in) :: observations

    integer :: i

    ! Initialize the length of the innovation vector
    constructor%nobs = observations%nobs

    ! Allocate space for the innovation vector values, positions, and times
    allocate(constructor%value(observations%nobs))
    allocate(constructor%position(observations%nobs))
    allocate(constructor%time(observations%nobs))

    ! Initialize the innovation vector positions and times
    constructor%position = observations%position
    constructor%time = observations%time

    ! Calculate the innovation vector values
    do i = 1, observations%nobs
       constructor%value(i) = observations%value(i) - background%state(observations%time(i), observations%position(i))
       write(*,'(A8,3I5,2F10.4)') "INO ", i, constructor%time(i), constructor%position(i), constructor%value(i), background%state(constructor%time(i), constructor%position(i))
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


end module innovation_vector
