module L96_State

  use kind,  only : r8kind
!  use State, only : state_type

  implicit none

  private

  public :: l96_state_type

  type :: l96_state_type
      private
      integer :: size
      real(r8kind), allocatable :: state(:)
  contains
      final :: destructor
      procedure :: get_size
      procedure :: get_state
  end type l96_state_type

  interface l96_state_type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized l96_state_type object
  !------------------------------------------------------------------
  type(l96_state_type) function constructor(state)

    real(r8kind) :: state(:)

    ! Initialize model state
    constructor%state(:) = state(:)
    constructor%size = size(state)

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a l96_state_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(l96_state_type), intent(inout) :: this

    ! No pointers in l96_state_type object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! get_size
  !
  ! Returns the size of the state array
  !------------------------------------------------------------------
  integer function get_size(this)

    class(L96_State_Type) :: this

    get_size = this%size

  end function get_size


  !------------------------------------------------------------------
  ! get_state
  !
  ! Returns the state array
  !------------------------------------------------------------------
  function get_state(this)

    class(L96_State_Type) :: this
    real(r8kind), dimension(this%size) :: get_state

    get_state = this%state

  end function get_state
 

end module L96_State
