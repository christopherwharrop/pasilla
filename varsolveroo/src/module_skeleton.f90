module module_skeleton

  use module_kind, only : r8kind

  implicit none

  private

  public :: object

  type object
      private
! instance variable
  contains
      final              :: destructor
! methods
  end type object

  interface object
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized object
  !------------------------------------------------------------------
  type(object) function constructor()

    ! Initialize object
!    constructor%variable = something

    ! Allocate object arrays
!    allocate(constructor%array(size))

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a lorenz96 object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(object), intent(inout) :: this

    ! No pointers in lorenz96 object so we do nothing

  end subroutine


end module module_skeleton
