module module_background

  use kind, only      : r8kind
  use module_constants, only : PI

  implicit none

  private

  public :: Background_Type

  type Background_Type
      private
      integer                           :: npoints
      integer                           :: ntimes
      real(r8kind), allocatable, public :: state(:,:)
      integer, allocatable, public      :: position(:,:)
      integer, allocatable, public      :: time(:)
! instance variable
  contains
      final              :: destructor
! methods
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
  type(Background_Type) function constructor(npoints, ntimes, method)

    integer, intent(in) :: npoints
    integer, intent(in) :: ntimes
    integer, intent(in) :: method

    integer :: i, t, tt

    constructor%npoints = npoints
    constructor%ntimes = ntimes

    ! Allocate object arrays
    allocate (constructor%state(ntimes,npoints))
    allocate (constructor%position(ntimes,npoints))
    allocate (constructor%time(ntimes))

    do t = 1, ntimes
       constructor%time(t) = t
       tt = t
       if (method .le. 2) tt = 2
       do i = 1, npoints
          constructor%position(t,i) = i
          constructor%state(t,i) = 50.0 + 50.0 * sin(((20.0 * float(tt - 2) + float(i)) / 1000.0) * PI)
       end do
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


end module module_background
