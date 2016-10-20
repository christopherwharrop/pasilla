module background

  use kind, only             : r8kind
  use module_constants, only : PI
  use lorenz96, only         : lorenz96_type

  implicit none

  private

  public :: Background_Type

  type Background_Type
      private
      integer, public                   :: npoints
      integer, public                   :: ntimes
      real(r8kind), allocatable, public :: state(:,:)
      integer, allocatable, public      :: position(:,:)
      integer, allocatable, public      :: time(:)
! instance variable
  contains
      final :: destructor
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
  type(Background_Type) function constructor(ntimes, method)

    integer, intent(in) :: ntimes
    integer, intent(in) :: method

    type(lorenz96_type) :: model
    integer             :: i, t, tt

    constructor%ntimes = ntimes

    ! Allocate time array
    allocate (constructor%time(ntimes))

    ! Read in model background for each time
    do t = 1, ntimes
       constructor%time(t) = t
       tt = t
       if (method .le. 2) tt = 2
       model = lorenz96_type(tt, "NETCDF")
       if (t==1) then
         constructor%npoints = model%size
         allocate (constructor%state(ntimes, model%size))
         allocate (constructor%position(ntimes, model%size))
       end if
       constructor%position(t,:) = model%location(:)
       constructor%state(t,:) = model%state(:)
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


end module background
