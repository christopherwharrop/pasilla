module background

  use kind, only             : r8kind
  use module_constants, only : PI
  use config, only           : Config_Type
  use lorenz96, only         : lorenz96_type

  implicit none

  private

  public :: Background_Type

  type Background_Type
      private
      ! instance variable
      integer, public                   :: npoints
      integer, public                   :: ntimes
      real(r8kind), allocatable, public :: state(:,:)
      integer, allocatable, public      :: position(:,:)
      integer, allocatable, public      :: time(:)
  contains
      ! methods
      final :: destructor
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

    type(lorenz96_type) :: model
    integer             :: t, tt

    constructor%ntimes = cfg%ntimes

    ! Allocate time array
    allocate (constructor%time(cfg%ntimes))

    ! Read in model background for each time
    do t = 1, cfg%ntimes
       constructor%time(t) = t
       tt = t
       if (cfg%method .le. 2) tt = 2
       model = lorenz96_type(tt, "NETCDF")
       if (t==1) then
         constructor%npoints = model%size
         allocate (constructor%state(cfg%ntimes, model%size))
         allocate (constructor%position(cfg%ntimes, model%size))
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
