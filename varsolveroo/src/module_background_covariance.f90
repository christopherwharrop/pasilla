module background_covariance

  use kind, only       : r8kind
  use background, only : Background_Type
  use config, only     : Config_Type

  implicit none

  private

  public :: Background_Covariance_Type

  type Background_Covariance_Type
      private
      ! instance variable
      integer                           :: size
      integer                           :: ntimes
      real(r8kind)                      :: var
      real(r8kind), allocatable, public :: covariance(:,:,:)
  contains
      ! methods
      final :: destructor
  end type Background_Covariance_Type

  interface Background_Covariance_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Background_Covariance
  !------------------------------------------------------------------
  type(Background_Covariance_Type) function constructor(background, cfg)

    class(Background_Type), intent(in) :: background
    class(Config_Type),     intent(in) :: cfg

    real(r8kind) :: var
    integer :: t, i, j, jj, rad 

    var = 3.61

    constructor%size = background%npoints
    constructor%ntimes = cfg%ntimes
    constructor%var = var

    ! Allocate object arrays
    allocate (constructor%covariance(cfg%ntimes, background%npoints, background%npoints))

    constructor%covariance(:,:,:) = 0.0
    rad = background%npoints / 10

    do t = 1, cfg%ntimes
       do i = 1, background%npoints
          do j = -rad, +rad
             jj = i + j
             if (jj > background%npoints) jj = jj - background%npoints
             if (jj < 1) jj = background%npoints + jj
             constructor%covariance(t, i, jj) = var * exp(-((float(j) * cfg%sigma)**2))
          end do
       end do
    end do

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Background_Covariance object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Background_Covariance_Type), intent(inout) :: this

    ! No pointers in Background_Covariance object so we do nothing

  end subroutine


end module background_covariance
