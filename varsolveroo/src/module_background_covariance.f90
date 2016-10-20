module background_covariance

  use kind, only      : r8kind
  use module_constants, only : PI
  use background, only : Background_Type
  implicit none

  private

  public :: Background_Covariance_Type

  type Background_Covariance_Type
      private
      integer                           :: size
      integer                           :: ntimes
      real(r8kind)                      :: var
      real(r8kind), allocatable, public :: covariance(:,:,:)
! instance variable
  contains
      final :: destructor
! methods
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
  type(Background_Covariance_Type) function constructor(background, ntimes, sigma)

    class(Background_Type), intent(in) :: background
    integer, intent(in) :: ntimes
    real(r8kind), intent(in) :: sigma

    real(r8kind) :: var
    integer :: t, i, j, jj, rad 

    var = 3.61

    constructor%size = background%npoints
    constructor%ntimes = ntimes
    constructor%var = var

    ! Allocate object arrays
    allocate (constructor%covariance(ntimes, background%npoints, background%npoints))

    constructor%covariance(:,:,:) = 0.0
    rad = background%npoints / 10

    do t = 1, ntimes
       do i = 1, background%npoints
          do j = -rad, +rad
             jj = i + j
             if (jj > background%npoints) jj = jj - background%npoints
             if (jj < 1) jj = background%npoints + jj
             constructor%covariance(t, i, jj) = var * exp(-((float(j) * sigma)**2))
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
