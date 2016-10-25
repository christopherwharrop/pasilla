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
      integer                   :: size
      integer                   :: ntimes
      real(r8kind)              :: var
      real(r8kind), allocatable :: covariance(:,:,:)
  contains
      ! methods
      final :: destructor
      procedure :: get_size
      procedure :: get_covariances_at_time
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

    constructor%size = background%get_npoints()
    constructor%ntimes = cfg%ntimes
    constructor%var = var

    ! Allocate object arrays
    allocate (constructor%covariance(cfg%ntimes, constructor%size, constructor%size))

    constructor%covariance(:,:,:) = 0.0
    rad = constructor%size / 10

    do t = 1, cfg%ntimes
       do i = 1, constructor%size
          do j = -rad, +rad
             jj = i + j
             if (jj > constructor%size) jj = jj - constructor%size
             if (jj < 1) jj = constructor%size + jj
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


  !------------------------------------------------------------------
  ! get_size
  !
  ! Return the size of the background covariance matrix
  !------------------------------------------------------------------
  integer function get_size(this)

    class(Background_Covariance_Type) :: this

    get_size = this%size

  end function get_size


  !------------------------------------------------------------------
  ! get_covariances_at_time
  !
  ! Return the covariance matrix for time t
  !------------------------------------------------------------------
  function get_covariances_at_time(this, t)

    class(Background_Covariance_Type) :: this
    integer                           :: t

    real(r8kind), dimension(this%size,this%size) :: get_covariances_at_time

    get_covariances_at_time = this%covariance(t,:,:)

  end function get_covariances_at_time


end module background_covariance
