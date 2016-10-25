module observation_covariance

  use kind, only         : r8kind
  use config, only       : Config_Type
  use observations, only : Observations_Type

  implicit none

  private

  public :: Observation_Covariance_Type

  type Observation_Covariance_Type
      private
      ! instance variable
      integer                   :: size
      integer                   :: ntimes
      real(r8kind), allocatable :: covariance(:,:,:)
  contains
      ! methods
      final :: destructor
      procedure :: get_size
      procedure :: get_covariances_at_time
  end type Observation_Covariance_Type

  interface Observation_Covariance_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Observation_Covariance
  !------------------------------------------------------------------
  type(Observation_Covariance_Type) function constructor(observations, cfg)

    class(Observations_Type), intent(in) :: observations
    class(Config_Type),       intent(in) :: cfg

    integer :: i, nobs

    nobs = observations%get_nobs()

    constructor%size = nobs
    constructor%ntimes = cfg%get_ntimes()
    
    allocate(constructor%covariance(cfg%get_ntimes(), nobs, nobs))

    constructor%covariance(:,:,:) = 0.0

    do i = 1, nobs
       constructor%covariance(observations%get_time_element(i), i, i) = 1.0
    end do

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Observation_Covariance object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Observation_Covariance_Type), intent(inout) :: this

    ! No pointers in Observation_Covariance object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! get_size
  !
  ! Get the size of the observation covariance matrix
  !------------------------------------------------------------------
  integer function get_size(this)

    class(Observation_Covariance_Type) :: this

    get_size = this%size

  end function get_size


  !------------------------------------------------------------------
  ! get_covariances_at_time
  !
  ! Get the covariance matrix at time t
  !------------------------------------------------------------------
  function get_covariances_at_time(this, t)

    class(Observation_Covariance_Type) :: this
    integer                            :: t

    real(r8kind), dimension(this%size, this%size) :: get_covariances_at_time

    get_covariances_at_time = this%covariance(t,:,:)

  end function get_covariances_at_time


end module observation_covariance
