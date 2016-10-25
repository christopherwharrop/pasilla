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
      integer                           :: size
      integer                           :: ntimes
      real(r8kind), allocatable, public :: covariance(:,:,:)
  contains
      ! methods
      final :: destructor
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
    constructor%ntimes = cfg%ntimes
    
    allocate(constructor%covariance(cfg%ntimes, nobs, nobs))

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


end module observation_covariance
