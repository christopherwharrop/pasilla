module observation_covariance

  use kind, only      : r8kind
  use module_constants, only : PI
  use observations, only : Observations_Type

  implicit none

  private

  public :: Observation_Covariance_Type

  type Observation_Covariance_Type
      private
      integer                           :: size
      integer                           :: ntimes
      real(r8kind), allocatable, public :: covariance(:,:,:)
! instance variable
  contains
      final :: destructor
! methods
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
  type(Observation_Covariance_Type) function constructor(observations, ntimes)

    class(Observations_Type), intent(in) :: observations
    integer, intent(in) :: ntimes

    integer :: i

    constructor%size = observations%nobs
    constructor%ntimes = ntimes
    
    allocate(constructor%covariance(ntimes, observations%nobs, observations%nobs))

    constructor%covariance(:,:,:) = 0.0

    do i = 1, observations%nobs
       constructor%covariance(observations%time(i), i, i) = 1.0
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
