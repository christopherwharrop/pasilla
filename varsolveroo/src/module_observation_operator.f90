module observation_operator

  use kind, only         : r8kind
  use background, only   : Background_Type
  use observations, only : Observations_Type

  implicit none

  private

  public :: Observation_Operator_Type

  type Observation_Operator_Type
      private
      integer, public                   :: nobs
      integer, public                   :: npoints
      integer, public                   :: ntimes
      real(r8kind), allocatable, public :: operator(:,:,:)
! instance variable
  contains
      final :: destructor
! methods
  end type Observation_Operator_Type

  interface Observation_Operator_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Observation_Operator
  !------------------------------------------------------------------
  type(Observation_Operator_Type) function constructor(bkg, obs)

    class(Background_Type), intent(in)   :: bkg
    class(Observations_Type), intent(in) :: obs

    integer :: i

    ! Initialize the dimensions of the operator
    constructor%nobs = obs%nobs
    constructor%npoints = bkg%npoints
    constructor%ntimes = bkg%ntimes

    ! Allocate space for the observation operator matrix
    allocate(constructor%operator(bkg%ntimes, obs%nobs, bkg%npoints))

    ! Initialize the matrix to zero
    constructor%operator(:,:,:) = 0.0

    ! Initialize non-zero values
    do i = 1, obs%nobs
       constructor%operator(obs%time(i), i, obs%position(i)) = 1.0
    end do

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Observation_Operator object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Observation_Operator_Type), intent(inout) :: this

    ! No pointers in Observation_Operator object so we do nothing

  end subroutine


end module observation_operator
