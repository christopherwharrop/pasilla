module observation_operator

  use kind, only         : r8kind
  use background, only   : Background_Type
  use observations, only : Observations_Type
  use L96_Model,  only   : l96_model_type

  implicit none

  private

  public :: Observation_Operator_Type

  type Observation_Operator_Type
      private
      ! instance variable
      integer                   :: nobs
      integer                   :: npoints
      integer                   :: ntimes
      real(r8kind), allocatable :: operator(:,:,:)
  contains
      ! methods
      final :: destructor
      procedure :: get_operator_at_time
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

    integer :: i, nobs
    integer :: lower_index, upper_index
    real(KIND=8) :: lower_weight, upper_weight
    type(L96_Model_Type) :: model

    nobs = obs%get_nobs()

    ! Initialize the dimensions of the operator
    constructor%nobs = nobs
    constructor%npoints = bkg%get_npoints()
    constructor%ntimes = bkg%get_ntimes()

    ! Allocate space for the observation operator matrix
    allocate(constructor%operator(constructor%ntimes, nobs, constructor%npoints))

    ! Initialize the matrix to zero
    constructor%operator(:,:,:) = 0.0

    ! Instantiate a L96 model to get interpolation weights from
    model = l96_model_type(bkg%get_config())

    ! Initialize non-zero values
    do i = 1, nobs
       call model%get_interpolation_weights(obs%get_position_element(i), lower_index, upper_index, lower_weight, upper_weight)
       constructor%operator(obs%get_time_element(i), i, lower_index) = lower_weight
       constructor%operator(obs%get_time_element(i), i, upper_index) = upper_weight
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


  !------------------------------------------------------------------
  ! get_operator_at_time
  !
  ! Get the operator matrix at time t
  !------------------------------------------------------------------
  function get_operator_at_time(this, t)

    class(Observation_Operator_Type), intent(in) :: this
    integer, intent(in)                          :: t

    real(r8kind), dimension(this%nobs, this%npoints) :: get_operator_at_time

    get_operator_at_time = this%operator(t,:,:)

  end function get_operator_at_time


end module observation_operator
