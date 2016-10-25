module config

  use kind, only : r8kind

  implicit none

  private

  public :: Config_Type

  type Config_Type
      private
      ! instance variable
      integer      :: method
      integer      :: ntimes
      real(r8kind) :: sigma
      real(r8kind) :: alpha
  contains
      ! methods
      final :: destructor
      procedure :: get_method
      procedure :: get_ntimes
      procedure :: get_sigma
      procedure :: get_alpha
  end type Config_Type

  interface Config_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Config
  !------------------------------------------------------------------
  type(Config_Type) function constructor(nlunit)

    integer :: nlunit

    ! Define namelists and default values
    integer          :: mthd = 4
    integer          :: tim_len = 3
    real(KIND=8)     :: alph = 1.0
    real(KIND=8)     :: sigma = 1.0
    namelist /control/  mthd, tim_len
    namelist /method1/  alph, sigma
    namelist /method2/  alph, sigma
    namelist /method3/  alph, sigma
    namelist /method4/  alph, sigma

    print *,"GET_METHOD"

    ! Read namelists
    read(nlunit,nml=control)

    select case (mthd)
      case(1)
        read(nlunit, nml=method1)
      case(2)
        read(nlunit, nml=method2)
      case(3)
        read(nlunit, nml=method3)
      case(4)
        read(nlunit, nml=method4)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: method "',mthd,'" is not supported!'
        stop
    end select

    print *,"METHOD = ", mthd
    print *,"NTIMES = ", tim_len
    print *,"ALPH = ", alph
    print *,"SIGMA = ", sigma

    ! Force ntimes=1 for 3DVAR
    if(mthd <= 2) tim_len = 1

    constructor%method = mthd
    constructor%ntimes = tim_len
    constructor%alpha = alph
    constructor%sigma = sigma

    print *,"GET_METHOD COMPLETE"


  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Config object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Config_Type), intent(inout) :: this

    ! No pointers in Config object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! get_method
  !
  ! Get the method from the configuration
  !------------------------------------------------------------------
  integer function get_method(this)

    class(Config_Type), intent(in) :: this

    get_method = this%method

  end function get_method


  !------------------------------------------------------------------
  ! get_ntimes
  !
  ! Get the ntimes from the configuration
  !------------------------------------------------------------------
  integer function get_ntimes(this)

    class(Config_Type), intent(in) :: this

    get_ntimes = this%ntimes

  end function get_ntimes


  !------------------------------------------------------------------
  ! get_alpha
  !
  ! Get the alpha from the configuration
  !------------------------------------------------------------------
  real(r8kind) function get_alpha(this)

    class(Config_Type), intent(in) :: this

    get_alpha = this%alpha

  end function get_alpha


  !------------------------------------------------------------------
  ! get_sigma
  !
  ! Get the sigma from the configuration
  !------------------------------------------------------------------
  real(r8kind) function get_sigma(this)

    class(Config_Type), intent(in) :: this

    get_sigma = this%sigma

  end function get_sigma


end module config
