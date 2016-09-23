module observations

  use kind, only      : r8kind
  use module_constants, only : PI

  implicit none

  private

  public :: Observations_Type

  type Observations_Type
      private
      integer, public                   :: nobs
      real(r8kind), allocatable, public :: value(:)
      integer, allocatable, public      :: position(:)
      integer, allocatable, public      :: time(:)
! instance variable
  contains
      final              :: destructor
! methods
  end type Observations_Type

  interface Observations_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Observations
  !------------------------------------------------------------------
  type(Observations_Type) function constructor(nobs, method)

    integer, intent(in) :: nobs
    integer, intent(in) :: method

    integer :: i, mod

    constructor%nobs = nobs

    ! Allocate object arrays
    allocate(constructor%value(nobs))
    allocate(constructor%position(nobs))
    allocate(constructor%time(nobs))

    do i = 1, nobs

        mod = modulo(i, 8)
        constructor%time(i) = 2
        constructor%position(i) = i * 200

        if (mod .gt. 3) then
            constructor%time(i) = 1
            constructor%position(i) = i * 200 - 35
        end if
        if (mod .gt. 5) then
            constructor%time(i) = 3
            constructor%position(i) = i * 200 - 75
        end if

        ! FOR MATCHED 3DVAR
        if (method .eq. 2) constructor%time(i) = 2

        constructor%value(i) = 50.0 + 50.0 * sin(((20.0 * float(constructor%time(i) - 1) + float(constructor%position(i))) / 1000.0) * PI)

        ! FOR 3DVAR
        if (method .le. 2) constructor%time(i) = 1

    end do

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Observations object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Observations_Type), intent(inout) :: this

    ! No pointers in Observations object so we do nothing

  end subroutine


end module observations
