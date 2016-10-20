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
  type(Observations_Type) function constructor(method)

    integer, intent(in) :: method

    integer            :: i
    character(len=128) :: filename   ! name of output file
    integer            :: fileunit


    ! Construct name of obs input file
    write(filename, '(A,I1,A)') 'lorenz96obs_', method, '.txt'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted', status='old')

    ! Read the number of obs
    read(fileunit, '(I)') constructor%nobs

    ! Allocate object arrays
    allocate(constructor%value(constructor%nobs))
    allocate(constructor%position(constructor%nobs))
    allocate(constructor%time(constructor%nobs))

    do i=1,constructor%nobs
      read(fileunit, '(2I,F)') constructor%time(i), constructor%position(i), constructor%value(i)
    end do

    close(fileunit)

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
