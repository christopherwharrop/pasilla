module NetCDF_Writer

  use netcdf
  use kind,      only  : r8kind

  implicit none

  private

  public :: netcdf_writer_type

  type :: netcdf_writer_type
      private
      integer :: ncFileID      ! netCDF file identifier
      integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
  contains
      final :: destructor
      procedure :: create
      procedure :: close
      generic   :: write_global_var => write_global_string, write_global_integer, write_global_real
      procedure, private :: write_global_string
      procedure, private :: write_global_integer
      procedure, private :: write_global_real
      procedure, private, nopass :: nc_check
  end type netcdf_writer_type

  interface netcdf_writer_type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized l96_writer_type object
  !------------------------------------------------------------------
  type(netcdf_writer_type) function constructor()

    constructor%ncFileID=0
    constructor%nDimensions=0
    constructor%nVariables=0
    constructor%nAttributes=0
    constructor%unlimitedDimID=0

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a l96_writer_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(netcdf_writer_type), intent(inout) :: this

    ! No pointers in netcdf_writer_type object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! Create
  !
  ! Open new file, overwriting previous contents if it exists
  !------------------------------------------------------------------
  subroutine create(this, filename)

    class(netcdf_writer_type), intent(inout) :: this
    character(len=*),          intent(in) :: filename

    integer :: ierr

    call nc_check(nf90_create(trim(filename) // '.nc', NF90_CLOBBER, this%ncFileID))
    call nc_check(nf90_Inquire(this%ncFileID, this%nDimensions, this%nVariables, this%nAttributes, this%unlimitedDimID))

  end subroutine create


  !------------------------------------------------------------------
  ! Write_global_string
  !
  ! Write a global string variable
  !------------------------------------------------------------------
  subroutine write_global_string(this, label, string_value)

    class(netcdf_writer_type), intent(in) :: this
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: string_value

    call nc_check(nf90_put_att(this%ncFileID, NF90_GLOBAL, label, string_value))

  end subroutine write_global_string


  !------------------------------------------------------------------
  ! Write_global_integer
  !
  ! Write a global integer variable
  !------------------------------------------------------------------
  subroutine write_global_integer(this, label, integer_value)

    class(netcdf_writer_type), intent(in) :: this
    character(len=*), intent(in) :: label
    integer, intent(in) :: integer_value

    call nc_check(nf90_put_att(this%ncFileID, NF90_GLOBAL, label, integer_value))

  end subroutine write_global_integer


  !------------------------------------------------------------------
  ! Write_global_real
  !
  ! Write a global real variable
  !------------------------------------------------------------------
  subroutine write_global_real(this, label, real_value)

    class(netcdf_writer_type), intent(in) :: this
    character(len=*), intent(in) :: label
    real(r8kind), intent(in) :: real_value

    call nc_check(nf90_put_att(this%ncFileID, NF90_GLOBAL, label, real_value))

  end subroutine write_global_real


  !------------------------------------------------------------------
  ! Close
  !
  ! Close the file
  !------------------------------------------------------------------
  subroutine close(this)

    class(netcdf_writer_type), intent(in) :: this

    call nc_check(nf90_enddef(this%ncfileID))
    call nc_check(nf90_sync(this%ncFileID))
    call nc_check(nf90_close(this%ncFileID))

  end subroutine close


  !------------------------------------------------------------------
  ! nc_check
  !
  ! Checks return status from a NetCDF API call.  If an error was
  ! returned, print the message and abort the program.
  !------------------------------------------------------------------
  subroutine nc_check(istatus)

    integer, intent (in) :: istatus

    character(len=512) :: error_msg

    ! If no error, nothing to do here.  we are done.
    if (istatus == nf90_noerr) return

    ! Otherwise, print the error and stop
    error_msg = nf90_strerror(istatus)
    print *,error_msg
    stop

  end subroutine nc_check


end module NetCDF_Writer
