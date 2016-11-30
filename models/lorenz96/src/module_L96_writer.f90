module L96_Writer

  use kind,      only  : r8kind
  use L96_Model, only  : l96_model_type
  use L96_Config, only : l96_config_type
!  use Writer, only : writer_type

  implicit none

  private

  public :: l96_writer_type

  type :: l96_writer_type
      private
      character(len=7) :: io_format
  contains
      final :: destructor
      procedure :: write
      procedure :: ascii_write
      procedure :: netcdf_write
!      procedure :: write_config
!      procedure :: write_state
!      procedure, private :: netcdf_write_state
!      procedure, private :: ascii_write_state
  end type l96_writer_type

  interface l96_writer_type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized l96_writer_type object
  !------------------------------------------------------------------
  type(l96_writer_type) function constructor(io_format)

    character(len=*) :: io_format

    select case (io_format)
      case('NETCDF')
        constructor%io_format = io_format
      case('ASCII')
        constructor%io_format = io_format
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "', io_format, '" is not supported!'
        stop
    end select

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a l96_writer_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(l96_writer_type), intent(inout) :: this

    ! No pointers in l96_writer_type object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! write
  !------------------------------------------------------------------
  integer function write(this, model)

    class(l96_writer_type), intent(in) :: this
    class(l96_model_type),  intent(in) :: model

    integer :: ierr          ! return value of function

    select case (this%io_format)
      case('NETCDF')
        ierr = this%netcdf_write(model)
      case('ASCII')
        ierr = this%ascii_write(model)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "', this%io_format, '" is not supported!'
        stop
    end select


    write = ierr

  end function write


  !------------------------------------------------------------------
  ! ascii_write_model_state
  !------------------------------------------------------------------
  integer function ascii_write(this, model)

    class(l96_writer_type), intent(in) :: this
    class(l96_model_type),  intent(in) :: model

    integer :: ierr          ! return value of function

    type(l96_config_type) :: config
    real(r8kind), allocatable :: state(:)

    character(len=128)    :: filename    ! name of output file
    integer :: fileunit
    integer :: i
    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=19) :: timestr

    ! Get the model configuration
    config = model%get_config()

    ! Get the model state
    allocate(state(config%get_nx()))
    state = model%get_state()

    ! Construct name of output file
    write(filename,'(A,I0.7,A)') 'lorenz96out_', model%get_step(), '.csv'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted')

    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)

    ! Write global data
    write(fileunit,'(3A)') 'creation_date', ',', timestr
    write(fileunit,'(3A)') 'model', ',', 'Lorenz_96'
    write(fileunit,'(2A,F12.7)') 'model_forcing', ',', config%get_forcing()
    write(fileunit,'(2A,F12.7)') 'model_delta_t', ',', config%get_time_step()
    write(fileunit,'(2A,F15.7)') 'model_t', ',', model%get_clock()
    write(fileunit,'(2A,I)') 'model_step', ',', model%get_step()
    write(fileunit,'(2A,I)') 'StateDim', ',', config%get_nx()

    ! Write record separator
    write(fileunit,*)
    write(fileunit,*)

    ! Write the coordinate, location, and state fields
    write(fileunit,'(5A)') 'Coordinates',',','Location',',','State'
    do i=1, config%get_nx()
!      write(fileunit,'(I,2(A,F12.7))') i,',',this%location(i),',',this%state(i)
      write(fileunit,'(I,2(A,F12.7))') i,',',i,',',state(i)
!      write(fileunit,'(I,A,F12.7)') i,',',state(i)
    end do

    ! Close the file
    close(fileunit)

    ascii_write = ierr

  end function ascii_write


 !------------------------------------------------------------------
  ! netcdf_write
  !
  ! Writes model state to NetCDF file
  !
  ! Typical sequence for adding new dimensions,variables,attributes:
  ! NF90_OPEN             ! open existing netCDF dataset
  !    NF90_redef         ! put into define mode
  !    NF90_def_dim       ! define additional dimensions (if any)
  !    NF90_def_var       ! define variables: from name, type, and dims
  !    NF90_put_att       ! assign attribute values
  ! NF90_ENDDEF           ! end definitions: leave define mode
  !    NF90_put_var       ! provide values for variable
  ! NF90_CLOSE            ! close: save updated netCDF dataset
  !------------------------------------------------------------------
  integer function netcdf_write(this, model)

    use netcdf

    class(l96_writer_type), intent(in) :: this
    class(l96_model_type),  intent(in) :: model

    integer :: ierr          ! return value of function

    type(l96_config_type) :: config
    real(r8kind), allocatable :: state(:)

    ! General netCDF variables
    integer :: ncFileID      ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: StateVarDimID, CoordinatesVarID, LocationVarID, StateVarID

    ! local variables
    integer               :: i           ! loop index variable
    character(len=128)    :: filename    ! name of output file
    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=19) :: timestr

    ! assume normal termination
    ierr = 0 

    ! Get the model configuration
    config = model%get_config()

    ! Get the model state
    allocate(state(config%get_nx()))
    state = model%get_state()

    ! Construct name of output file
    write(filename,'(A,I0.7,A)') 'lorenz96out_', model%get_step(), '.nc'

    ! Open new file, overwriting previous contents
    call nc_check(nf90_create(trim(filename), NF90_CLOBBER, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Write Global Attributes 
    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",timestr))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Lorenz_96"))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_forcing", config%get_forcing() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_delta_t", config%get_time_step() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_t", model%get_clock() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_step", model%get_step() ))

    ! Define the model size
    call nc_check(nf90_def_dim(ncid=ncFileID, name="StateDim", &
                               len=config%get_nx(), dimid = StateVarDimID))

    ! Define the state vector coordinates
    call nc_check(nf90_def_var(ncid=ncFileID,name="Coordinates", xtype=nf90_int, &
                  dimids=StateVarDimID, varid=CoordinatesVarID))
    call nc_check(nf90_put_att(ncFileID, CoordinatesVarID, "long_name", "Model State Coordinates"))
    call nc_check(nf90_put_att(ncFileID, CoordinatesVarID, "units",     "Indexical"))
    call nc_check(nf90_put_att(ncFileID, CoordinatesVarID, "valid_range", (/ 1, config%get_nx() /)))

    ! Define the state vector locations
    call nc_check(NF90_def_var(ncFileID, name="Location", xtype=nf90_double, &
                  dimids = StateVarDimID, varid=LocationVarID))
    call nc_check(nf90_put_att(ncFileID, LocationVarID, "long_name", "Model State Location"))
    call nc_check(nf90_put_att(ncFileID, LocationVarID, "units", "Nondimensional"))
    call nc_check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8kind, 1.0_r8kind /)))

    ! Define the actual state vector
    call nc_check(nf90_def_var(ncid=ncFileID, name="State", xtype=nf90_double, &
               dimids=StateVarDimID, varid=StateVarID))
    call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "Model State"))
    call nc_check(nf90_put_att(ncFileID, StateVarID, "units", "Nondimensional"))

    ! Leave define mode so we can fill
    call nc_check(nf90_enddef(ncfileID))

    ! Fill the state coordinate variable
    call nc_check(nf90_put_var(ncFileID, CoordinatesVarID, (/ (i,i=1, config%get_nx()) /) ))

    ! Fill the location variable
!    call nc_check(nf90_put_var(ncFileID, LocationVarID, (/ (this%location(i),i=1, config%nx()) /) ))

    ! Fill the state variable
    call nc_check(nf90_put_var(ncFileID, StateVarID, (/ (state(i),i=1, config%get_nx()) /) ))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    netcdf_write = ierr

  end function netcdf_write


  !------------------------------------------------------------------
  ! nc_check
  ! 
  ! Checks return status from a NetCDF API call.  If an error was
  ! returned, print the message and abort the program.
  !------------------------------------------------------------------
  subroutine nc_check(istatus)

    use netcdf

    integer, intent (in)                   :: istatus
  
    character(len=512) :: error_msg
  
    ! if no error, nothing to do here.  we are done.
    if( istatus == nf90_noerr) return

    error_msg = nf90_strerror(istatus)
  
    print *,error_msg
    stop  

  end subroutine nc_check


end module L96_Writer
