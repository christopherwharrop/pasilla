module QG_Reader

  use kind,      only  : r8kind
  use QG_Model, only  : qg_model_type
  use QG_Config, only : qg_config_type
!  use Reader, only : reader_type

  implicit none

  private

  public :: qg_reader_type

  type :: qg_reader_type
      private
      character(len=7) :: io_format
      character(len=16) :: format_extension
  contains
      final :: destructor
      procedure :: read
      procedure, private :: ascii_read
      procedure, private :: netcdf_read
  end type qg_reader_type

  interface qg_reader_type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized qg_reader_type object
  !------------------------------------------------------------------
  type(qg_reader_type) function constructor(io_format)

    character(len=*) :: io_format

    select case (io_format)
      case('NETCDF')
        constructor%io_format = io_format
        constructor%format_extension = '.nc'
      case('ASCII')
        constructor%io_format = io_format
        constructor%format_extension = '.csv'
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "', io_format, '" is not supported!'
        stop
    end select

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a qg_reader_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(qg_reader_type), intent(inout) :: this

    ! No pointers in qg_reader_type object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! read
  !------------------------------------------------------------------
  subroutine read(this, model, filename)

    class(qg_reader_type),    intent(   in) :: this
    class(qg_model_type),     intent(inout) :: model
    character(len=*),          intent(   in) :: filename

    select case (this%io_format)
      case('NETCDF')
        call this%netcdf_read(model, filename)
      case('ASCII')
        call this%ascii_read(model, filename)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "', this%io_format, '" is not supported!'
        stop
    end select

  end subroutine read


  !------------------------------------------------------------------
  ! ascii_read
  !------------------------------------------------------------------
  subroutine ascii_read(this, model, filename)

    class(qg_reader_type),    intent(   in) :: this
    type(qg_model_type),      intent(inout) :: model
    character(len=*),          intent(   in) :: filename

    
    integer :: fileunit
    character(len=80) :: line
    character(len=16) :: linefmt
    character(len=64) :: attr_name
    integer :: position
    integer :: ignore
    integer :: i
    integer               :: nx
    real(r8kind)          :: forcing
    real(r8kind)          :: time_step
    type(qg_config_type) :: config
    integer               :: step
    real(r8kind)          :: clock
    real(r8kind), allocatable :: state(:)
    real(r8kind), allocatable :: location(:)

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename) // trim(this%format_extension), form='formatted', status='old')

    ! Read global attributes
    read(fileunit, '(A)') line
    position = index(line, ',')
    do while (position /= 0)

      ! Read global attribute name
      write(linefmt, '(A,I0,A)') '(A', position - 1, ')'
      read(line, linefmt) attr_name

      ! Read in global attribute value
      select case (attr_name)
        case('model_forcing')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) forcing
        case('model_delta_t')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) time_step
        case('model_t')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) clock
        case('model_step')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',I)'
          read(line, linefmt) step
        case('StateDim')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',I)'
          read(line, linefmt) nx
        case DEFAULT
          ! Ignore gloval settings we don't need
          read(line,*)
      end select

      ! Get the next line and position of the comma
      read(fileunit, '(A)') line
      position = index(line, ',')

    end do

    ! Read record separator
    read(fileunit, '(A)') line

    ! Read field header
    read(fileunit, '(A)') line

    ! Allocate space for the model state and location
    allocate(state(nx))
    allocate(location(nx))

    ! Read the coordinate, location, and state fields
    do i=1, nx
      read(fileunit, *) ignore, location(i), state(i)
    end do

    ! Close the file
    close(fileunit)

    ! Create the model configuration
!    config = qg_config_type(nx, time_step, forcing)

    ! Instantiate a model from the configuration
!    model = qg_model_type(config, state=state, step=step)

  end subroutine ascii_read


  !------------------------------------------------------------------
  ! netcdf_read
  !
  ! Reads model state to NetCDF file
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
  subroutine netcdf_read(this, model, filename)

    use netcdf

    class(qg_reader_type), intent(   in) :: this
    type(qg_model_type),   intent(inout) :: model
    character(len=*),      intent(   in) :: filename

    ! Config variables
    type(QG_config_type) :: config    ! Model configuration
    integer              :: resolution ! Model resolution
    integer              :: time_step  ! Model time step (in seconds)
    character(len=32)    :: obsfile    ! Name of observation file
    logical              :: inf        ! If .true. then artificial PV forcing read from file
    logical              :: obsf       ! If .true. PV forcing is calculated from observations in routine artiforc
    real(r8kind)         :: tdis       ! Ekman dissipation timescale in days at lower level
    real(r8kind)         :: addisl     ! Parameter used in the computation of the dissipation timescale at the lower level over land
    real(r8kind)         :: addish     ! Parameter used in the computation of the dissipation timescale at the lower level as a function of topography
    real(r8kind)         :: trel       ! Relaxation time scale in days of the temperature
    real(r8kind)         :: tdif       ! Dissipation timescale of scale-selective diffusion in days for wavenumber nm
    integer              :: idif       ! Determines scale-selectivity of hyperviscosity; power of laplace operator
    real(r8kind)         :: h0         ! Scale factor for the topographically induced upward motion at the lower level
    real(r8kind)         :: rrdef1     ! Rossby radius of deformation of 200-500 thickness
    real(r8kind)         :: rrdef2     ! Rossby radius of deformation of 500-800 thickness

    ! Model variables
    integer :: step                        ! Current model step
    integer :: nsh2, nvl                   ! Model dimensions in spectral space
    real(r8kind), allocatable :: psi(:,:)  ! Stream function at the nvl levela
    real(r8kind), allocatable :: for(:,:)  ! Constant potential vorticity forcing at the nvl levels

    ! General netCDF variables
    integer :: ncFileID      ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: Nsh2DimID, NvlDimID
    integer :: PsiVarID, ForVarID

    ! Open file for read only
    call nc_check(nf90_open(trim(filename) // trim(this%format_extension), NF90_NOWRITE, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Read Global Attributes
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "resolution", resolution ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "time_step", time_step ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "obsfile", obsfile ))
!    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "inf", inf ))
!    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "obsf", obsf ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "tdis", tdis ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "addisl", addisl ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "addish", addish ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "trel", trel ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "tdif", tdif ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "idif", idif ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "h0", h0 ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "rrdef1", rrdef1 ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "rrdef2", rrdef2 ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "step", step ))

    ! Read the spectral model dimensions
    call nc_check(nf90_inq_dimid(ncFileID, "level", NvlDimID))
    call nc_check(nf90_inquire_dimension(ncFileID, NvlDimID, len=nvl))
    call nc_check(nf90_inq_dimid(ncFileID, "nsh2", Nsh2DimID))
    call nc_check(nf90_inquire_dimension(ncFileID, Nsh2DimID, len=nsh2))

    ! Allocate space for the model streamfunction and forcing
    allocate(psi(nsh2, nvl))
    allocate(for(nsh2, nvl))

    ! Get the streamfunction ID
    call nc_check(nf90_inq_varid(ncFileID, "Psi", PsiVarID))

    ! Get the forcing ID
    call nc_check(nf90_inq_varid(ncFileID, "For", ForVarID))

    ! Get the streamfunction variable
    call nc_check(nf90_get_var(ncFileID, PsiVarID, psi))

    ! Get the forcing variable
    call nc_check(nf90_get_var(ncFileID, ForVarID, for))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    ! Create the model configuration
    config = qg_config_type(resolution, time_step, obsfile, .false., .false., tdis, addisl, addish, trel, tdif, idif, h0, rrdef1, rrdef2)

    ! Instantiate a model from the configuration
    model = qg_model_type(config, state=psi, for=for, step=step)

  end subroutine netcdf_read


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


end module QG_Reader
