module QG_Writer

  use kind,          only : r8kind
  use QG_Model,      only : qg_model_type
  use QG_Config,     only : qg_config_type
  use gptl
!  use NetCDF_Writer, only : netcdf_writer_type

  implicit none

  private

  public :: qg_writer_type

  type :: qg_writer_type
      private
      character(len=7)  :: io_format
      character(len=16) :: format_extension
!      type(netcdf_writer_type) :: io_writer
  contains
      final :: destructor
      procedure :: write
!      procedure :: generic_write
      procedure, private :: ascii_write
      procedure, private :: netcdf_write
  end type qg_writer_type

  interface qg_writer_type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized qg_writer_type object
  !------------------------------------------------------------------
  type(qg_writer_type) function constructor(io_format)

    character(len=*), intent(in) :: io_format

    constructor%io_format = io_format
    select case (io_format)
      case('NETCDF')
!        constructor%io_writer = netcdf_writer_type()
        constructor%format_extension = '.nc'
      case('ASCII')
!        constructor%io_writer = ascii_writer_type()
        constructor%format_extension = '.csv'
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "', io_format, '" is not supported!'
        stop
    end select

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a qg_writer_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(qg_writer_type), intent(inout) :: this

    ! No pointers in qg_writer_type object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! write
  !------------------------------------------------------------------
  subroutine write(this, model, filename)

    class(qg_writer_type), intent(inout) :: this
    class(qg_model_type),  intent(   in) :: model
    character(len=*),      intent(   in) :: filename


    select case (this%io_format)
      case('NETCDF')
        call this%netcdf_write(model, filename)
      case('ASCII')
        call this%ascii_write(model, filename)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "', this%io_format, '" is not supported!'
        stop
    end select


  end subroutine write


  !------------------------------------------------------------------
  ! generic_write
  !------------------------------------------------------------------
!  subroutine generic_write(this, model)
!
!    class(qg_writer_type), intent(inout) :: this
!    class(qg_model_type),  intent(   in) :: model
!
!    type(QG_config_type)  :: config      ! model configuration
!    character(len=128)    :: filename    ! name of output file
!    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
!    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
!    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
!    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
!    character(len=19)     :: timestr     ! date/time string
!
!    ! Get the model configuration
!    config = model%get_config()
!
!    ! Construct name of output file
!    write(filename,'(A,I0.7)') 'qgout_', model%get_step()
!
!    ! Create the output file
!    call this%io_writer%create(filename)
!
!    ! Construct a string containing the current date and time
!    call DATE_AND_TIME(crdate, crtime, crzone, values)
!    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
!          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)
!
!    ! Write the global header data
!    call this%io_writer%write_global_var('creation_date', timestr)
!    call this%io_writer%write_global_var('model', 'QG')
!    call this%io_writer%write_global_var('model_forcing', config%get_forcing())
!    call this%io_writer%write_global_var('model_time_step', config%get_time_step())
!    call this%io_writer%write_global_var('model_clock', model%get_clock())
!    call this%io_writer%write_global_var('model_step', model%get_step())
!
!    ! Close the output file
!    call this%io_writer%close()
!
!  end subroutine generic_write


  !------------------------------------------------------------------
  ! ascii_write
  !------------------------------------------------------------------
  subroutine ascii_write(this, model, filename)

    class(qg_writer_type), intent(in) :: this
    class(qg_model_type),  intent(in) :: model
    character(len=*),      intent(in) :: filename

    integer :: ierr ! return value of function

    type(QG_config_type) :: config
    real(r8kind), allocatable :: location(:)
    real(r8kind), allocatable :: state(:)

    integer :: fileunit
    integer :: i
    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=19)     :: timestr

    ! Get the model configuration
    config = model%get_config()

    ! Get the model location
!    allocate(location(config%get_nx()))
!    location = model%get_location()

    ! Get the model state
!    allocate(state(config%get_nx()))
!    state = model%get_state()

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename) // trim(this%format_extension), form='formatted')

    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)

    ! Write global data
    write(fileunit,'(3A)') 'creation_date', ',', timestr
    write(fileunit,'(3A)') 'model', ',', 'QG'
!    write(fileunit,'(2A,F12.7)') 'model_forcing', ',', config%get_forcing()
!    write(fileunit,'(2A,F12.7)') 'model_delta_t', ',', config%get_time_step()
!    write(fileunit,'(2A,F15.7)') 'model_t', ',', model%get_clock()
    write(fileunit,'(2A,I)') 'model_step', ',', model%get_step()
!    write(fileunit,'(2A,I)') 'StateDim', ',', config%get_nx()

    ! Write record separator
    write(fileunit,*)
    write(fileunit,*)

    ! Write the coordinate, location, and state fields
    write(fileunit,'(5A)') 'Coordinates',',','Location',',','State'
!    do i=1, config%get_nx()
!      write(fileunit,'(I,2(A,F12.7))') i,',',location(i),',',state(i)
!    end do

    ! Close the file
    close(fileunit)

  end subroutine ascii_write


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
  subroutine netcdf_write(this, model, filename)

    use netcdf

    class(qg_writer_type), intent(in) :: this
    class(qg_model_type),  intent(in) :: model
    character(len=*),      intent(in) :: filename

    type(QG_config_type) :: config
    integer :: nlat, nlon, nsh2, nvl

    ! Output grid fields
    real(r8kind), allocatable :: lat(:)      ! Grid Latitude
    real(r8kind), allocatable :: lon(:)      ! Grid Longitude
    real(r8kind), allocatable :: lvl(:)      ! Grid Level
    real(r8kind), allocatable :: geopg(:,:,:)! Geopotential on the grid
    real(r8kind), allocatable :: psig(:,:,:) ! Grid values of dimensional streamfunction at the three levels
    real(r8kind), allocatable :: forg(:,:,:) ! Grid values of dimensional forcing at the three levels
    real(r8kind), allocatable :: qgpv(:,:,:) ! Grid values of dimensional pv at the three levels
    real(r8kind), allocatable :: ug(:,:,:)   ! Grid values of zonal velocity at the three levels in m/s
    real(r8kind), allocatable :: vg(:,:,:)   ! Grid values of meridional velocity at the three levels in m/s

    ! General netCDF variables
    integer :: ncFileID      ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: LatDimID, LonDimID, Nsh2DimID, NvlDimID
    integer :: LatVarID, LonVarID, LvlVarID, GeopgVarID, PsiVarID, PsigVarID, ForVarID, ForgVarID, QgpvVarID, UgVarID, VgVarID

    ! local variables
    integer               :: i, j, k     ! loop index variable
    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=19)     :: timestr

    ! Get the model configuration
    config = model%get_config()

    ! Get model dimensions
    nlat = model%get_nlat()
    nlon = model%get_nlon()
    nsh2 = model%get_nsh2()
    nvl = model%get_nvl()

    ! Allocate space for output fields
    allocate(lat(nlat))
    allocate(lon(nlon))
    allocate(lvl(nvl))
    allocate(geopg(nlat,nlon,nvl))
    allocate(psig(nlat,nlon,nvl))
    allocate(forg(nlat,nlon,nvl))
    allocate(qgpv(nlat,nlon,nvl))
    allocate(ug(nlat,nlon,nvl))
    allocate(vg(nlat,nlon,nvl))

    ! Get model streamfunction and derived fields on gaussian grid
    call model%gridfields(lat, lon, lvl, geopg, psig, forg, qgpv, ug, vg)

    ! Open new file, overwriting previous contents
    call nc_check(nf90_create(trim(filename) // trim(this%format_extension), NF90_CLOBBER, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Write Global Attributes 
    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",timestr))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "QG"))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "resolution", config%get_resolution() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "time_step", config%get_time_step() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "obsfile", config%get_obsfile() ))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "inf", config%get_inf() ))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "obsf", config%get_obsf() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "tdis", config%get_tdis() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "addisl", config%get_addisl() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "addish", config%get_addish() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "trel", config%get_trel() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "tdif", config%get_tdif() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "idif", config%get_idif() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "h0", config%get_h0() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "rrdef1", config%get_rrdef1() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "rrdef2", config%get_rrdef2() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "clock", model%get_clock() ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "step", model%get_step() ))

    ! Define the lat/lon dimensions
    call nc_check(nf90_def_dim(ncid=ncFileID, name="lat", len=nlat, dimid = LatDimID))
    call nc_check(nf90_def_dim(ncid=ncFileID, name="lon", len=nlon, dimid = LonDimID))
    call nc_check(nf90_def_dim(ncid=ncFileID, name="level", len=nvl, dimid = NvlDimID))

    ! Define the spectral dimensions
    call nc_check(nf90_def_dim(ncid=ncFileID, name="nsh2", len=nsh2, dimid = Nsh2DimID))

    ! Define the lat field
    call nc_check(nf90_def_var(ncid=ncFileID,name="lat", xtype=nf90_double, &
                  dimids=(/LatDimID/), varid=LatVarID))
    call nc_check(nf90_put_att(ncFileID, LatVarID, "long_name", "latitude"))
    call nc_check(nf90_put_att(ncFileID, LatVarID, "units",     "degrees_north"))
    call nc_check(nf90_put_att(ncFileID, LatVarID, "valid_range", (/ -90.0_r8kind, 90.0_r8kind /)))

    ! Define the lon field
    call nc_check(NF90_def_var(ncFileID, name="lon", xtype=nf90_double, &
                  dimids=(/LonDimID/), varid=LonVarID))
    call nc_check(nf90_put_att(ncFileID, LonVarID, "long_name", "longitude"))
    call nc_check(nf90_put_att(ncFileID, LonVarID, "units", "degrees_east"))
    call nc_check(nf90_put_att(ncFileID, LonVarID, "valid_range", (/ 0.0_r8kind, 360.0_r8kind /)))

    ! Define the level field
    call nc_check(NF90_def_var(ncFileID, name="level", xtype=nf90_double, &
                  dimids=(/NvlDimID/), varid=LvlVarID))
    call nc_check(nf90_put_att(ncFileID, LvlVarID, "long_name", "pressure_level"))
    call nc_check(nf90_put_att(ncFileID, LvlVarID, "units", "millibar"))
    call nc_check(nf90_put_att(ncFileID, LvlVarID, "valid_range", (/ 800, 200, 200 /)))

    ! Define the geopg variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="Zg", xtype=nf90_double, &
                  dimids=(/LonDimID, LatDimID, NvlDimID/), varid=GeopgVarID))
    call nc_check(nf90_put_att(ncFileID, GeopgVarID, "long_name", "Geopotential Height"))
    call nc_check(nf90_put_att(ncFileID, GeopgVarID, "units", "m2 / s2"))

    ! Define the spectral streamfunction variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="Psi", xtype=nf90_double, &
                  dimids=(/Nsh2DimID, NvlDimID/), varid=PsiVarID))
    call nc_check(nf90_put_att(ncFileID, PsiVarID, "long_name", "Spectral Streamfunction"))
    call nc_check(nf90_put_att(ncFileID, PsiVarID, "units", "m2 / s"))

    ! Define the gaussian grid streamfunction variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="Psig", xtype=nf90_double, &
                  dimids=(/LonDimID, LatDimID, NvlDimID/), varid=PsigVarID))
    call nc_check(nf90_put_att(ncFileID, PsigVarID, "long_name", "Gaussian Streamfunction"))
    call nc_check(nf90_put_att(ncFileID, PsigVarID, "units", "m2 / s"))

    ! Define the spectral forcing variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="For", xtype=nf90_double, &
                  dimids=(/Nsh2DimID, NvlDimID/), varid=ForVarID))
    call nc_check(nf90_put_att(ncFileID, ForVarID, "long_name", "Spectral Forcing"))
    call nc_check(nf90_put_att(ncFileID, ForVarID, "units", "Nondimensional"))

    ! Define the gaussian grid forcing variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="Forg", xtype=nf90_double, &
                  dimids=(/LonDimID, LatDimID, NvlDimID/), varid=ForgVarID))
    call nc_check(nf90_put_att(ncFileID, ForgVarID, "long_name", "Gaussian Forcing"))
    call nc_check(nf90_put_att(ncFileID, ForgVarID, "units", "Nondimensional"))

    ! Define the qgpv variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="Q", xtype=nf90_double, &
                  dimids=(/LonDimID, LatDimID, NvlDimID/), varid=QgpvVarID))
    call nc_check(nf90_put_att(ncFileID, QgpvVarID, "long_name", "Potential Vorticity"))
    call nc_check(nf90_put_att(ncFileID, QgpvVarID, "units", "m2 / s"))

    ! Define the ug variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="U", xtype=nf90_double, &
                  dimids=(/LonDimID, LatDimID, NvlDimID/), varid=UgVarID))
    call nc_check(nf90_put_att(ncFileID, UgVarID, "long_name", "Zonal Wind"))
    call nc_check(nf90_put_att(ncFileID, UgVarID, "units", "m / s"))

    ! Define the vg variable
    call nc_check(nf90_def_var(ncid=ncFileID, name="V", xtype=nf90_double, &
                  dimids=(/LonDimID, LatDimID, NvlDimID/), varid=VgVarID))
    call nc_check(nf90_put_att(ncFileID, VgVarID, "long_name", "Meridional Wind"))
    call nc_check(nf90_put_att(ncFileID, VgVarID, "units", "m / s"))

    ! Leave define mode so we can fill
    call nc_check(nf90_enddef(ncfileID))

    ! Fill the lat variable
    call nc_check(nf90_put_var(ncFileID, LatVarID, lat))

    ! Fill the lon variable
    call nc_check(nf90_put_var(ncFileID, LonVarID, lon))

    ! Fill the level variable
    call nc_check(nf90_put_var(ncFileID, LvlVarID, lvl))

    ! Fill the geopotential height variable
    call nc_check(nf90_put_var(ncFileID, GeopgVarID, geopg, count=(/nlon, nlat, nvl/), map=(/nlat, 1, nlat*nlon/) ))

    ! Fill the spectral streamfunction variable
    call nc_check(nf90_put_var(ncFileID, PsiVarID, model%get_psi()))

    ! Fill the gaussian grid streamfunction variable
    call nc_check(nf90_put_var(ncFileID, PsigVarID, psig, count=(/nlon, nlat, nvl/), map=(/nlat, 1, nlat*nlon/) ))

    ! Fill the spectral forcing variable
    call nc_check(nf90_put_var(ncFileID, ForVarID, model%get_for()))

    ! Fill the gaussian grid forcing variable
    call nc_check(nf90_put_var(ncFileID, ForgVarID, forg, count=(/nlon, nlat, nvl/), map=(/nlat, 1, nlat*nlon/) ))

    ! Fill the potential vorticticy
    call nc_check(nf90_put_var(ncFileID, QgpvVarID, qgpv, count=(/nlon, nlat, nvl/), map=(/nlat, 1, nlat*nlon/) ))

    ! Fill the wind variables
    call nc_check(nf90_put_var(ncFileID, UgVarID, ug, count=(/nlon, nlat, nvl/), map=(/nlat, 1, nlat*nlon/) ))
    call nc_check(nf90_put_var(ncFileID, VgVarID, vg, count=(/nlon, nlat, nvl/), map=(/nlat, 1, nlat*nlon/) ))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

  end subroutine netcdf_write


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


end module QG_Writer
