module sine

  use kind, only : r8kind
  use module_constants, only : PI

  implicit none

  private

  public :: sine_type, sine_TL_type, sine_ADJ_type

  type sine_type
      private
      integer      :: size
      real(r8kind) :: amplitude
      real(r8kind) :: bias
      real(r8kind) :: frequency
      real(r8kind) :: phase
      real(r8kind) :: delta_t
      real(r8kind) :: t
      integer, public      :: step
      real(r8kind), allocatable, public :: state(:)
      real(r8kind), allocatable :: location(:)
  contains
      final              :: destructor
      procedure          :: adv_nsteps
      procedure          :: interpolate
      procedure          :: read_model_state
      procedure, private :: netcdf_read_model_state
      procedure, private :: ascii_read_model_state
      procedure          :: write_model_state
      procedure, private :: netcdf_write_model_state
      procedure, private :: ascii_write_model_state
  end type sine_type


  type, extends(sine_type) :: sine_TL_type
      private
      real(r8kind), allocatable, public :: trajectory(:)
  contains
      final              :: destructor_TL
      procedure          :: adv_nsteps => adv_nsteps_d
  end type sine_TL_type


  type, extends(sine_type) :: sine_ADJ_type
      private
      real(r8kind), allocatable, public :: trajectory(:)
  contains
      final              :: destructor_ADJ
      procedure          :: adv_nsteps => adv_nsteps_b
  end type sine_ADJ_type


  interface sine_type
    procedure constructor_parm
    procedure constructor_file
  end interface

  interface sine_TL_type
    procedure constructor_TL_parm
    procedure constructor_TL_file
  end interface

  interface sine_ADJ_type
    procedure constructor_ADJ_parm
    procedure constructor_ADJ_file
  end interface

contains

  !------------------------------------------------------------------
  ! constructor_parm
  !
  ! Returns an initialized sine object
  !------------------------------------------------------------------
  type(sine_type) function constructor_parm(size, amplitude, bias, frequency, phase, delta_t)

    integer, intent(in)      :: size
    real(r8kind), intent(in) :: amplitude
    real(r8kind), intent(in) :: bias
    real(r8kind), intent(in) :: frequency
    real(r8kind), intent(in) :: phase
    real(r8kind), intent(in) :: delta_t

    integer :: j
    real(r8kind) :: t

    ! Initialize model parameters    
    constructor_parm%size = size
    constructor_parm%amplitude = amplitude
    constructor_parm%bias = bias
    constructor_parm%frequency = frequency
    constructor_parm%phase = phase
    constructor_parm%delta_t = delta_t
    constructor_parm%t = 0
    constructor_parm%step = 0

    ! Allocate model variables
    allocate(constructor_parm%state(size))
    allocate(constructor_parm%location(size))

    ! Initialize model variables
    t = 0.0
    do j = 1, size
      constructor_parm%state(j) = bias + amplitude * sin(((frequency * (t + phase) + float(j)) / size) * PI * 4)
    end do

    ! Localize the domain
    do j = 1, size
      constructor_parm%location(j) = (j - 1.0_r8kind) / size
    end do

  end function


  !------------------------------------------------------------------
  ! constructor_file
  !
  ! Returns an initialized sine object
  !------------------------------------------------------------------
  type(sine_type) function constructor_file(read_step, format)

    integer, intent(in)      :: read_step
    character(*), intent(in) :: format

    integer :: ierr

    select case (format)
      case('NETCDF')
        ! Read the header
        ierr = netcdf_read_model_header(read_step, constructor_file%size, constructor_file%amplitude, &
                                 & constructor_file%bias, constructor_file%frequency, constructor_file%phase, &
                                 & constructor_file%delta_t, constructor_file%t, constructor_file%step)
      case('ASCII')
        ! Read the header
        ierr = ascii_read_model_header(read_step, constructor_file%size, constructor_file%amplitude, &
                                 & constructor_file%bias, constructor_file%frequency, constructor_file%phase, &
                                 & constructor_file%delta_t, constructor_file%t, constructor_file%step)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "',format,'" is not supported!'
        stop
    end select

    ! Allocate space for the data
    allocate(constructor_file%state(constructor_file%size))
    allocate(constructor_file%location(constructor_file%size))

    ! Read the data
    select case (format)
      case('NETCDF')
        ierr = netcdf_read_model_data(read_step, constructor_file%size, constructor_file%location, constructor_file%state)
      case('ASCII')
        ierr = ascii_read_model_data(read_step, constructor_file%size, constructor_file%location, constructor_file%state)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "',format,'" is not supported!'
        stop
    end select

  end function


  !------------------------------------------------------------------
  ! constructor_TL_parm
  !
  ! Returns an initialized sine_TL object
  !------------------------------------------------------------------
  type(sine_TL_type) function constructor_TL_parm(size, amplitude, bias, frequency, phase, delta_t)

    integer, intent(in)      :: size
    real(r8kind), intent(in) :: amplitude
    real(r8kind), intent(in) :: bias
    real(r8kind), intent(in) :: frequency
    real(r8kind), intent(in) :: phase
    real(r8kind), intent(in) :: delta_t

    ! Call constructor for superclass
    constructor_TL_parm%sine_type = sine_type(size, amplitude, bias, frequency, phase, delta_t)

    ! Allocate model variables
    allocate(constructor_TL_parm%trajectory(size))

    ! Initialize model variables
    constructor_TL_parm%trajectory = constructor_TL_parm%state

  end function


  !------------------------------------------------------------------
  ! constructor_TL_file
  !
  ! Returns an initialized sine_TL object
  !------------------------------------------------------------------
  type(sine_TL_type) function constructor_TL_file(read_step, format)

    integer, intent(in)      :: read_step
    character(*), intent(in) :: format

    ! Call constructor for superclass
    constructor_TL_file%sine_type = sine_type(read_step, format)

    ! Allocate model variables
    allocate(constructor_TL_file%trajectory(constructor_TL_file%size))

    ! Initialize model variables
    constructor_TL_file%trajectory = constructor_TL_file%state

  end function


  !------------------------------------------------------------------
  ! constructor_ADJ_parm
  !
  ! Returns an initialized sine_ADJ object
  !------------------------------------------------------------------
  type(sine_ADJ_type) function constructor_ADJ_parm(size, amplitude, bias, frequency, phase, delta_t)

    integer, intent(in)      :: size
    real(r8kind), intent(in) :: amplitude
    real(r8kind), intent(in) :: bias
    real(r8kind), intent(in) :: frequency
    real(r8kind), intent(in) :: phase
    real(r8kind), intent(in) :: delta_t

    ! Call constructor for superclass
    constructor_ADJ_parm%sine_type = sine_type(size, amplitude, bias, frequency, phase, delta_t)

    ! Allocate model variables
    allocate(constructor_ADJ_parm%trajectory(size))

    ! Initialize model variables
    constructor_ADJ_parm%trajectory = constructor_ADJ_parm%state

  end function


  !------------------------------------------------------------------
  ! constructor_ADJ_file
  !
  ! Returns an initialized sine_ADJ object
  !------------------------------------------------------------------
  type(sine_ADJ_type) function constructor_ADJ_file(read_step, format)

    integer, intent(in)      :: read_step
    character(*), intent(in) :: format

    ! Call constructor for superclass
    constructor_ADJ_file%sine_type = sine_type(read_step, format)

    ! Allocate model variables
    allocate(constructor_ADJ_file%trajectory(constructor_ADJ_file%size))

    ! Initialize model variables
    constructor_ADJ_file%trajectory = constructor_ADJ_file%state

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a sine object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(sine_type), intent(inout) :: this

    ! No pointers in sine object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! destructor_TL
  !
  ! Deallocates pointers used by a sine_TL object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_TL(this)

    type(sine_TL_type), intent(inout) :: this

    ! No pointers in sine_TL object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a sine_ADJ object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_ADJ(this)

    type(sine_ADJ_type), intent(inout) :: this

    ! No pointers in sine_ADJ object so we do nothing

  end subroutine

  !------------------------------------------------------------------
  ! adv_nsteps
  !
  !------------------------------------------------------------------
  subroutine adv_nsteps(this, nsteps)

    class(sine_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    integer :: step, j
    real(r8kind) :: phaset, phasedt
    
    do step = 1, nsteps

      ! Increment time step
      this%t = this%t + this%delta_t
      this%step = this%step + 1

      do j = 1, this%size
        this%state(j) = this%bias + this%amplitude * sin((this%frequency * (this%t + this%phase) + float(j)) / this%size * PI * 4)
      end do

    end do

  end subroutine adv_nsteps


  !------------------------------------------------------------------
  ! adv_nsteps_d
  !
  !------------------------------------------------------------------
  subroutine adv_nsteps_d(this, nsteps)

    class(sine_TL_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    integer :: step, j

    do step = 1, nsteps

      ! Increment time step
      this%t = this%t + this%delta_t
      this%step = this%step + 1

      do j=1, this%size
        this%state(j) = this%trajectory(j) + ((this%amplitude * PI * 4 * this%frequency) / this%size) * cos((this%frequency * (this%t + this%phase) + float(j)) / this%size * PI * 4)
        this%trajectory(j) = this%bias + this%amplitude * sin((this%frequency * (this%t + this%phase) + float(j)) / this%size * PI * 4)
      end do

    end do

  end subroutine adv_nsteps_d


  !------------------------------------------------------------------
  ! adv_nsteps_b
  !
  !------------------------------------------------------------------
  subroutine adv_nsteps_b(this, nsteps)

    class(sine_ADJ_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    integer :: step, j

    do step = 1, nsteps

      ! Increment time step
      this%t = this%t - this%delta_t
      this%step = this%step - 1

      do j = 1, this%size
        this%state(j) = this%trajectory(j) - ((this%amplitude * PI * 4 * this%frequency) / this%size) * cos((this%frequency * (this%t + this%phase) + float(j)) / this%size * PI * 4)
        this%trajectory(j) = this%bias + this%amplitude * sin((this%frequency * (this%t + this%phase) + float(j)) / this%size * PI * 4)
      end do

    end do

  end subroutine adv_nsteps_b


  !------------------------------------------------------------------  
  ! Interpolates from state vector x to the location. 
  !------------------------------------------------------------------  
  subroutine interpolate(this, location, state_val)

    class(sine_type), intent(in) :: this
    real(r8kind), intent(in)    :: location
    real(r8kind), intent(out)   :: state_val

    integer :: lower_index, upper_index, i
    real(r8kind) :: lctn, lctnfrac

    ! Scale the location to the size of the domain
    lctn = this%size * location

    ! Compute grid indices bounding the location
    lower_index = int(lctn) + 1
    upper_index = lower_index + 1
    if(lower_index > this%size) lower_index = lower_index - this%size
    if(upper_index > this%size) upper_index = upper_index - this%size

    ! Interpolate model value at the location
    lctnfrac = lctn - int(lctn)
    state_val = (1.0_r8kind - lctnfrac) * this%state(lower_index) + lctnfrac * this%state(upper_index)
 
  end subroutine interpolate


  !------------------------------------------------------------------
  ! write_model_state
  !------------------------------------------------------------------
  integer function write_model_state(this, format)

    class(sine_type), intent(in) :: this
    character(*), intent(in)    :: format

    integer :: ierr          ! return value of function

    select case (format)
      case('NETCDF')
        ierr = this%netcdf_write_model_state()
      case('ASCII')
        ierr = this%ascii_write_model_state()
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "',format,'" is not supported!'
        stop
    end select

    write_model_state = ierr

  end function write_model_state


  !------------------------------------------------------------------
  ! netcdf_write_model_state
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
  integer function netcdf_write_model_state(this)

    use netcdf

    class(sine_type), intent(in) :: this

    integer :: ierr          ! return value of function

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

    ! Construct name of output file
    write(filename,'(A,I0.7,A)') 'sineout_', this%step, '.nc'

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
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Sine Wave"))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_amplitude", this%amplitude ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_bias", this%bias ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_frequency", this%frequency ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_phase", this%phase ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_delta_t", this%delta_t ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_t", this%t ))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_step", this%step ))

    ! Define the model size
    call nc_check(nf90_def_dim(ncid=ncFileID, name="StateDim", &
                               len=this%size, dimid = StateVarDimID))

    ! Define the state vector coordinates
    call nc_check(nf90_def_var(ncid=ncFileID,name="Coordinates", xtype=nf90_int, &
                  dimids=StateVarDimID, varid=CoordinatesVarID))
    call nc_check(nf90_put_att(ncFileID, CoordinatesVarID, "long_name", "Model State Coordinates"))
    call nc_check(nf90_put_att(ncFileID, CoordinatesVarID, "units",     "Indexical"))
    call nc_check(nf90_put_att(ncFileID, CoordinatesVarID, "valid_range", (/ 1, this%size /)))

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
    call nc_check(nf90_put_var(ncFileID, CoordinatesVarID, (/ (i,i=1,this%size) /) ))

    ! Fill the location variable
    call nc_check(nf90_put_var(ncFileID, LocationVarID, (/ (this%location(i),i=1,this%size) /) ))

    ! Fill the state variable
    call nc_check(nf90_put_var(ncFileID, StateVarID, (/ (this%state(i),i=1,this%size) /) ))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    netcdf_write_model_state = ierr

  end function netcdf_write_model_state


  !------------------------------------------------------------------
  ! ascii_write_model_state
  !------------------------------------------------------------------
  integer function ascii_write_model_state(this)

    class(sine_type), intent(in) :: this

    integer :: ierr          ! return value of function

    character(len=128)    :: filename    ! name of output file
    integer :: fileunit
    integer :: i
    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=19) :: timestr

    ! Construct name of output file
    write(filename,'(A,I0.7,A)') 'sineout_', this%step, '.csv'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted')

    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)

    ! Write global data
    write(fileunit,'(3A)') 'creation_date', ',', timestr
    write(fileunit,'(3A)') 'model', ',', 'Sine Wave'
    write(fileunit,'(2A,F12.7)') 'model_amplitude', ',', this%amplitude
    write(fileunit,'(2A,F12.7)') 'model_bias', ',', this%bias
    write(fileunit,'(2A,F12.7)') 'model_frequency', ',', this%frequency
    write(fileunit,'(2A,F12.7)') 'model_phase', ',', this%phase
    write(fileunit,'(2A,F12.7)') 'model_delta_t', ',', this%delta_t
    write(fileunit,'(2A,F15.7)') 'model_t', ',', this%t
    write(fileunit,'(2A,I)') 'model_step', ',', this%step
    write(fileunit,'(2A,I)') 'StateDim', ',', this%size

    ! Write record separator
    write(fileunit,*)
    write(fileunit,*)

    ! Write the coordinate, location, and state fields
    write(fileunit,'(5A)') 'Coordinates',',','Location',',','State'
    do i=1, this%size
      write(fileunit,'(I,2(A,F12.7))') i,',',this%location(i),',',this%state(i)
    end do

    ! Close the file
    close(fileunit)

    ascii_write_model_state = ierr

  end function ascii_write_model_state


  !------------------------------------------------------------------
  ! read_model_state
  !------------------------------------------------------------------
  integer function read_model_state(this, read_step, format)

    class(sine_type), intent(inout) :: this
    integer, intent(in)            :: read_step
    character(*), intent(in)       :: format

    integer :: ierr          ! return value of function

    select case (format)
      case('NETCDF') 
        ierr = this%netcdf_read_model_state(read_step)
      case('ASCII')
        ierr = this%ascii_read_model_state(read_step)
      case DEFAULT
        write(*,'(A,A,A)') 'ERROR: IO Format "',format,'" is not supported!'
        stop
    end select

    ! Initialize trajectories for TL and ADJ
    select type(this)
      class is (sine_TL_type)
        this%trajectory = this%state
      class is (sine_ADJ_type)
        this%trajectory = this%state
      class default
        ! Do nothing
    end select

    read_model_state = ierr

  end function read_model_state


  !------------------------------------------------------------------
  ! netcdf_read_model_header
  !------------------------------------------------------------------
  integer function netcdf_read_model_header(read_step, size, amplitude, bias, frequency, phase, delta_t, t, step)

    use netcdf

    integer, intent(in)       :: read_step ! Read in data for this time step
    integer, intent(out)      :: size
    real(r8kind), intent(out) :: amplitude
    real(r8kind), intent(out) :: bias
    real(r8kind), intent(out) :: frequency
    real(r8kind), intent(out) :: phase
    real(r8kind), intent(out) :: delta_t
    real(r8kind), intent(out) :: t
    integer, intent(out)      :: step

    ! General netCDF variables
    integer :: ncFileID  ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: StateVarDimID, CoordinatesVarID, LocationVarID, StateVarID

    ! local variables
    character(len=128) :: filename
    integer :: ierr

    ! assume normal termination
    ierr = 0

    ! Calculate name of file based on time step requested
    write(filename,'(A,I0.7,A)') 'sineout_', read_step, '.nc'

    ! Open file for read only
    call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Read Global Attributes
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_amplitude", amplitude ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_bias", bias ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_frequency", frequency ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_phase", phase ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_delta_t", delta_t ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_t", t ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_step", step ))

    ! Read the model size
    call nc_check(nf90_inq_dimid(ncFileID, "StateDim", StateVarDimID))
    call nc_check(nf90_inquire_dimension(ncFileID, StateVarDimID, len=size))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    netcdf_read_model_header = ierr

  end function netcdf_read_model_header


  !------------------------------------------------------------------
  ! netcdf_read_model_data
  !------------------------------------------------------------------
  integer function netcdf_read_model_data(read_step, size, location, state)

    use netcdf

    integer, intent(in) :: read_step ! Read in data for this time step
    integer, intent(in)         :: size
    real(r8kind), intent(inout) :: location(:)
    real(r8kind), intent(inout) :: state(:)

    integer :: ierr  ! return value of function

    ! General netCDF variables
    integer :: ncFileID  ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: StateVarDimID, CoordinatesVarID, LocationVarID, StateVarID

    ! local variables
    character(len=128) :: filename

    ! assume normal termination
    ierr = 0

    ! Calculate name of file based on time step requested
    write(filename,'(A,I0.7,A)') 'sineout_', read_step, '.nc'

    ! Open file for read only
    call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Get the state vector location ID
     call nc_check(nf90_inq_varid(ncFileID, "Location", LocationVarID))

    ! Get the actual state vector ID
    call nc_check(nf90_inq_varid(ncFileID, "State", StateVarID))

    ! Get the location variable
    call nc_check(nf90_get_var(ncFileID, LocationVarID, location))

    ! Get the state variable
    call nc_check(nf90_get_var(ncFileID, StateVarID, state))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    netcdf_read_model_data = ierr

  end function netcdf_read_model_data


  !------------------------------------------------------------------
  ! netcdf_read_model_state
  !------------------------------------------------------------------
  integer function netcdf_read_model_state(this,read_step)

    use netcdf

    class(sine_type), intent(inout) :: this
    integer, intent(in) :: read_step ! Read in data for this time step

    integer :: ierr  ! return value of function

    ! local variables
    integer      :: size
    real(r8kind) :: amplitude
    real(r8kind) :: bias
    real(r8kind) :: frequency
    real(r8kind) :: phase
    real(r8kind) :: delta_t
    real(r8kind) :: t
    integer      :: step
    character(len=128) :: filename

    ! assume normal termination
    ierr = 0 

    ! Calculate name of file based on time step requested
    write(filename,'(A,I0.7,A)') 'sineout_', read_step, '.nc'

    ! Read the model header
    ierr = netcdf_read_model_header(read_step, size, amplitude, bias, frequency, phase, delta_t, t, step)

    ! Validate the input
    if (amplitude /= this%amplitude) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file amplitude =',amplitude,', expecting ',this%amplitude
      stop
    end if
    if (bias /= this%bias) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file bias =',bias,', expecting ',this%bias
      stop
    end if
    if (frequency /= this%frequency) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file frequency =',frequency,', expecting ',this%frequency
      stop
    end if
    if (phase /= this%phase) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file phase =',phase,', expecting ',this%phase
      stop
    end if
    if (delta_t /= this%delta_t) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file delta_t =',delta_t,', expecting ',this%delta_t
      stop
    end if

    ! Set the time and the step
    this%t = t
    this%step = step

    ! Read the model data
    ierr = netcdf_read_model_data(read_step, this%size, this%location, this%state)

    netcdf_read_model_state = ierr

  end function netcdf_read_model_state


  !------------------------------------------------------------------
  ! ascii_read_model_header
  !------------------------------------------------------------------
  integer function ascii_read_model_header(read_step, size, amplitude, bias, frequency, phase, delta_t, t, step)

    integer, intent(in)       :: read_step ! Read in data for this time step
    integer, intent(out)      :: size
    real(r8kind), intent(out) :: amplitude
    real(r8kind), intent(out) :: bias
    real(r8kind), intent(out) :: frequency
    real(r8kind), intent(out) :: phase
    real(r8kind), intent(out) :: delta_t
    real(r8kind), intent(out) :: t
    integer, intent(out)      :: step

    integer :: ierr                  ! return value of function

    character(len=128) :: filename   ! name of output file
    integer :: fileunit
    character(len=80) :: line
    character(len=16) :: linefmt
    character(len=64) :: attr_name
    integer :: position
    integer :: ignore
    integer :: i

    ! assume normal termination
    ierr = 0

    ! Construct name of input file
    write(filename, '(A,I0.7,A)') 'sineout_', read_step, '.csv'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted', status='old')

    ! Read global attributes
    read(fileunit, '(A)') line
    position = index(line, ',')
    do while (position /= 0)

      ! Read global attribute name
      write(linefmt, '(A,I0,A)') '(A', position - 1, ')'
      read(line, linefmt) attr_name

      ! Read in global attribute value
      select case (attr_name)
        case('model_amplitude')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) amplitude
        case('model_bias')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) bias
        case('model_frequency')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) frequency
        case('model_phase')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) phase
        case('model_delta_t')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) delta_t
        case('model_t')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) t
        case('model_step')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',I)'
          read(line, linefmt) step
        case('StateDim')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',I)'
          read(line, linefmt) size
        case DEFAULT
          ! Ignore global settings we don't need
          read(line,*)
      end select

      ! Get the next line and position of the comma
      read(fileunit, '(A)') line
      position = index(line, ',')

    end do

    ! Close the file
    close(fileunit)

    ascii_read_model_header = ierr

  end function ascii_read_model_header


  !------------------------------------------------------------------
  ! ascii_read_model_data
  !------------------------------------------------------------------
  integer function ascii_read_model_data(read_step, size, location, state)

    integer, intent(in)         :: read_step ! Read in data for this time step
    integer, intent(in)         :: size
    real(r8kind), intent(inout) :: location(:)
    real(r8kind), intent(inout) :: state(:)

    integer :: ierr                  ! return value of function

    character(len=128) :: filename   ! name of output file
    integer :: fileunit
    character(len=80) :: line
    character(len=64) :: attr_name
    integer :: position
    integer :: ignore
    integer :: i

    ! Construct name of input file
    write(filename, '(A,I0.7,A)') 'sineout_', read_step, '.csv'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted', status='old')

    ! Read past the global attributes
    read(fileunit, '(A)') line
    position = index(line, ',')
    do while (position /= 0)

      ! Get the next line and position of the comma
      read(fileunit, '(A)') line
      position = index(line, ',')

    end do

    ! Read record separator
    read(fileunit, '(A)') line

    ! Read field header
    read(fileunit, '(A)') line

    ! Read the coordinate, location, and state fields
    do i=1, size
      read(fileunit, *) ignore, location(i), state(i)
    end do

    ! Close the file
    close(fileunit)

    ascii_read_model_data = ierr

  end function ascii_read_model_data


  !------------------------------------------------------------------
  ! ascii_read_model_state
  !------------------------------------------------------------------
  integer function ascii_read_model_state(this,read_step)

    class(sine_type), intent(inout) :: this
    integer, intent(in) :: read_step ! Read in data for this time step

    integer :: ierr                  ! return value of function

    character(len=128) :: filename   ! name of output file
    integer      :: size
    real(r8kind) :: amplitude
    real(r8kind) :: bias
    real(r8kind) :: frequency
    real(r8kind) :: phase
    real(r8kind) :: delta_t
    real(r8kind) :: t
    integer      :: step

    ! assume normal termination
    ierr = 0

    ! Construct name of input file
    write(filename, '(A,I0.7,A)') 'sineout_', read_step, '.csv'

    ! Read the header data
    ierr = ascii_read_model_header(read_step, size, amplitude, bias, frequency, phase, delta_t, t, step)

    ! Validate the input
    if (amplitude /= this%amplitude) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file amplitude =',amplitude,', expecting ',this%amplitude
      stop
    end if
    if (bias /= this%bias) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file bias =',bias,', expecting ',this%bias
      stop
    end if
    if (frequency /= this%frequency) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file frequency =',frequency,', expecting ',this%frequency
      stop
    end if
    if (phase /= this%phase) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file phase =',phase,', expecting ',this%phase
      stop
    end if
    if (delta_t /= this%delta_t) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file delta_t =',delta_t,', expecting ',this%delta_t
      stop
    end if
    if (size /= this%size) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,I,A,I)') '       Input file size =',size,', expecting ',this%size
      stop
    end if

    ! Set the time and the step
    this%t = t
    this%step = step

    ! Read the model data
    ierr = ascii_read_model_data(read_step, this%size, this%location, this%state)

    ascii_read_model_state = ierr

  end function ascii_read_model_state

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


end module sine
