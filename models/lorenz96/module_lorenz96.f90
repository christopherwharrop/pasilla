module module_lorenz96

  use module_kind, only : r8kind

  implicit none

  private

  public :: lorenz96

  type lorenz96
!    private
      integer  :: size
      real(r8kind) :: forcing
      real(r8kind) :: delta_t
      real(r8kind) :: t
      integer :: step
      real(r8kind), allocatable :: x(:)
      real(r8kind), allocatable :: location(:)
  contains
      final              :: destructor
      procedure, private :: comp_dt
      procedure          :: adv_nsteps
      procedure          :: interpolate
      procedure          :: nc_read_model_state
      procedure          :: nc_write_model_state
  end type lorenz96

  interface lorenz96
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized lorenz96 object
  !------------------------------------------------------------------
  type(lorenz96) function constructor(size, forcing, delta_t)

    integer      :: size
    real(r8kind) :: forcing
    real(r8kind) :: delta_t

    integer :: j

    ! Initialize model parameters    
    constructor%size = size
    constructor%forcing = forcing
    constructor%delta_t = delta_t
    constructor%t = 0
    constructor%step = 0

    ! Allocate model variables
    allocate(constructor%x(size))
    allocate(constructor%location(size))

    ! Initialize model variables
    constructor%x = forcing
    constructor%x(1) = 1.001_r8kind * forcing

    ! Localize the domain
    do j = 1, size
      constructor%location(j) = (j - 1.0_r8kind) / size
    end do

  end function

  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a lorenz96 object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(lorenz96), intent(inout) :: this

    ! No pointers in lorenz96 object so we do nothing

  end subroutine

  !------------------------------------------------------------------
  ! comp_dt
  !
  ! Private routine to compute the time tendency of a lorenz96 object 
  ! given a state, x, and return it in dt.
  !------------------------------------------------------------------
  subroutine comp_dt(this, x, dt)

    class(lorenz96), intent(in) :: this
    real(r8kind), intent( in)   :: x(:)
    real(r8kind), intent(out)   :: dt(:)

    integer :: j, jp1, jm1, jm2

    do j = 1, this%size
       jp1 = j + 1
       if(jp1 > this%size) jp1 = 1
       jm2 = j - 2
       if(jm2 < 1) jm2 = this%size + jm2
       jm1 = j - 1
       if(jm1 < 1) jm1 = this%size

       dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + this%forcing
    end do

  end subroutine comp_dt


  !------------------------------------------------------------------
  ! adv_nsteps
  !
  ! Does n time step advances for lorenz 96 model
  ! using four-step rk time step
  !------------------------------------------------------------------
  subroutine adv_nsteps(this,steps)

    class(lorenz96), intent(inout) :: this
    integer, intent(in) :: steps

    real(r8kind), dimension(size(this%x)) :: x1, x2, x3, x4, dx, inter
    integer :: step
    
    do step = 1, steps

      call this%comp_dt(this%x, dx)   !  Compute the first intermediate step
      x1    = this%delta_t * dx
      inter = this%x + x1 / 2.0_r8kind

      call this%comp_dt(inter, dx)    !  Compute the second intermediate step
      x2    = this%delta_t * dx
      inter = this%x + x2 / 2.0_r8kind

      call this%comp_dt(inter, dx)    !  Compute the third intermediate step
      x3    = this%delta_t * dx
      inter = this%x + x3

      call this%comp_dt(inter, dx)    !  Compute fourth intermediate step
      x4 = this%delta_t * dx

      !  Compute new value for x
      this%x = this%x + x1/6.0_r8kind + x2/3.0_r8kind + x3/3.0_r8kind + x4/6.0_r8kind

      ! Increment time step
      this%t = this%t + this%delta_t
      this%step = this%step + 1

    end do

  end subroutine adv_nsteps


  !------------------------------------------------------------------  
  ! Interpolates from state vector x to the location. 
  !------------------------------------------------------------------  
  subroutine interpolate(this, location, obs_val)

    class(lorenz96), intent(in) :: this
    real(r8kind), intent(in)    :: location
    real(r8kind), intent(out)   :: obs_val

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
    obs_val = (1.0_r8kind - lctnfrac) * this%x(lower_index) + lctnfrac * this%x(upper_index)
 
  end subroutine interpolate


  !------------------------------------------------------------------
  ! nc_write_model_state
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
  integer function nc_write_model_state(this)

    use netcdf

    class(lorenz96), intent(in) :: this

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
    character(len=NF90_MAX_NAME) :: timestr

    ! assume normal termination
    ierr = 0 

    ! Construct name of output file
    write(filename,'(A,I0.6,A)') 'lorenz96out_', this%step, '.nc'

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
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_forcing", this%forcing ))
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
    call nc_check(nf90_put_var(ncFileID, StateVarID, (/ (this%x(i),i=1,this%size) /) ))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    nc_write_model_state = ierr

  end function nc_write_model_state


  !------------------------------------------------------------------
  ! nc_read_model_state
  !------------------------------------------------------------------
  integer function nc_read_model_state(this,filename)

    use netcdf

    class(lorenz96), intent(inout) :: this
    character(*), intent(in)       :: filename  ! name of input file

    integer :: ierr  ! return value of function

    ! General netCDF variables
    integer :: ncFileID  ! netCDF file identifier
    integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
    integer :: StateVarDimID, CoordinatesVarID, LocationVarID, StateVarID

    ! local variables
    integer      :: i  ! loop index variable
    integer      :: size
    real(r8kind) :: forcing
    real(r8kind) :: delta_t

    ! assume normal termination
    ierr = 0 

    ! Open file for read only
    call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Read Global Attributes 
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_forcing", forcing ))
    if (forcing /= this%forcing) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file forcing =',forcing,', expecting ',this%forcing
      stop
    end if
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_delta_t", delta_t ))
    if (delta_t /= this%delta_t) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file delta_t =',delta_t,', expecting ',this%delta_t
      stop
    end if
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_t", this%t ))
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_step", this%step ))

    ! Read the model size
    call nc_check(nf90_inq_dimid(ncFileID, "StateDim", StateVarDimID))
    call nc_check(nf90_inquire_dimension(ncFileID, StateVarDimID, len=size))
    if (size /= this%size) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,I,A,I)') '       Input file size =',size,', expecting ',this%size
      stop
    end if

    ! Get the state vector location ID
     call nc_check(nf90_inq_varid(ncFileID, "Location", LocationVarID))

    ! Get the actual state vector ID
    call nc_check(nf90_inq_varid(ncFileID, "State", StateVarID))

    ! Get the location variable
    call nc_check(nf90_get_var(ncFileID, LocationVarID, this%location))

    ! Get the state variable
    call nc_check(nf90_get_var(ncFileID, StateVarID, this%x))

    ! Flush buffers
    call nc_check(nf90_sync(ncFileID))

    ! Close the NetCDF file
    call nc_check(nf90_close(ncFileID))

    nc_read_model_state = ierr

  end function nc_read_model_state



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




end module module_lorenz96
