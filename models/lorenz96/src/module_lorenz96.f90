module lorenz96

  use kind, only : r8kind
  use model, only : model_type

  implicit none

  private

  public :: lorenz96_type, lorenz96_TL_type, lorenz96_ADJ_type

  type, extends(model_type) :: lorenz96_type
      private
      integer, public      :: size
      real(r8kind) :: forcing
      real(r8kind) :: delta_t
      real(r8kind) :: t
      integer, public      :: step
      real(r8kind), allocatable, public :: state(:)
      real(r8kind), allocatable, public :: location(:)
  contains
      final              :: destructor
      procedure, private :: comp_dt
      procedure          :: adv_nsteps
      procedure          :: interpolate
      procedure          :: read_model_state
      procedure, private :: netcdf_read_model_state
      procedure, private :: ascii_read_model_state
      procedure          :: write_model_state
      procedure, private :: netcdf_write_model_state
      procedure, private :: ascii_write_model_state
  end type lorenz96_type


  type, extends(lorenz96_type) :: lorenz96_TL_type
      private
      real(r8kind), allocatable, public :: trajectory(:)
  contains
      final              :: destructor_TL
      procedure          :: adv_nsteps => adv_nsteps_d
      procedure, private :: comp_dt_d
  end type lorenz96_TL_type


  type, extends(lorenz96_type) :: lorenz96_ADJ_type
      private
      real(r8kind), allocatable, public :: trajectory(:)
  contains
      final              :: destructor_ADJ
      procedure          :: adv_nsteps => adv_nsteps_b
      procedure, private :: comp_dt_b
  end type lorenz96_ADJ_type


  interface lorenz96_type
    procedure constructor_parm
    procedure constructor_file
  end interface

  interface lorenz96_TL_type
     procedure constructor_TL_parm
     procedure constructor_TL_file
  end interface

  interface lorenz96_ADJ_type
    procedure constructor_ADJ_parm
    procedure constructor_ADJ_file
  end interface

contains

  !------------------------------------------------------------------
  ! constructor_parm
  !
  ! Returns an initialized lorenz96 object
  !------------------------------------------------------------------
  type(lorenz96_type) function constructor_parm(size, forcing, delta_t)

    integer, intent(in)      :: size
    real(r8kind), intent(in) :: forcing
    real(r8kind), intent(in) :: delta_t

    integer :: j

    ! Initialize model parameters    
    constructor_parm%size = size
    constructor_parm%forcing = forcing
    constructor_parm%delta_t = delta_t
    constructor_parm%t = 0
    constructor_parm%step = 0

    ! Allocate model variables
    allocate(constructor_parm%state(size))
    allocate(constructor_parm%location(size))

    ! Initialize model variables
    constructor_parm%state = forcing
    constructor_parm%state(1) = 1.001_r8kind * forcing

    ! Localize the domain
    do j = 1, size
      constructor_parm%location(j) = (j - 1.0_r8kind) / size
    end do

  end function


  !------------------------------------------------------------------
  ! constructor_file
  !
  ! Returns an initialized lorenz96 object
  !------------------------------------------------------------------
  type(lorenz96_type) function constructor_file(read_step, format)

    integer, intent(in)      :: read_step
    character(*), intent(in) :: format

    integer :: ierr

    select case (format)
      case('NETCDF')
        ! Read the header
        ierr = netcdf_read_model_header(read_step, constructor_file%size, constructor_file%forcing, &
                                 & constructor_file%delta_t, constructor_file%t, constructor_file%step)
      case('ASCII')
        ! Read the header
        ierr = ascii_read_model_header(read_step, constructor_file%size, constructor_file%forcing, &
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
  ! Returns an initialized lorenz96_TL object
  !------------------------------------------------------------------
  type(lorenz96_TL_type) function constructor_TL_parm(size, forcing, delta_t)

    integer, intent(in)      :: size
    real(r8kind), intent(in) :: forcing
    real(r8kind), intent(in) :: delta_t

    ! Call constructor for superclass
    constructor_TL_parm%lorenz96_type = lorenz96_type(size, forcing, delta_t)

    ! Allocate model variables
    allocate(constructor_TL_parm%trajectory(size))

    ! Initialize model variables
    constructor_TL_parm%trajectory = constructor_TL_parm%state

  end function


  !------------------------------------------------------------------
  ! constructor_TL_file
  !
  ! Returns an initialized lorenz96 object
  !------------------------------------------------------------------
  type(lorenz96_TL_type) function constructor_TL_file(read_step, format)

    integer, intent(in)      :: read_step
    character(*), intent(in) :: format

    ! Call constructor for superclass
    constructor_TL_file%lorenz96_type = lorenz96_type(read_step, format)

    ! Allocate model variables
    allocate(constructor_TL_file%trajectory(constructor_TL_file%size))

    ! Initialize model variables
    constructor_TL_file%trajectory = constructor_TL_file%state


  end function


  !------------------------------------------------------------------
  ! constructor_ADJ_parm
  !
  ! Returns an initialized lorenz96_ADJ object
  !------------------------------------------------------------------
  type(lorenz96_ADJ_type) function constructor_ADJ_parm(size, forcing, delta_t)

    integer, intent(in)      :: size
    real(r8kind), intent(in) :: forcing
    real(r8kind), intent(in) :: delta_t

    ! Call constructor for superclass
    constructor_ADJ_parm%lorenz96_type = lorenz96_type(size, forcing, delta_t)

    ! Allocate model variables
    allocate(constructor_ADJ_parm%trajectory(size))

    ! Initialize model variables
    constructor_ADJ_parm%trajectory = constructor_ADJ_parm%state

  end function


  !------------------------------------------------------------------
  ! constructor_ADJ_file
  !
  ! Returns an initialized lorenz96_ADJ object
  !------------------------------------------------------------------
  type(lorenz96_ADJ_type) function constructor_ADJ_file(read_step, format)

    integer, intent(in)      :: read_step
    character(*), intent(in) :: format

    ! Call constructor for superclass
    constructor_ADJ_file%lorenz96_type = lorenz96_type(read_step, format)

    ! Allocate model variables
    allocate(constructor_ADJ_file%trajectory(constructor_ADJ_file%size))

    ! Initialize model variables
    constructor_ADJ_file%trajectory = constructor_ADJ_file%state

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a lorenz96 object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(lorenz96_type), intent(inout) :: this

    ! No pointers in lorenz96 object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! destructor_TL
  !
  ! Deallocates pointers used by a lorenz96_TL object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_TL(this)

    type(lorenz96_TL_type), intent(inout) :: this

    ! No pointers in lorenz96_TL object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a lorenz96_ADJ object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_ADJ(this)

    type(lorenz96_ADJ_type), intent(inout) :: this

    ! No pointers in lorenz96_ADJ object so we do nothing

  end subroutine


  !------------------------------------------------------------------
  ! comp_dt
  !
  ! Private routine to compute the time tendency of a lorenz96_type object 
  ! given a state, x, and return it in dt.
  !------------------------------------------------------------------
  subroutine comp_dt(this, x, dt)

    class(lorenz96_type), intent(in) :: this
    real(r8kind), intent( in) ::  x(:)
    real(r8kind), intent(out) :: dt(:)

    integer :: j, jp1, jm1, jm2

    ! compute dt(1)
    dt(1) = (x(2) - x(this%size - 1)) * x(this%size) - x(1) + this%forcing

    ! compute dt(2)
    dt(2) = (x(3) - x(this%size)) * x(1) - x(2) + this%forcing 

    ! compute dt(3) thru dt(size -1)
    do j = 3, this%size - 1
       jp1 = j + 1
       jm2 = j - 2
       jm1 = j - 1
       dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + this%forcing
    end do

    ! compute dt(size)
    dt(this%size) = (x(1) - x(this%size - 2)) * x(this%size - 1) - x(this%size) + this%forcing

  end subroutine comp_dt


  !------------------------------------------------------------------
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.11 (r6148) - 16 Aug 2016 14:18
  !
  !  Differentiation of comp_dt in forward (tangent) mode:
  !   variations   of useful results: dt
  !   with respect to varying inputs: dt x
  !------------------------------------------------------------------
  ! comp_dt_d
  !
  ! Private routine to compute the time tendency of a lorenz96_TL_type object
  ! given a state, x, and return it in dt.
  !------------------------------------------------------------------
  subroutine comp_dt_d(this, x, xd, dt, dtd)

    class(lorenz96_TL_type), intent(in) :: this
    real(r8kind), intent( in) ::   x(this%size)
    real(r8kind), intent( in) ::  xd(this%size)
    real(r8kind), intent(out) ::  dt(this%size)
    real(r8kind), intent(out) :: dtd(this%size)

    integer :: j

    ! compute dtd(1)
    dtd(1) = (xd(2) - xd(this%size - 1)) * x(this%size) + (x(2) - x(this%size - 1)) * xd(this%size) - xd(1)
    dt(1) =  ( x(2) -  x(this%size - 1)) * x(this%size) - x(1) + this%forcing

    ! compute dtd(2)
    dtd(2) = (xd(3) - xd(this%size)) * x(1) + (x(3) - x(this%size)) * xd(1) - xd(2)
    dt(2) =  (x(3) -   x(this%size)) * x(1) - x(2) + this%forcing

    ! compute dtd(3) thru dtd(size -1)
    do j = 3, this%size - 1
       dtd(j) = (xd(j + 1) - xd(j - 2)) * x(j - 1) + (x(j + 1) - x(j - 2)) * xd(j - 1) - xd(j)
        dt(j) = ( x(j + 1) -  x(j - 2)) * x(j - 1) - x(j) + this%forcing
    end do

    ! compute dtd(size)
    dtd(this%size) = (xd(1) - xd(this%size - 2)) * x(this%size - 1) + (x(1) - x(this%size - 2)) * xd(this%size - 1) - xd(this%size)
    dt(this%size) =   (x(1) -  x(this%size - 2)) * x(this%size - 1) - x(this%size) + this%forcing

  end subroutine comp_dt_d


  !------------------------------------------------------------------
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.11 (r6148) - 16 Aug 2016 14:18
  !
  !  Differentiation of comp_dt in reverse (adjoint) mode:
  !   gradient     of useful results: dt x
  !   with respect to varying inputs: dt x
  !------------------------------------------------------------------
  ! comp_dt
  !
  ! Private routine to compute the time tendency of a lorenz96_type object 
  ! given a state, x, and return it in dt.
  !------------------------------------------------------------------
  subroutine comp_dt_b(this, x, xb, dt, dtb)

    class(lorenz96_ADJ_type), intent(in) :: this
    real(r8kind), intent(   in) ::   x(this%size)
    real(r8kind), intent(inout) ::  xb(this%size)
    real(r8kind), intent(   in) ::  dt(this%size)
    real(r8kind), intent(inout) :: dtb(this%size)

    integer :: j
    real(r8kind) :: tempb

    xb(1) = xb(1) + x(this%size - 1) * dtb(this%size)
    xb(this%size - 2) = xb(this%size - 2) - x(this%size - 1) * dtb(this%size)
    xb(this%size - 1) = xb(this%size - 1) + (x(1) - x(this%size - 2)) * dtb(this%size)
    xb(this%size) = xb(this%size) - dtb(this%size)
    dtb(this%size) = 0.0

    do j = this%size - 1, 3, -1
      tempb = x(j-1) * dtb(j)
      xb(j + 1) = xb(j + 1) + tempb
      xb(j - 2) = xb(j - 2) - tempb
      xb(j - 1) = xb(j - 1) + (x(j + 1) - x(j - 2)) * dtb(j)
      xb(j) = xb(j) - dtb(j)
      dtb(j) = 0.0
    end do

    xb(3) = xb(3) + x(1) * dtb(2)
    xb(this%size) = xb(this%size) - x(1) * dtb(2)
    xb(1) = xb(1) + (x(3) - x(this%size)) * dtb(2)
    xb(2) = xb(2) - dtb(2)
    dtb(2) = 0.0

    xb(2) = xb(2) + x(this%size) * dtb(1)
    xb(this%size - 1) = xb(this%size - 1) - x(this%size) * dtb(1)
    xb(this%size) = xb(this%size) + (x(2) - x(this%size - 1)) * dtb(1)
    xb(1) = xb(1) - dtb(1)
    dtb(1) = 0.0

  end subroutine comp_dt_b

  !------------------------------------------------------------------
  ! adv_nsteps
  !
  ! Does n time step advances for lorenz 96 model
  ! using four-step rk time step
  !------------------------------------------------------------------
  subroutine adv_nsteps(this, nsteps)

    class(lorenz96_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    real(r8kind), dimension(this%size) :: x1, x2, x3, x4, dx, inter
    integer :: step
    
    do step = 1, nsteps

      call this%comp_dt(this%state, dx)   !  Compute the first intermediate step
      x1    = this%delta_t * dx
      inter = this%state + x1 / 2.0_r8kind

      call this%comp_dt(inter, dx)        !  Compute the second intermediate step
      x2    = this%delta_t * dx
      inter = this%state + x2 / 2.0_r8kind

      call this%comp_dt(inter, dx)        !  Compute the third intermediate step
      x3    = this%delta_t * dx
      inter = this%state + x3

      call this%comp_dt(inter, dx)        !  Compute fourth intermediate step
      x4 = this%delta_t * dx

      !  Compute new value for state
      this%state = this%state + x1 / 6.0_r8kind + x2 / 3.0_r8kind + x3 / 3.0_r8kind + x4 / 6.0_r8kind

      ! Increment time step
      this%t = this%t + this%delta_t
      this%step = this%step + 1

    end do

  end subroutine adv_nsteps


  !------------------------------------------------------------------
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.11 (r6148) - 16 Aug 2016 14:18
  !
  !  Differentiation of adv_nsteps in forward (tangent) mode:
  !   variations   of useful results: state
  !   with respect to varying inputs: state
  !   RW status of diff variables: state:in-out
  !------------------------------------------------------------------
  ! adv_nsteps_d
  !
  ! Does n time step advances for lorenz 96 model
  ! using four-step rk time step
  !------------------------------------------------------------------
  subroutine adv_nsteps_d(this, nsteps)

    class(lorenz96_TL_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    real(r8kind), dimension(this%size) :: x1, x2, x3, x4, dx, inter
    real(r8kind), dimension(this%size) :: x1d, x2d, x3d, x4d, dxd, interd
    real(r8kind), dimension(this%size,this%size) :: mprime

    integer :: step

    do step = 1, nsteps

! TAPENADE
!       dxd = 0.0
!       call this%comp_dt_d(this%trajectory, this%state, dx, dxd)  !  Compute the first intermediate step
!       x1d = this%delta_t * dxd
!       x1 = this%delta_t * dx
!       interd = this%state + x1d / 2.0
!       inter = this%trajectory + x1 / 2.0
!
!       call this%comp_dt_d(inter, interd, dx, dxd)            !  Compute the second intermediate step
!       x2d = this%delta_t * dxd
!       x2 = this%delta_t * dx
!       interd = this%state + x2d / 2.0
!       inter = this%trajectory + x2 / 2.0
!
!       call this%comp_dt_d(inter, interd, dx, dxd)            !  Compute the third intermediate step
!       x3d = this%delta_t * dxd
!       x3 = this%delta_t * dx
!       interd = this%state + x3d
!       inter = this%trajectory + x3
!
!       call this%comp_dt_d(inter, interd, dx, dxd)            !  Compute fourth intermediate step
!       x4d = this%delta_t * dxd
!       x4 = this%delta_t * dx
!
!       !  Compute new value for x
!       this%state = this%state + x1d / 6.0 + x2d / 3.0 + x3d / 3.0 + x4d / 6.0
!       this%trajectory = this%trajectory + x1 / 6.0 + x2 / 3.0 + x3 / 3.0 + x4 / 6.0
!

! Lidia's matrix formulation with correct trajectory
44    FORMAT(A8,6F10.6)
      call buildMprime(this%state, mprime)
      this%state = this%state + this%trajectory
      this%trajectory = this%trajectory + this%delta_t * matmul(mprime, this%trajectory)
      ! Increment time step
      this%t = this%t + this%delta_t
      this%step = this%step + 1

    end do

  end subroutine adv_nsteps_d


  !------------------------------------------------------------------
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.11 (r6148) - 16 Aug 2016 14:18
  !
  !  Differentiation of adv_nsteps in reverse (adjoint) mode:
  !   gradient     of useful results: state
  !   with respect to varying inputs: state
  !   RW status of diff variables: state:in-out
  !------------------------------------------------------------------
  ! adv_nsteps_b
  !
  ! Does n time step advances for lorenz 96 model
  ! using four-step rk time step
  !------------------------------------------------------------------
  subroutine adv_nsteps_b(this, nsteps)

    class(lorenz96_ADJ_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    real(r8kind), dimension(this%size) :: x1, x2, x3, x4, dx, inter
    real(r8kind), dimension(this%size) :: x1b, x2b, x3b, x4b, dxb, interb

    integer :: step
    real(r8kind), dimension(this%size,this%size) :: mprime

    do step = 1, nsteps

! TAPENADE code
!      call this%comp_dt(this%trajectory, dx)
!      x1 = this%delta_t * dx
!      inter = this%trajectory + x1 / 2.0
!
!      call this%comp_dt(inter, dx)
!      x2 = this%delta_t * dx
!      inter = this%trajectory + x2 / 2.0
!
!      call this%comp_dt(inter, dx)
!      x3 = this%delta_t * dx
!      inter = this%trajectory + x3
!
!      x1b = this%state / 6.0
!      x2b = this%state / 3.0
!      x4b = this%state / 6.0
!      dxb = this%delta_t * x4b
!      interb = 0.0
!
!      call this%comp_dt_b(inter, interb, dx, dxb)
!      x3b = interb + this%state / 3.0
!      this%state = this%state + interb
!      dxb = dxb + this%delta_t * x3b
!      inter = this%trajectory + x2 / 2.0
!      interb = 0.0
!
!      call this%comp_dt_b(inter, interb, dx, dxb)
!      this%state = this%state + interb
!      x2b = x2b + interb / 2.0
!      dxb = dxb + this%delta_t * x2b
!      inter = this%trajectory + x1 / 2.0
!      interb = 0.0
!
!      call this%comp_dt_b(inter, interb, dx, dxb)
!      this%state = this%state + interb
!      x1b = x1b + interb / 2.0
!      dxb = dxb + this%delta_t * x1b
!
!      call this%comp_dt_b(this%trajectory, this%state, dx, dxb)
!


! Lidia's MATRIX formulation with correct trajectory
44    FORMAT (A8,6F10.6)
      call buildMprime(this%state, mprime)
      x2 = -this%delta_t * matmul(mprime,this%trajectory)
      x1 = -this%delta_t * matmul(transpose(mprime), x2)
      this%trajectory = this%trajectory + x2 + x1  
      this%state = this%state + this%trajectory
      ! Increment time step
      this%t = this%t - this%delta_t
      this%step = this%step - 1

    end do

  end subroutine adv_nsteps_b


  !------------------------------------------------------------------
  ! buildMprime
  !
  ! Constructs matrix for computing lorenz96 derivatives.
  ! Used to advance the TL and ADJ models.
  !------------------------------------------------------------------
  subroutine buildMprime(x,mm)

    real(r8kind) :: x(:)
    real(r8kind) :: mm(:,:)

    real(r8kind), allocatable :: tmp(:)

    integer :: n,i

    n=size(x)
    allocate(tmp(n))
    mm(:,:)=0.0
    do i=1,n
      tmp = 0.0
      tmp(1) = -1.0
      tmp(2) = x(mmod(i-1, n))
      tmp(n-1) = -x(mmod(i-1,n))
      tmp(n) = x(mmod(i+1,n))-x(mmod(i-2,n))
      tmp =cshift(tmp,-1*(i-1))
      mm(i,:) = tmp
    end do

  end subroutine buildMprime


  !------------------------------------------------------------------
  ! mmod
  !
  ! Helper function used by buildMprime to compute modified mod function
  !------------------------------------------------------------------
  integer function mmod(x,n)

    integer :: x,n

    mmod=modulo(x,n)
    if (mmod==0) mmod=n

  end function mmod


  !------------------------------------------------------------------  
  ! Interpolates from state vector x to the location. 
  !------------------------------------------------------------------  
  subroutine interpolate(this, location, state_val)

    class(lorenz96_type), intent(in) :: this
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

    class(lorenz96_type), intent(in) :: this
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

    class(lorenz96_type), intent(in) :: this

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
    write(filename,'(A,I0.7,A)') 'lorenz96out_', this%step, '.nc'

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

    class(lorenz96_type), intent(in) :: this

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
    write(filename,'(A,I0.7,A)') 'lorenz96out_', this%step, '.csv'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted')

    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(timestr,'(i4,2(a,i2.2),1x,i2.2,2(a,i2.2))') &
          values(1), '/', values(2), '/', values(3), values(5), ':', values(6), ':', values(7)

    ! Write global data
    write(fileunit,'(3A)') 'creation_date', ',', timestr
    write(fileunit,'(3A)') 'model', ',', 'Lorenz_96'
    write(fileunit,'(2A,F12.7)') 'model_forcing', ',', this%forcing
    write(fileunit,'(2A,F12.7)') 'model_delta_t', ',', this%delta_t
    write(fileunit,'(2A,F15.7)') 'model_t', ',', this%t
    write(fileunit,'(2A,I)') 'model_step', ',', this%step
    write(fileunit,'(2A,I)') 'StateDim', ',', this%size

    ! Write record separator
    write(fileunit,*)
    write(fileunit,*)

    ! Write the coordinate, location, and state fields
    write(fileunit,'(5A)') 'Coordinates',',','Location',',','State'
    do i=1, size(this%state)
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

    class(lorenz96_type), intent(inout) :: this
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
      class is (lorenz96_TL_type)
        this%trajectory = this%state
      class is (lorenz96_ADJ_type)
        this%trajectory = this%state
      class default
        ! Do nothing
    end select

    read_model_state = ierr

  end function read_model_state


  !------------------------------------------------------------------
  ! netcdf_read_model_header
  !------------------------------------------------------------------
  integer function netcdf_read_model_header(read_step, size, forcing, delta_t, t, step)

    use netcdf

    integer, intent(in)       :: read_step ! Read in data for this time step
    integer, intent(out)      :: size
    real(r8kind), intent(out) :: forcing
    real(r8kind), intent(out) :: delta_t
    real(r8kind), intent(out) :: t
    integer, intent(out)      :: step

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
    write(filename,'(A,I0.7,A)') 'lorenz96out_', read_step, '.nc'

    ! Open file for read only
    call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncFileID))
    call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

    ! Read Global Attributes 
    call nc_check(nf90_get_att(ncFileID, NF90_GLOBAL, "model_forcing", forcing ))
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
    integer      :: i  ! loop index variable
    character(len=128) :: filename

    ! assume normal termination
    ierr = 0

    ! Calculate name of file based on time step requested
    write(filename,'(A,I0.7,A)') 'lorenz96out_', read_step, '.nc'

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

    class(lorenz96_type), intent(inout) :: this
    integer, intent(in) :: read_step ! Read in data for this time step

    integer :: ierr  ! return value of function

    ! local variables
    integer      :: size
    real(r8kind) :: forcing
    real(r8kind) :: delta_t
    real(r8kind) :: t
    integer      :: step
    character(len=128) :: filename

    ! assume normal termination
    ierr = 0

    ! Calculate name of file based on time step requested
    write(filename,'(A,I0.7,A)') 'lorenz96out_', read_step, '.nc'

    ! Read the model header
    ierr = netcdf_read_model_header(read_step, size, forcing, delta_t, t, step)

    ! Validate the input
    if (forcing /= this%forcing) then
      write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
      write(*,'(A,F7.3,A,F7.3)') '       Input file forcing =',forcing,', expecting ',this%forcing
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
    ierr = netcdf_read_model_data(read_step, this%size, this%location, this%state)

    netcdf_read_model_state = ierr

  end function netcdf_read_model_state


  !------------------------------------------------------------------
  ! ascii_read_model_header
  !------------------------------------------------------------------
  integer function ascii_read_model_header(read_step, size, forcing, delta_t, t, step)

    integer, intent(in)       :: read_step ! Read in data for this time step
    integer, intent(out)      :: size
    real(r8kind), intent(out) :: forcing
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
    write(filename, '(A,I0.7,A)') 'lorenz96out_', read_step, '.csv'

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
        case('model_forcing')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) forcing
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
          ! Ignore gloval settings we don't need
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
    integer :: position
    integer :: ignore
    integer :: i

    ! assume normal termination
    ierr = 0

    ! Construct name of input file
    write(filename, '(A,I0.7,A)') 'lorenz96out_', read_step, '.csv'

    ! Open the output csv file
    open(newunit=fileunit, file=trim(filename), form='formatted', status='old')

    ! Read global attributes
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

    class(lorenz96_type), intent(inout) :: this
    integer, intent(in) :: read_step ! Read in data for this time step

    integer :: ierr                  ! return value of function

    character(len=128) :: filename   ! name of output file
    integer :: fileunit
    character(len=80) :: line
    character(len=16) :: linefmt
    character(len=64) :: attr_name
    integer :: position
    integer :: ignore
    integer :: i
    integer      :: size
    real(r8kind) :: forcing
    real(r8kind) :: delta_t

    ! assume normal termination
    ierr = 0

    ! Construct name of input file
    write(filename, '(A,I0.7,A)') 'lorenz96out_', read_step, '.csv'

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
        case('model_forcing')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) forcing
          if (forcing /= this%forcing) then
            write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
            write(*,'(A,F7.3,A,F7.3)') '       Input file forcing =',forcing,', expecting ',this%forcing
            stop
          end if
        case('model_delta_t')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) delta_t
          if (delta_t /= this%delta_t) then
            write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
            write(*,'(A,F7.3,A,F7.3)') '       Input file delta_t =',delta_t,', expecting ',this%delta_t
            stop
          end if
        case('model_t')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',F)'
          read(line, linefmt) this%t
        case('model_step')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',I)'
          read(line, linefmt) this%step
        case('StateDim')
          write(linefmt, '(A2,I0,A3)') '(T', position + 1, ',I)'
          read(line, linefmt) size
          if (size /= this%size) then
            write(*,'(A,A)') 'ERROR: Incompatible input file: ', filename
            write(*,'(A,I,A,I)') '       Input file size =',size,', expecting ',this%size
            stop
          end if
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

    ! Read the coordinate, location, and state fields
    do i=1, this%size
      read(fileunit, *) ignore, this%location(i), this%state(i)
    end do

    ! Close the file
    close(fileunit)

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


end module lorenz96
