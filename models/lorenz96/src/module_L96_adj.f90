module L96_ADJ

  use kind,       only : r8kind
  use model,      only : model_type
  use L96_Config, only : l96_config_type

  implicit none

  private

  public :: l96_adj_type

  type, extends(model_type) :: l96_adj_type
    private
    type(l96_config_type)     :: config
    integer                   :: step
    real(r8kind)              :: clock
    real(r8kind), allocatable :: state(:)
    real(r8kind), allocatable :: trajectory(:)
  contains
    final              :: destructor_adj
    procedure, private :: comp_dt
    procedure          :: adv_nsteps
    procedure          :: get_config
!    procedure          :: get_location
    procedure          :: get_state
    procedure          :: get_step
    procedure          :: get_clock
    procedure          :: print
!    procedure          :: interpolate
  end type l96_adj_type

  interface l96_adj_type
    procedure constructor_adj
  end interface


contains


  !------------------------------------------------------------------
  ! constructor_adj
  !
  ! Returns an initialized l96_adj_type object
  !------------------------------------------------------------------
  type(l96_adj_type) function constructor_adj(config, state, trajectory, step)

    class(l96_config_type), intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:)
    real(r8kind), optional, intent(in) :: trajectory(:)
    integer,      optional, intent(in) :: step

    constructor_adj%config = config

    ! Initialize model step
    if (present(step)) then
      constructor_adj%step = step
    else
      constructor_adj%step = 0
    end if

    ! Initialize model clock
    constructor_adj%clock = constructor_adj%step * config%get_time_step()

    ! Initialize model state
    allocate(constructor_adj%state(config%get_nx()))
    if (present(state)) then
      constructor_adj%state(:) = state(:)
    else
      constructor_adj%state(:) = config%get_forcing()
      constructor_adj%state(1) = 1.001_r8kind * config%get_forcing()
    end if

    ! Initialize model trajectory
    allocate(constructor_adj%trajectory(config%get_nx()))
    if (present(trajectory)) then
      constructor_adj%trajectory(:) = trajectory(:)
    else
      constructor_adj%trajectory(:) = constructor_adj%state(:)
    end if

  end function constructor_adj


  !------------------------------------------------------------------
  ! destructor_adj
  !
  ! Deallocates pointers used by a l96_adj_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_adj(this)

    type(l96_adj_type), intent(inout) :: this

    ! No pointers in l96_adj_type object so we do nothing

  end subroutine destructor_adj


  !------------------------------------------------------------------
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.11 (r6148) - 16 Aug 2016 14:18
  !
  !  Differentiation of comp_dt in reverse (adjoint) mode:
  !   gradient     of useful results: dt x
  !   with respect to varying inputs: dt x
  !------------------------------------------------------------------
  ! comp_dt_b
  !
  ! Private routine to compute the time tendency of a lorenz96_type object
  ! given a state, x, and return it in dt.
  !------------------------------------------------------------------
  subroutine comp_dt(this, state, stateb, dt, dtb)

    class(l96_adj_type), intent(   in) :: this
    real(r8kind),        intent(   in) :: state(:)
    real(r8kind),        intent(inout) :: stateb(size(state))
    real(r8kind),        intent(   in) :: dt(:)
    real(r8kind),        intent(inout) :: dtb(size(stateb))

    integer      :: j
    integer      :: nelements
    real(r8kind) :: tempb

    ! Get the size of the array
    nelements = size(state)

    stateb(1) = stateb(1) + state(nelements - 1) * dtb(nelements)
    stateb(nelements - 2) = stateb(nelements - 2) - state(nelements - 1) * dtb(nelements)
    stateb(nelements - 1) = stateb(nelements - 1) + (state(1) - state(nelements - 2)) * dtb(nelements)
    stateb(nelements) = stateb(nelements) - dtb(nelements)
    dtb(nelements) = 0.0

    do j = nelements - 1, 3, -1
      tempb = state(j-1) * dtb(j)
      stateb(j + 1) = stateb(j + 1) + tempb
      stateb(j - 2) = stateb(j - 2) - tempb
      stateb(j - 1) = stateb(j - 1) + (state(j + 1) - state(j - 2)) * dtb(j)
      stateb(j) = stateb(j) - dtb(j)
      dtb(j) = 0.0
    end do

    stateb(3) = stateb(3) + state(1) * dtb(2)
    stateb(nelements) = stateb(nelements) - state(1) * dtb(2)
    stateb(1) = stateb(1) + (state(3) - state(nelements)) * dtb(2)
    stateb(2) = stateb(2) - dtb(2)
    dtb(2) = 0.0

    stateb(2) = stateb(2) + state(nelements) * dtb(1)
    stateb(nelements - 1) = stateb(nelements - 1) - state(nelements) * dtb(1)
    stateb(nelements) = stateb(nelements) + (state(2) - state(nelements - 1)) * dtb(1)
    stateb(1) = stateb(1) - dtb(1)
    dtb(1) = 0.0

  end subroutine comp_dt


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
  subroutine adv_nsteps(this, nsteps)

    class(l96_adj_type), intent(inout) :: this
    integer, intent(in) :: nsteps

    real(r8kind), dimension(this%config%get_nx()) :: x1, x2, x3, x4, dx, inter
    real(r8kind), dimension(this%config%get_nx()) :: x1b, x2b, x3b, x4b, dxb, interb

    integer :: step
    real(r8kind), dimension(this%config%get_nx(),this%config%get_nx()) :: mprime

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
      x2 = -this%config%get_time_step() * matmul(mprime,this%trajectory)
      x1 = -this%config%get_time_step() * matmul(transpose(mprime), x2)
      this%trajectory = this%trajectory + x2 + x1
      this%state = this%state + this%trajectory

      ! Increment time step
      this%clock = this%clock - this%config%get_time_step()
      this%step = this%step - 1

    end do

  end subroutine adv_nsteps


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
  ! get_config
  !
  ! Return the model configuration
  !------------------------------------------------------------------  
  pure type(l96_config_type) function get_config(this)

    class(l96_adj_type), intent(in) :: this

    get_config = this%config

  end function get_config


  !------------------------------------------------------------------  
  ! get_location
  !
  ! Return the model location vector
  !------------------------------------------------------------------  
!  pure function get_location(this)
!
!    class(l96_adj_type), intent(in) :: this
!    real(r8kind), dimension(this%config%get_nx()) :: get_location
!
!    get_location = this%location    
!
!  end function get_location


  !------------------------------------------------------------------  
  ! get_state
  !
  ! Return the model state vector
  !------------------------------------------------------------------  
  pure function get_state(this)

    class(l96_adj_type), intent(in) :: this
    real(r8kind), dimension(this%config%get_nx()) :: get_state

    get_state = this%state    

  end function get_state


  !------------------------------------------------------------------  
  ! get_step
  !
  ! Return the model step
  !------------------------------------------------------------------  
  pure integer function get_step(this)

    class(l96_adj_type), intent(in) :: this

    get_step = this%step

  end function get_step



  !------------------------------------------------------------------  
  ! get_clock
  !
  ! Return the model clock
  !------------------------------------------------------------------  
  pure real(r8kind) function get_clock(this)

    class(l96_adj_type), intent(in) :: this

    get_clock = this%clock

  end function get_clock


  !------------------------------------------------------------------  
  ! print
  !
  ! Print the contents of this model object
  !------------------------------------------------------------------  
  subroutine print(this)

    class(l96_adj_type), intent(in) :: this

    character(len=32) :: numstr, varstr, indexstr
    integer :: i

    call this%config%print()

    write(numstr,'(I)') this%step
    write(*,'(A20,A)') 'step = ', adjustl(numstr)
    write(numstr,'(F16.13)') this%clock
    write(*,'(A20,A)') 'clock = ', adjustl(numstr)
    do i=1, this%config%get_nx()
      write(numstr,'(F16.13)') this%state(i)
      write(indexstr,'(I)') i
      write(*,'(A20,A)') 'state(' // trim(adjustl(indexstr)) // ') = ', numstr
    end do

  end subroutine print


  !------------------------------------------------------------------  
  ! Interpolates from state vector state to the location. 
  !------------------------------------------------------------------  
!  subroutine interpolate(this, location, state_val)
!
!    class(l96_model_type), intent(in) :: this
!    real(r8kind), intent(in)    :: location
!    real(r8kind), intent(out)   :: state_val
!
!    integer :: lower_index, upper_index, i
!    real(r8kind) :: lctn, lctnfrac
!
!    ! Scale the location to the size of the domain
!    lctn = this%size * location
!
!    ! Compute grid indices bounding the location
!    lower_index = int(lctn) + 1
!    upper_index = lower_index + 1
!    if(lower_index > this%size) lower_index = lower_index - this%size
!    if(upper_index > this%size) upper_index = upper_index - this%size
!
!    ! Interpolate model value at the location
!    lctnfrac = lctn - int(lctn)
!    state_val = (1.0_r8kind - lctnfrac) * this%get_state()(lower_index) + lctnfrac * this%get_state()(upper_index)
! 
!  end subroutine interpolate


end module L96_ADJ
