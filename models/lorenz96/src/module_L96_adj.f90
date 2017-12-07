module L96_ADJ

  use kind,           only : r8kind
  use Abstract_Model, only : abstract_model_type
  use L96_Config,     only : l96_config_type
  use L96_TL,         only : l96_tl_type

  implicit none

  private

  public :: l96_adj_type

  type, extends(abstract_model_type) :: l96_adj_type
    private
    type(l96_config_type)     :: config
    integer                   :: step
    real(r8kind)              :: clock
    real(r8kind), allocatable :: state(:)
    real(r8kind), allocatable :: location(:)
    real(r8kind), allocatable :: trajectory(:)
  contains
    final              :: destructor_adj
    procedure, private :: comp_dt
    procedure, private :: comp_dt_b
    procedure          :: adv_nsteps
    procedure          :: get_config
    procedure          :: get_location
    procedure          :: get_state
    procedure          :: get_trajectory
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

    integer :: j

    constructor_adj%config = config

    ! Initialize model step
    if (present(step)) then
      constructor_adj%step = step
    else
      constructor_adj%step = 0
    end if

    ! Initialize model clock
    constructor_adj%clock = constructor_adj%step * config%get_time_step()

    ! Localize the model domain
    allocate(constructor_adj%location(config%get_nx()))
    do j = 1, config%get_nx()
      constructor_adj%location(j) = (j - 1.0_r8kind) / config%get_nx()
    end do

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
  ! comp_dt
  !
  ! Private function to compute the time tendency for a l96_model_type
  ! object given a state, state.
  !------------------------------------------------------------------
  pure function comp_dt(this, state)

    class(l96_adj_type), intent(in)    :: this
    real(r8kind),          intent(in)    :: state(:)
    real(r8kind), dimension(size(state)) :: comp_dt

    integer :: nelements
    integer :: j

    ! Get the size of the array
    nelements = size(state)

    ! compute comp_dt(1)
    comp_dt(1) = (state(2) - state(nelements - 1)) * state(nelements) - state(1) + this%config%get_forcing()

    ! compute comp_dt(2)
    comp_dt(2) = (state(3) - state(nelements)) * state(1) - state(2) + this%config%get_forcing()

    ! compute comp_dt(3) thru comp_dt(size -1)
    do j = 3, nelements - 1
       comp_dt(j) = (state(j + 1) - state(j - 2)) * state(j - 1) - state(j) + this%config%get_forcing()
    end do

    ! compute comp_dt(size)
    comp_dt(nelements) = (state(1) - state(nelements - 2)) * state(nelements - 1) - state(nelements) + this%config%get_forcing()

  end function comp_dt


  !  Differentiation of comp_dt in reverse (adjoint) mode:
  !   gradient     of useful results: comp_dt state
  !   with respect to varying inputs: state
  !------------------------------------------------------------------
  ! comp_dt_b
  !
  ! Private function to compute the time tendency for a l96_model_type
  ! object given a state, state.
  !------------------------------------------------------------------
  subroutine comp_dt_b(this, state, stateb, comp_dtb)

    class(l96_adj_type), intent(inout)   :: this
    real(r8kind),        intent(   in)   :: state(:)
    real(r8kind),        intent(inout)   :: stateb(:)
    real(r8kind), dimension(size(state)) :: comp_dt
    real(r8kind), dimension(size(state)) :: comp_dtb

    integer      :: nelements
    integer      :: j
    real(r8kind) :: tempb
    real(r8kind) :: tempb0

    ! Get the size of the array
    nelements = size(state)

    tempb0 = state(nelements - 1) * comp_dtb(nelements)
    stateb(1) = stateb(1) + tempb0
    stateb(nelements - 2) = stateb(nelements - 2) - tempb0
    stateb(nelements - 1) = stateb(nelements - 1) + (state(1) - state(nelements - 2)) * comp_dtb(nelements)
    stateb(nelements) = stateb(nelements) - comp_dtb(nelements)
    comp_dtb(nelements) = 0.0_8
    do j = nelements - 1, 3, -1
      tempb = state(j - 1) * comp_dtb(j)
      stateb(j + 1) = stateb(j + 1) + tempb
      stateb(j - 2) = stateb(j - 2) - tempb
      stateb(j - 1) = stateb(j - 1) + (state(j + 1) - state(j - 2)) * comp_dtb(j)
      stateb(j) = stateb(j) - comp_dtb(j)
      comp_dtb(j) = 0.0_8
    end do
    stateb(3) = stateb(3) + state(1) * comp_dtb(2)
    stateb(nelements) = stateb(nelements) - state(1) * comp_dtb(2)
    stateb(1) = stateb(1) + (state(3) - state(nelements)) * comp_dtb(2)
    stateb(2) = stateb(2) - comp_dtb(2)
    comp_dtb(2) = 0.0_8
    stateb(2) = stateb(2) + state(nelements) * comp_dtb(1)
    stateb(nelements - 1) = stateb(nelements - 1) - state(nelements) * comp_dtb(1)
    stateb(nelements) = stateb(nelements) + (state(2) - state(nelements - 1)) * comp_dtb(1)
    stateb(1) = stateb(1) - comp_dtb(1)

  end subroutine comp_dt_b


  !------------------------------------------------------------------
  !        Generated by TAPENADE     (INRIA, Ecuador team)
  !  Tapenade 3.11 (r6148) - 16 Aug 2016 14:18
  !
  !  Differentiation of adv_nsteps in reverse (adjoint) mode:
  !   gradient     of useful results: state
  !   with respect to varying inputs: state
  !   RW status of diff variables: state:in-out
  !------------------------------------------------------------------
  ! adv_nsteps
  !
  ! Does n time step advances for lorenz 96 model
  ! using four-step rk time step
  !------------------------------------------------------------------
  subroutine adv_nsteps(this, nsteps)

    class(l96_adj_type), intent(inout) :: this
    integer            , intent(   in) :: nsteps

    real(r8kind), dimension(this%config%get_nx()) :: x1, x2, x3, x4, dx, inter
    real(r8kind), dimension(this%config%get_nx()) :: x1b, x2b, x3b, x4b, interb
    real(r8kind), dimension(this%config%get_nx()) :: result1
    real(r8kind), dimension(this%config%get_nx()) :: result1b
    real(r8kind), dimension(this%config%get_nx()) :: save_trajectory

    type(l96_tl_type) :: l96_tl

    integer :: step

    real(r8kind), dimension(this%config%get_nx(),this%config%get_nx()) :: mprime
    real(r8kind) :: time_step

    time_step = this%config%get_time_step()

    do step = 1, nsteps

      ! Initialize a L96 TL model with the current state and trajectory
      l96_tl = l96_tl_type(this%get_config(), state=this%state, trajectory=this%trajectory)

      ! Advance the L96 TL model one step
      call l96_tl%adv_nsteps(1)

      ! Save the current trajectory
      save_trajectory = this%trajectory

      ! Set the adjoint trajectory to the initial perturbation
      this%trajectory = l96_tl%get_trajectory() - this%trajectory

      ! Compute the first intermediate step
      result1 = this%comp_dt(this%state)
      x1 = time_step * result1
      inter = this%state + x1 / 2.0_r8kind

      ! Compute the second intermediate step
      result1 = this%comp_dt(inter)
      x2 = time_step * result1
      inter = this%state + x2 / 2.0_r8kind

      ! Compute the third intermediate step
      result1 = this%comp_dt(inter)
      x3 = time_step * result1
      inter = this%state + x3

      ! Update trajectory
      x1b = 0.0_8
      x2b = 0.0_8
      x3b = 0.0_8
      x4b = 0.0_8
      x1b = this%trajectory / 6.0_r8kind
      x2b = this%trajectory / 3.0_r8kind
      x4b = this%trajectory / 6.0_r8kind
      result1b = 0.0_8
      result1b = time_step * x4b
      interb = 0.0_8
      call this%comp_dt_b(inter, interb, result1b)

      x3b = interb + this%trajectory / 3.0_r8kind
      this%trajectory = this%trajectory + interb
      result1b = 0.0_8
      result1b = time_step * x3b
      inter = this%state + x2 / 2.0_r8kind
      interb = 0.0_8
      call this%comp_dt_b(inter, interb, result1b)

      this%trajectory = this%trajectory + interb
      x2b = x2b + interb / 2.0_r8kind
      result1b = 0.0_8
      result1b = time_step * x2b
      inter = this%state + x1 / 2.0_r8kind
      interb = 0.0_8
      call this%comp_dt_b(inter, interb, result1b)

      this%trajectory = this%trajectory + interb
      x1b = x1b + interb / 2.0_r8kind
      result1b = 0.0_8
      result1b = time_step * x1b
      call this%comp_dt_b(this%state, this%trajectory, result1b)

      ! Update the trajectory
      this%trajectory = save_trajectory - this%trajectory

      ! Update state with original trajectory
      this%state = this%state - this%trajectory

! Lidia's MATRIX formulation with compensation for an incorrect input trajectory
!44    FORMAT (A8,6F10.6)
!      call buildMprime(this%state, mprime)
!      x2 = -this%config%get_time_step() * matmul(mprime,this%trajectory)
!      x1 = -this%config%get_time_step() * matmul(transpose(mprime), x2)
!      this%trajectory = this%trajectory + x2 + x1
!      this%state = this%state + this%trajectory

! MATRIX formulation that assumes a correct input trajectory is available for the adjoint
! The above MATRIX forumulation compensates for an input trajectory valid at time t+1
! of the correct one.  This formulation assumes the correct trajectory is available
!44    FORMAT (A8,6F10.6)
!      call buildMprime(this%state, mprime)
!      this%state = this%state + this%trajectory
!      this%trajectory = this%trajectory + this%config%get_time_step() * matmul(transpose(mprime), this%trajectory) 

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
  pure function get_location(this)

    class(l96_adj_type), intent(in) :: this
    real(r8kind), dimension(this%config%get_nx()) :: get_location

    get_location = this%location

  end function get_location


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
  ! get_trajectory
  !
  ! Return the model trajectory vector
  !------------------------------------------------------------------
  pure function get_trajectory(this)

    class(l96_adj_type), intent(in) :: this
    real(r8kind), dimension(this%config%get_nx()) :: get_trajectory

    get_trajectory = this%trajectory

  end function get_trajectory


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
