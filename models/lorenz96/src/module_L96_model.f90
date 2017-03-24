module L96_Model

  use Kind,           only : r8kind
  use Abstract_Model, only : abstract_model_type
  use L96_Config,     only : l96_config_type

  implicit none

  private

  public :: l96_model_type

  type, extends(abstract_model_type) :: l96_model_type
    private
    type(l96_config_type)     :: config
    real(r8kind), allocatable :: state(:)
    real(r8kind), allocatable :: location(:)
    integer                   :: step
    real(r8kind)              :: clock
  contains
    final              :: destructor_model
    procedure, private :: comp_dt
    procedure          :: adv_nsteps
    procedure          :: get_config
    procedure          :: get_location
    procedure          :: get_state
    procedure          :: get_step
    procedure          :: get_clock
    procedure          :: print
    procedure          :: interpolate
  end type l96_model_type

  interface l96_model_type
    procedure constructor_model
  end interface

contains

  !------------------------------------------------------------------
  ! constructor_model
  !
  ! Returns an initialized l96_model_type object
  !------------------------------------------------------------------
  type(l96_model_type) function constructor_model(config, state, step)

    class(l96_config_type), intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:)
    integer, optional, intent(in) :: step

    integer :: j

    ! Initialize model parameters    
    constructor_model%config = config

    ! Localize the model domain
    allocate(constructor_model%location(config%get_nx()))
    do j = 1, config%get_nx()
      constructor_model%location(j) = (j - 1.0_r8kind) / config%get_nx()
    end do

    ! Initialize model state
    allocate(constructor_model%state(config%get_nx()))
    if (present(state)) then
      constructor_model%state(:) = state(:)
    else
      constructor_model%state(:) = config%get_forcing()
      constructor_model%state(1) = 1.001_r8kind * config%get_forcing()
    end if

    ! Initialize model step
    if (present(step)) then
      constructor_model%step = step
    else
      constructor_model%step = 0
    end if

    ! Initialize model clock
    constructor_model%clock = constructor_model%step * config%get_time_step()

  end function constructor_model


  !------------------------------------------------------------------
  ! destructor_model
  !
  ! Deallocates pointers used by a l96_model_type object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor_model(this)

    type(l96_model_type), intent(inout) :: this

    ! No pointers in l96_model_type object so we do nothing

  end subroutine destructor_model


  !------------------------------------------------------------------
  ! comp_dt
  !
  ! Private function to compute the time tendency for a l96_model_type 
  ! object given a state, state.
  !------------------------------------------------------------------
  pure function comp_dt(this, state)

    class(l96_model_type), intent(in)    :: this
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


  !------------------------------------------------------------------
  ! adv_nsteps
  !
  ! Does n time step advances for lorenz 96 model
  ! using four-step rk time step
  !------------------------------------------------------------------
  subroutine adv_nsteps(this, nsteps)

    class(l96_model_type),                         intent(inout) :: this
    integer,                                       intent(   in) :: nsteps

    real(r8kind), dimension(this%config%get_nx()) :: x1, x2, x3, x4, dx, inter
    integer :: step
    
    do step = 1, nsteps

      !  Compute the first intermediate step
      x1    = this%config%get_time_step() * this%comp_dt(this%state)
      inter = this%state + x1 / 2.0_r8kind

      !  Compute the second intermediate step
      x2    = this%config%get_time_step() * this%comp_dt(inter)
      inter = this%state + x2 / 2.0_r8kind

      !  Compute the third intermediate step
      x3    = this%config%get_time_step() * this%comp_dt(inter)
      inter = this%state + x3

      !  Compute fourth intermediate step
      x4 = this%config%get_time_step() * this%comp_dt(inter)

      !  Compute new value for state
      this%state = this%state + x1 / 6.0_r8kind + x2 / 3.0_r8kind + x3 / 3.0_r8kind + x4 / 6.0_r8kind

      ! Increment time step
      this%clock = this%clock + this%config%get_time_step()
      this%step = this%step + 1

    end do

  end subroutine adv_nsteps


  !------------------------------------------------------------------
  ! get_config
  !
  ! Return the model configuration
  !------------------------------------------------------------------
  pure type(l96_config_type) function get_config(this)

    class(l96_model_type), intent(in) :: this

    get_config = this%config

  end function get_config


  !------------------------------------------------------------------
  ! get_location
  !
  ! Return the model location vector
  !------------------------------------------------------------------
  pure function get_location(this)

    class(l96_model_type), intent(in) :: this
    real(r8kind), dimension(this%config%get_nx()) :: get_location

    get_location = this%location

  end function get_location


  !------------------------------------------------------------------
  ! get_state
  !
  ! Return the model state vector
  !------------------------------------------------------------------
  pure function get_state(this)

    class(l96_model_type), intent(in) :: this
    real(r8kind), dimension(this%config%get_nx()) :: get_state

    get_state = this%state

  end function get_state


  !------------------------------------------------------------------  
  ! get_step
  !
  ! Return the model step
  !------------------------------------------------------------------  
  pure integer function get_step(this)

    class(l96_model_type), intent(in) :: this

    get_step = this%step

  end function get_step



  !------------------------------------------------------------------  
  ! get_clock
  !
  ! Return the model clock
  !------------------------------------------------------------------  
  pure real(r8kind) function get_clock(this)

    class(l96_model_type), intent(in) :: this

    get_clock = this%clock

  end function get_clock


  !------------------------------------------------------------------  
  ! print
  !
  ! Print the contents of this model object
  !------------------------------------------------------------------  
  subroutine print(this)

    class(l96_model_type), intent(in) :: this

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
  subroutine interpolate(this, location, state_val)

    class(l96_model_type), intent(in) :: this
    real(r8kind), intent(in)    :: location
    real(r8kind), intent(out)   :: state_val

    integer      :: lower_index, upper_index, i, size
    real(r8kind) :: lctn, lctnfrac

    ! Get the number of grid points
    size = this%config%get_nx()

    ! Scale the location to the size of the domain
    lctn = size * location

    ! Compute grid indices bounding the location
    lower_index = int(lctn) + 1
    upper_index = lower_index + 1
    if(lower_index > size) lower_index = lower_index - size
    if(upper_index > size) upper_index = upper_index - size

    ! Interpolate model value at the location
    lctnfrac = lctn - int(lctn)
    state_val = (1.0_r8kind - lctnfrac) * this%state(lower_index) + lctnfrac * this%state(upper_index)

  end subroutine interpolate


end module L96_Model
