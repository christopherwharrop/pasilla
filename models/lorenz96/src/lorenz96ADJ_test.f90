program lorenz96ADJ_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use L96_Model, only : l96_model_type
  use L96_Config, only : l96_config_type
  use L96_TL, only    : l96_tl_type
  use L96_ADJ, only    : l96_adj_type
  use L96_Reader
  use kind, only     : r8kind, i8kind

  implicit none


  type(l96_model_type)  :: l96
  type(l96_tl_type)     :: l96tl
  type(l96_adj_type)    :: l96adj
  type(l96_reader_type) :: reader
  type(l96_config_type) :: config

  real(r8kind), allocatable :: state(:), delta(:), x(:), y(:)

  ! Define namelists and default values
  integer          :: nx    = 40
  real(r8kind)     :: forcing = 8.00_r8kind
  real(r8kind)     :: time_step = 0.05_r8kind
  namelist /params/  nx, forcing, time_step
  integer          :: start_step = 0
  integer          :: spinup_steps = 1000000
  integer          :: run_steps   = 1000
  integer          :: output_interval_steps = 100
  character(len=6) :: io_format = 'NETCDF'
  namelist /runtime/ start_step, run_steps, output_interval_steps, io_format

  integer :: t            ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: filename

  ! Read namelists from stdin
!  read(stdin,nml=params)
!  read(stdin,nml=runtime)

  ! If this is a restart, initialize model from a restart file
  if (start_step /= 0) then

    ! Initialize a reader to read restart file
    reader = l96_reader_type("NETCDF")

    ! Construct name of restart file
    write(filename,'(A,I0.7)') 'l96out_', start_step

    ! Load restart file into new model instance
    call reader%read(l96, filename)

  ! Otherwise create a new model from the namelist
  else
    ! Create a configuration as specified by the namelist
    config = l96_config_type(stdin)

    ! Create a model with the namelist configuration
    l96 = l96_model_type(config)

  end if

  nx = config%get_nx()
  allocate(state(nx))
  allocate(x(nx))
  allocate(y(nx))
  allocate(delta(nx))

  ! Spinup the forward model (if required)
  if (spinup_steps > 0) then

    write(*,*) 'Integrating forward model spinup steps: ', spinup_steps
    call l96%adv_nsteps(spinup_steps)

  end if

  ! Advance the forward model 1 step to t=t+1 so we can compute tl trajectory
  state = l96%get_state()
  call l96%adv_nsteps(10000)
  x = (l96%get_state() - state)
  x = x * 0.001
  call l96%adv_nsteps(10000)
  y = (l96%get_state() - state)
  y = y * 0.001

  l96tl = l96_tl_type(config, state=state, trajectory=x)
  call l96tl%adv_nsteps(1)

  l96adj = l96_adj_type(config, state=state, trajectory=y)
  call l96adj%adv_nsteps(1)

  write(*,'(F20.15)') dot_product(l96tl%get_trajectory(), y) / dot_product(x, l96adj%get_trajectory())

end program lorenz96ADJ_Test
