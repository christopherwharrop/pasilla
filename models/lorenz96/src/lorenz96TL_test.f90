program lorenz96TL_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use L96_Model, only : l96_model_type
  use L96_Config, only : l96_config_type
  use L96_TL, only    : l96_tl_type
  use L96_Reader
  use kind, only     : r8kind, i8kind

  implicit none


  type(l96_model_type)  :: l96, l96_delta
  type(l96_tl_type)     :: l96tl
  type(l96_reader_type) :: reader
  type(l96_config_type) :: config

  real(r8kind), allocatable :: state(:), delta(:), m(:), m_delta(:), mprime_delta(:)

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

  real(r8kind) :: lambda  ! scaling factor
  integer :: d            ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: filename
  integer, parameter :: digits = 12  ! Number of significant digits to check

  ! Read namelists from stdin
!  read(stdin,nml=params)
!  read(stdin,nml=runtime)

  ! Initialize lambda
  lambda = 1.0_r8kind

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
  allocate(delta(nx))
  allocate(m(nx))
  allocate(m_delta(nx))
  allocate(mprime_delta(nx))

  ! Spinup the forward model (if required)
  if (spinup_steps > 0) then

    write(*,*) 'Integrating forward model spinup steps: ', spinup_steps
    call l96%adv_nsteps(spinup_steps)

  end if

  ! Advance the forward model 1 step to t=t+1 so we can compute tl trajectory
  state = l96%get_state()
  call l96%adv_nsteps(10000)
  delta = (l96%get_state() - state)
  delta = delta * 0.001
  write(*,'(2A15)') 'min dx', 'max dx'
  write(*,'(2F15.10)') delta(minloc(delta)), delta(maxloc(delta))

  ! Write column headers
  write(*,'(A12,1x,A12,19x,A)') "Step", "Lambda", "( M(x + lambda * dx) - M(x) ) / M'(lambda * dx)"

  ! Loop over digits of precision
  do d = 1, digits

    ! Create a Lorenz96 model configured as specified in namelist for M(x)
    l96 = l96_model_type(config, state=state)

    ! Create a Lorenz96 model configured as specified in namelist for M(x + lambda * dx)
    l96_delta = l96_model_type(config, state=state + lambda * delta)

    ! Create a Lorenz96TL model configured as specified in namelist to calculate dx
    l96tl = l96_tl_type(config, state=state, trajectory=lambda * delta)

    ! Advance L96, L96delta, and TLdelta
    call l96%adv_nsteps(1)
    call l96_delta%adv_nsteps(1)
    call l96tl%adv_nsteps(1)

    ! Calculate and print test metric
    m = l96%get_state()
    m_delta = l96_delta%get_state()
    mprime_delta = l96tl%get_trajectory()
    write(*,'(I,4F)') l96%get_step(), lambda, (m_delta(20) - m(20)) / mprime_delta(20)

    ! Increase precision
    lambda = lambda / 10.0_r8kind

  end do

end program
