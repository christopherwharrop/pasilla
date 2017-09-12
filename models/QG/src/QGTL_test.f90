program QGTL_test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use QG_Model
  use QG_Config
  use QG_Model_TL
  use QG_Reader
  use kind

  implicit none

  integer              :: nsh, nsh2, nvl
  type(qg_model_type)  :: qg, qg_delta
  type(qg_tl_type)     :: qgtl
  type(qg_reader_type) :: reader
  type(qg_config_type) :: config

  real(r8kind), allocatable :: state(:,:), delta(:,:), m(:,:), m_delta(:,:), mprime_delta(:,:)

  ! Define namelists and default values
  integer :: start_step = 0
  integer :: spinup_steps = 720
  integer :: run_steps = 1440
  integer :: output_interval_steps = 3
  logical :: readstart = .true.
  namelist /runtime/ start_step, spinup_steps, run_steps, output_interval_steps, readstart

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
    reader = qg_reader_type("NETCDF")

    ! Construct name of restart file
    write(filename,'(A,I0.7)') 'qgout_', start_step

    ! Load restart file into new model instance
    call reader%read(qg, filename)

  ! Otherwise create a new model from the namelist
  else
    ! Create a configuration as specified by the namelist
    config = qg_config_type(stdin)

    ! Create a model with the namelist configuration
    qg = qg_model_type(config)

  end if

  nsh2 = qg%get_nsh2()
  nsh = nsh2 / 1
  nvl = qg%get_nvl()

  allocate(state(nsh2, nvl))
  allocate(delta(nsh2, nvl))
  allocate(m(nsh2, nvl))
  allocate(m_delta(nsh2, nvl))
  allocate(mprime_delta(nsh2, nvl))

  ! Spinup the forward model (if required)
  if (spinup_steps > 0) then

    write(*,*) 'Integrating forward model spinup steps: ', spinup_steps
    call qg%adv_nsteps(spinup_steps)

  end if

  ! Advance the forward model 1 step to t=t+1 so we can compute tl trajectory
  state = qg%get_psi()
  call qg%adv_nsteps(1000)
  delta = (qg%get_psi() - state)
!  call random_number(delta)
  delta = delta * 0.001
  write(*,'(2A15)') 'min dx', 'max dx'
  write(*,'(2F15.10)') minval(delta), maxval(delta)

  ! Write column headers
  write(*,'(A12,1x,A12,19x,A)') "Step", "Lambda", "( M(x + lambda * dx) - M(x) ) / M'(lambda * dx)"

  ! Loop over digits of precision
  do d = 1, digits

    ! Create a Lorenz96 model configured as specified in namelist for M(x)
    qg = qg_model_type(config, state=state)

    ! Create a Lorenz96 model configured as specified in namelist for M(x + lambda * dx)
    qg_delta = qg_model_type(config, state=state + lambda * delta)

    ! Create a Lorenz96TL model configured as specified in namelist to calculate dx
    qgtl = qg_tl_type(config, state=state, trajectory=lambda * delta)

    ! Advance QG, QGdelta, and TLdelta
    call qg%adv_nsteps(1)
    call qg_delta%adv_nsteps(1)
    call qgtl%adv_nsteps(1)

    ! Calculate and print test metric
    m = qg%get_psi()
    m_delta = qg_delta%get_psi()
    mprime_delta = qgtl%get_trajectory()
    write(*,'(I,4F)') qg%get_step(), lambda, (m_delta(nsh,1) - m(nsh,1)) / mprime_delta(nsh,1)

    ! Increase precision
    lambda = lambda / 10.0_r8kind

  end do
 
end program QGTL_test
