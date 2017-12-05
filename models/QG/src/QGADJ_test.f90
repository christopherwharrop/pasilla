program QGADJ_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use QG_Model, only : qg_model_type, qg_tl_type, qg_adj_type
  use QG_Config, only : qg_config_type
  use QG_Reader
  use kind, only     : r8kind, i8kind

  implicit none


  integer              :: nsh, nsh2, nvl
  type(qg_model_type)  :: qg
  type(qg_tl_type)     :: qgtl
  type(qg_adj_type)    :: qgadj
  type(qg_reader_type) :: reader
  type(qg_config_type) :: config

  real(r8kind), allocatable :: state(:,:), delta(:,:), x(:,:), y(:,:)

  ! Define namelists and default values
  integer :: start_step = 0
  integer :: spinup_steps = 720
  integer :: run_steps = 1440
  integer :: output_interval_steps = 3
  logical :: readstart = .true.
  namelist /runtime/ start_step, spinup_steps, run_steps, output_interval_steps, readstart

  integer :: t            ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: filename

  ! Read namelists from stdin
!  read(stdin,nml=params)
!  read(stdin,nml=runtime)

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
  nvl = qg%get_nvl()
  allocate(state(nsh2,nvl))
  allocate(x(nsh2,nvl))
  allocate(y(nsh2,nvl))
  allocate(delta(nsh2,nvl))

  ! Spinup the forward model (if required)
  if (spinup_steps > 0) then

    write(*,*) 'Integrating forward model spinup steps: ', spinup_steps
    call qg%adv_nsteps(spinup_steps)

  end if

  ! Advance the forward model 1 step to t=t+1 so we can compute tl trajectory
  state = qg%get_psi()
  call qg%adv_nsteps(1000)
  x = (qg%get_psi() - state)
  x = x * 0.001
  call qg%adv_nsteps(1000)
  y = (qg%get_psi() - state)
  y = y * 0.001

  qgtl = qg_tl_type(config, state=state, trajectory=x)
  call qgtl%adv_nsteps(1)

  qgadj = qg_adj_type(config, state=state, trajectory=y)
  call qgadj%adv_nsteps(1)

!  write(*,'(F20.15)') dot_product(qgtl%get_trajectory(), y) / dot_product(x, qgadj%get_trajectory())
  write(*,'(F20.15)') dot_product(reshape(qgtl%get_trajectory(),(/nsh2*nvl/)), reshape(y,(/nsh2*nvl/))) / dot_product(reshape(x,(/nsh2*nvl/)), reshape(qgadj%get_trajectory(),(/nsh2*nvl/)))

end program QGADJ_Test
