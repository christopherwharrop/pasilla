program QG

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use kind
  use QG_Config
  use QG_Model
  use QG_Reader
  use QG_Writer
  use QGTL_Writer
  use QGADJ_Writer
!  use gptl

  implicit none

  integer, parameter :: increment=1  ! Controls the number of steps to move the models each iteration
  integer :: start_step = 0
  integer :: spinup_steps = 720
  integer :: run_steps = 1440
  integer :: output_interval_steps = 3
  logical :: readstart = .true.

  namelist /runtime/ start_step, spinup_steps, run_steps, output_interval_steps, readstart

  type(qg_config_type) :: config
  type(qg_model_type)  :: model
  type(qg_tl_type)     :: model_tl
  type(qg_adj_type)    :: model_adj
  type(qg_reader_type) :: reader
  type(qg_writer_type) :: writer
  type(qgtl_writer_type) :: writer_tl
  type(qgadj_writer_type) :: writer_adj
  character(len=64)    :: filename  ! output/input filename

  real(r8kind), allocatable :: old_state(:,:)
  real(r8kind), allocatable :: new_state(:,:)
  real(r8kind), allocatable :: tplusone_state(:,:)
  real(r8kind), allocatable :: input_trajectory(:,:)
  real(r8kind), allocatable :: output_trajectory(:,:)
  integer                   :: nsh2, nvl
  integer                   :: step, ret

!  ret = gptlsetutr (gptlnanotime)
!  ret = gptlinitialize ()
!  ret = gptlstart ('QG')

  ! Read namelist from stdin
  read(stdin,nml=runtime)
  rewind(stdin)

  ! If this is a restart, initialize model from a restart file
  if (start_step /= 0) then

    ! Initialize a reader to read restart file
    reader = qg_reader_type("NETCDF")

    ! Construct name of restart file
    write(filename,'(A,I0.7)') 'qgout_', start_step

    ! Load restart file into new model instance
    call reader%read(model, filename)

  ! Otherwise create a new model from the namelist
  else
    ! Create a configuration as specified by the namelist
    config = qg_config_type(stdin)

    ! Create a model with the namelist configuration
    model = qg_model_type(config)

  end if

  ! Get the model's spectral dimensions
  nsh2 = model%get_nsh2()
  nvl = model%get_nvl()

  ! Allocate arrays for tracking state and computing trajectories
  allocate(old_state(nsh2, nvl))
  allocate(new_state(nsh2, nvl))
  allocate(tplusone_state(nsh2,nvl))
  allocate(input_trajectory(nsh2, nvl))
  allocate(output_trajectory(nsh2, nvl))

  ! Create a writer to write model state
  writer = qg_writer_type('NETCDF')
  writer_tl = qgtl_writer_type('NETCDF')
  writer_adj = qgadj_writer_type('NETCDF')

  ! Write the initial state
  write(filename,'(A,I0.7)') 'qg_model_out_', model%get_step()
  call writer%write(model, filename)

  ! Spinup the forward model (if required)
  if (spinup_steps > 0) then

    write(*,*) 'Integrating forward model spinup steps: ', spinup_steps
    call model%adv_nsteps(spinup_steps)

    ! Write the post-spinup state
    write(filename,'(A,I0.7)') 'qg_model_out_', model%get_step()
    call writer%write(model, filename)
  end if

  ! Advance the forward model 1 step to t=t+1 so we can compute tl trajectory
  old_state = model%get_psi()
  call model%adv_nsteps(1)
  new_state = model%get_psi()

  ! Run the model
  write(*,*) 'Integrating forward model trajectory steps: ', run_steps

  do step = 1, run_steps * increment, increment

    ! Instantiate a tangent linear model with the forward model's state at t=t
    model_tl = qg_tl_type(config, state = old_state, trajectory = new_state - old_state, step = step - 1)

    ! Advance the tangent linear model trajectory forward three steps to t=t+increment
    call model_tl%adv_nsteps(increment)

    ! Advance the forward model increment-1 steps to t=t+increment so we can output in increments of increment steps
    call model%adv_nsteps(increment - 1)
    old_state = model%get_psi()

    ! Output forward model state at t=t+increment
    write(filename,'(A,I0.7)') 'qg_model_out_', model%get_step()
    call writer%write(model, filename)

    ! Advance the forward model 1 step to t=t+increment+1 so we can compute trajectory for adj
    call model%adv_nsteps(1)
    new_state = model%get_psi()

    ! Instantiate an adjoint with the new state at t=t+increment
!    model_adj = qg_adj_type(config, state = old_state, trajectory = -(new_state - old_state), step = step - 1 + increment)
    model_adj = qg_adj_type(config, state = new_state, trajectory = old_state - new_state, step = step + increment)

    ! Advance adjoint trajectory (backward) increment steps to t=t
!    call model_adj%adv_nsteps(increment)
    call model_adj%adv_nsteps(increment+1)

    ! Output tangent linear state at t=t+increment
    write(filename,'(A,I0.7)') 'qg_tl_out_', model_tl%get_step()
    call writer_tl%write(model_tl, filename)

    ! Output adjoint state at t=t
    write(filename,'(A,I0.7)') 'qg_adj_out_', model_adj%get_step()
    call writer_adj%write(model_adj, filename)

  enddo

!  ret = gptlstop ('QG')
!  ret = gptlpr (0)
!  ret = gptlfinalize ()

 
end program QG
