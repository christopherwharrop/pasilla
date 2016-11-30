program L96

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use kind,      only  : r8kind, i8kind
  use L96_Model, only  : l96_model_type
  use L96_Config, only : l96_config_type
  use L96_Writer, only : l96_writer_type

  implicit none

  integer          :: start_step = 0
  integer          :: run_steps   = 1000
  integer          :: output_interval_steps = 100
  character(len=6) :: io_format = 'NETCDF'
  namelist /runtime/ start_step, run_steps, output_interval_steps, io_format

  type(L96_config_type) :: config ! Lorenz96 configuration
  type(L96_model_type)  :: model  ! Lorenz96 model
  type(L96_writer_type) :: writer ! Lorenz96 writer

  integer :: t     ! Loop index
  integer :: ierr  ! Error code 
  character(len=21) :: restart_file

  ! Create a Lorenz96 model configuration as specified in namelist
!  config = l96_config_type(size, delta_t, forcing)
!  config = l96_config_type('lorenz96.namelist')
  config = l96_config_type(stdin)

  ! Read namelist from stdin
  read(stdin,nml=runtime)

  ! Create and configure a Lorenz96 model 
  model = l96_model_type(config)

call model%print()

  ! Create a writer to write model state
  writer = l96_writer_type(io_format)

  ! Write out model state if needed
  if (output_interval_steps <= run_steps) ierr = writer%write(model)

  ! Run the model
  do t = 1, run_steps, output_interval_steps

    ! Advance the model to next output interval
    call model%adv_nsteps(min(output_interval_steps, run_steps))

    ! Write out model state if needed
    if (output_interval_steps <= run_steps) ierr = writer%write(model)

call model%print()

  end do

end program
