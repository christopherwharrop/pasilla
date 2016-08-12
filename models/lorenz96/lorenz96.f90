program lorenz96Model

  ! Get unit numbers for stdin, stdout, stderr
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use module_lorenz96, only : lorenz96  
  use module_kind, only     : r8kind, i8kind

  implicit none

  ! Lorenz96 model object
  type(lorenz96) :: L96

  ! Define namelists and default values
  integer         :: size    = 40
  real(r8kind)    :: forcing = 8.00_r8kind
  real(r8kind)    :: delta_t = 0.05_r8kind
  integer         :: start_step = 1
  integer         :: run_steps   = 1000
  integer         :: output_interval_steps = 100
  namelist /params/  size, forcing, delta_t 
  namelist /runtime/ start_step, run_steps, output_interval_steps

  integer :: t     ! Loop index
  integer :: ierr  ! Error code 
  character(len=21) :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Create a Lorenz96 model configured as specified in namelist
  L96 = lorenz96(size, forcing, delta_t)

  ! Read initial model state from restart file if start_t is not 1
  if (start_step /= 1 ) then
    write(restart_file,'(A,I0.6,A)') 'lorenz96out_', start_step, '.nc'
    ierr = L96%nc_read_model_state(trim(restart_file))
  end if

  ! Write out initial model state
  if (output_interval_steps < run_steps) ierr = L96%nc_write_model_state()

  ! Run the model
  do t = 1, run_steps, output_interval_steps

    ! Advance the model to next output interval
    call L96%adv_nsteps(min(output_interval_steps, run_steps))

    ! Write out model state if needed
    if (output_interval_steps < run_steps) ierr = L96%nc_write_model_state()

  end do

end program
