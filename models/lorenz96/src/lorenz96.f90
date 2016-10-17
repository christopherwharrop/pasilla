program lorenz96Model

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use lorenz96, only : lorenz96_type
  use kind, only     : r8kind, i8kind

  implicit none

  ! Lorenz96 model object
  type(lorenz96_type) :: L96

  ! Define namelists and default values
  integer          :: size    = 40
  real(r8kind)     :: forcing = 8.00_r8kind
  real(r8kind)     :: delta_t = 0.05_r8kind
  namelist /params/  size, forcing, delta_t
  integer          :: start_step = 0
  integer          :: run_steps   = 1000
  integer          :: output_interval_steps = 100
  character(len=6) :: io_format = 'NETCDF'
  namelist /runtime/ start_step, run_steps, output_interval_steps, io_format

  integer :: t     ! Loop index
  integer :: ierr  ! Error code 
  character(len=21) :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Create a Lorenz96 model configured as specified in namelist
  L96 = lorenz96_type(size, forcing, delta_t)

  ! Read initial model state from previous output file if start_t is not 0
  if (start_step /= 0 ) then
!    ierr = L96%read_model_state(start_step, io_format)
    L96 = lorenz96_type(start_step, io_format)
  end if

  ! Write out initial model state
  if (output_interval_steps <= run_steps) ierr = L96%write_model_state(io_format)

  ! Run the model
  do t = 1, run_steps, output_interval_steps

    ! Advance the model to next output interval
    call L96%adv_nsteps(min(output_interval_steps, run_steps))

    ! Write out model state if needed
    if (output_interval_steps <= run_steps) ierr = L96%write_model_state(io_format)

  end do

end program
