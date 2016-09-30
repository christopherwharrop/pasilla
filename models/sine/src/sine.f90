program sineModel

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use sine, only : sine_type
  use kind, only : r8kind, i8kind

  implicit none

  ! Sine model object
  type(sine_type) :: sine_wave

  ! Define namelists and default values
  integer          :: size    = 40
  real(r8kind)     :: amplitude = 50.0_r8kind
  real(r8kind)     :: bias = 50.0_r8kind
  real(r8kind)     :: frequency = 20.0_r8kind
  real(r8kind)     :: phase = 0.0_r8kind
  real(r8kind)     :: delta_t = 1.0_r8kind
  namelist /params/  size, amplitude, bias, frequency, phase, delta_t
  integer          :: start_step = 0
  integer          :: run_steps   = 3
  integer          :: output_interval_steps = 1
  character(len=6) :: io_format = 'NETCDF'
  namelist /runtime/ start_step, run_steps, output_interval_steps, io_format

  integer :: t     ! Loop index
  integer :: ierr  ! Error code 
  character(len=21) :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Create a Sine model configured as specified in namelist
  sine_wave = sine_type(size, amplitude, bias, frequency, phase, delta_t)

  ! Read initial model state from previous output file if start_t is not 0
  if (start_step /= 0 ) then
    ierr = sine_wave%read_model_state(start_step, io_format)
  end if

  ! Write out initial model state
  if (output_interval_steps <= run_steps) ierr = sine_wave%write_model_state(io_format)

  ! Run the model
  do t = 1, run_steps, output_interval_steps

    ! Advance the model to next output interval
    call sine_wave%adv_nsteps(min(output_interval_steps, run_steps))

    ! Write out model state if needed
    if (output_interval_steps <= run_steps) ierr = sine_wave%write_model_state(io_format)

  end do

end program
