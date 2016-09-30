program sineTL_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use sine, only : sine_type, sine_TL_type
  use kind, only : r8kind, i8kind

  implicit none

  ! Sine model objects for M(t)
  type(sine_type) :: S

  ! SineTL model object to calculate M'(t + delta_t)
  type(sine_TL_type) :: STL 

  ! Define namelists and default values
  integer          :: size    = 40
  real(r8kind)     :: amplitude = 50.0_r8kind
  real(r8kind)     :: bias = 50.0_r8kind
  real(r8kind)     :: frequency = 20.0_r8kind
  real(r8kind)     :: phase = 0.0_r8kind
  real(r8kind)     :: delta_t = 0.05_r8kind
  namelist /params/  size, amplitude, bias, frequency, phase, delta_t
  integer          :: start_step = 0
  integer          :: run_steps   = 1000
  integer          :: output_interval_steps = 100
  character(len=6) :: io_format = 'NETCDF'
  namelist /runtime/ start_step, run_steps, output_interval_steps, io_format

  integer :: step         ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Write column headers
  write(*,'(A12,1x,A12,18x,A5,21x,A16)') "Step", "M(t+dt)", "M'(t)", "M(t+dt) / M'(dt)"

  ! Create a Sine model configured as specified in namelist for M(t)
  S = sine_type(size, amplitude, bias, frequency, phase, delta_t)

  ! Create a SineTL model configured as specified in namelist to calculate dx
  STL = sine_TL_type(size, amplitude, bias, frequency, phase, delta_t)

  ! Read initial model state from previous output file if start_t is not 0
  if (start_step /= 0 ) then
    ierr = S%read_model_state(start_step, io_format)
    ierr = STL%read_model_state(start_step, io_format)
  end if

  ! Loop over digits of precision
  do step = 1, run_steps

    ! Advance the models one step
    call S%adv_nsteps(1)
    call STL%adv_nsteps(1)

    write(*,'(I,3F)')  step, S%state(20), STL%state(20), S%state(20) / STL%state(20)

  end do

end program
