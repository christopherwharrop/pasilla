program sineADJ_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use sine, only : sine_type, sine_TL_type, sine_ADJ_type
  use kind, only     : r8kind, i8kind

  implicit none

  ! Sine model object 
  type(sine_type) :: S

  ! SineTL model object to go forward
  type(sine_TL_type) :: STL 

  ! SineADJ model object to go backward
  type(sine_ADJ_type) :: SADJ

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

  integer :: t            ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Write column headers
  write(*,'(A12,A7,24x,A,17x,A)') "Step", "x", "M*(M'(x))", "x / M*(M'(x))"

  ! Create a Sine model configured as specified in namelist
  S = sine_type(size, amplitude, bias, frequency, phase, delta_t)

  ! Create a SineTL model configured as specified in namelist
  STL = sine_TL_type(size, amplitude, bias, frequency, phase, delta_t)

  ! Create a SineADJ model configured as specified in namelist
  SADJ = sine_ADJ_type(size, amplitude, bias, frequency, phase, delta_t)

  ! Read initial model state from previous output file if start_t is not 0
  if (start_step /= 0 ) then
     ierr = S%read_model_state(start_step, io_format)
     ierr = STL%read_model_state(start_step, io_format)
     ierr = SADJ%read_model_state(start_step, io_format)
  end if

  ! Loop over digits of precision
  do t = 1, run_steps

    ! Initialize tangent linear state to real model state
    STL%state = S%state
    STL%trajectory = S%state

    ! Advance the tangent linear model 1 step to go forward
    call STL%adv_nsteps(1)

    ! Set the adjoint state to the the tangent linear state
    SADJ%state = STL%state
    SADJ%trajectory = STL%state

    ! Advance the adjoint model 1 step to go backward
    call SADJ%adv_nsteps(1)

    ! Print out ratio of norms
    write (*,'(I,4F)') S%step, S%state(20), SADJ%state(20), S%state(20) / SADJ%state(20)

    ! Advance the real model to repeat for next step
    call S%adv_nsteps(1)

  end do

end program
