program lorenz96TL_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use lorenz96, only : lorenz96_type, lorenz96_TL_type, lorenz96_ADJ_type
  use kind, only     : r8kind, i8kind

  implicit none

  ! Lorenz96 model object 
  type(lorenz96_type) :: L96

  ! Lorenz96TL model object to go forward
  type(lorenz96_TL_type) :: L96TL 

  ! Lorenz96ADJ model object to go backward
  type(lorenz96_ADJ_type) :: L96ADJ

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

  integer :: t            ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Write column headers
  write(*,'(A12,1x,A7,24x,A,16x,A)') "Step", "x", "M*(M'(x))", "x / M*(M'(x))"

  ! Create a Lorenz96 model configured as specified in namelist
  L96 = lorenz96_type(size, forcing, delta_t)

  ! Create a Lorenz96TL model configured as specified in namelist
  L96TL = lorenz96_TL_type(size, forcing, delta_t)

  ! Create a Lorenz96ADJ model configured as specified in namelist
  L96ADJ = lorenz96_ADJ_type(size, forcing, delta_t)

  ! Read initial model state from previous output file if start_t is not 0
  if (start_step /= 0 ) then
     ierr = L96%read_model_state(start_step, io_format)
     ierr = L96TL%read_model_state(start_step, io_format)
     ierr = L96ADJ%read_model_state(start_step, io_format)
  end if

  ! Loop over digits of precision
  do t = 1, run_steps

    ! Initialize tangent linear state to real model state
    L96TL%state = L96%state
    L96TL%trajectory = L96%state

    ! Advance the tangent linear model 1 step to go forward
    call L96TL%adv_nsteps(1)

    ! Set the adjoint state to the the tangent linear state
    L96ADJ%state = L96TL%state
    L96ADJ%trajectory = L96TL%state

    ! Advance the adjoint model 1 step to go backward
    call L96ADJ%adv_nsteps(1)

    ! Print out ratio of norms
    write (*,'(I,4F)') L96%step, L96%state(20), L96ADJ%state(20), L96%state(20) / L96ADJ%state(20)

    ! Advance the real model to repeat for next step
    call L96%adv_nsteps(1)

  end do

end program
