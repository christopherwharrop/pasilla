program lorenz96TL_Test

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use lorenz96, only : lorenz96_type, lorenz96_TL_type
  use kind, only     : r8kind, i8kind

  implicit none

  ! Lorenz96 model objects for M(x) and M(x + lambda * dx)
  type(lorenz96_type) :: L96, L96delta

  ! Lorenz96TL model object to calculate dx
  type(lorenz96_TL_type) :: L96TL 

  ! Lorenz96TL model object for M'(lambda * dx)
  type(lorenz96_TL_type) :: L96TLdelta

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

  real(r8kind) :: lambda  ! scaling factor
  integer :: d            ! Loop index
  integer :: ierr         ! Error code 
  character(len=21)  :: restart_file
  integer, parameter :: digits = 12  ! Number of significant digits to check

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Initialize lambda
  lambda = 1.0_r8kind

  ! Write column headers
  write(*,'(A12,1x,A12,19x,A)') "Step", "Lambda", "( M(x + lambda * dx) - M(x) ) / M'(lambda * dx)"

  ! Loop over digits of precision
  do d = 1, digits

    ! Create a Lorenz96 model configured as specified in namelist for M(x)
    L96 = lorenz96_type(size, forcing, delta_t)

    ! Create a Lorenz96 model configured as specified in namelist for M(x + lambda * dx)
    L96delta = lorenz96_type(size, forcing, delta_t)

    ! Create a Lorenz96TL model configured as specified in namelist to calculate dx
    L96TL = lorenz96_TL_type(size, forcing, delta_t)

    ! Create a Lorenz96TL model configured as specified in namelist to calculate M'(lambda * dx)
    L96TLdelta = lorenz96_TL_type(size, forcing, delta_t)

    ! Read initial model state from previous output file if start_t is not 0
    if (start_step /= 0 ) then
      ierr = L96%read_model_state(start_step, io_format)
      ierr = L96delta%read_model_state(start_step, io_format)
      ierr = L96TL%read_model_state(start_step, io_format)
      ierr = L96TLdelta%read_model_state(start_step, io_format)
    end if

    ! Advance the L96TL one step to obtain x + dx
    call L96TL%adv_nsteps(1)

    ! Initialize L96delta to x + dx * lambda
    L96delta%state = L96%state + ((L96TL%state - L96%state) * lambda)

    ! Initialize L96TLdelta to dx * lambda 
    L96TLdelta%state = (L96TL%state - L96%state) * lambda

    ! Advance L96, L96delta, and TLdelta
    call L96%adv_nsteps(1)
    call L96delta%adv_nsteps(1)
    call L96TLdelta%adv_nsteps(1)

    ! Calculate and print test metric
    write(*,'(I,2F)') L96%step, lambda, (L96delta%state(20) - L96%state(20)) / L96TLdelta%state(20)

    ! Increase precision
    lambda = lambda / 10.0_r8kind

  end do

end program
