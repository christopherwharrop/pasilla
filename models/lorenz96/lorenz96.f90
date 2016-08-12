program lorenz96Model

  ! Get unit numbers for stdin, stdout, stderr
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use module_kind, only : r8kind, i8kind
  use module_lorenz96, only : lorenz96  

  implicit none

  ! Lorenz96 model object
  type(lorenz96) :: L96

  ! Define namelists and default values
  integer         :: size    = 40
  real(r8kind)    :: forcing = 8.00_r8kind
  real(r8kind)    :: delta_t = 0.05_r8kind
  integer(i8kind) :: start_t = 0
  integer(i8kind) :: end_t   = 1000
  namelist /params/ size, forcing, delta_t 
  namelist /runtime/ start_t, end_t

  integer :: t     ! Loop index
  integer :: ierr  ! Error code 
  character(len=256) :: restart_file

  ! Read namelists from stdin
  read(stdin,nml=params)
  read(stdin,nml=runtime)

  ! Configure a lorenz96 object
  L96 = lorenz96(size, forcing, delta_t)

  ! Read from restart file if start_t is not 0
  if (start_t /= 0 ) then
    write(restart_file,'(A,I0.6,A)') 'lorenz96out_', start_t, '.nc'
    ierr = L96%nc_read_model_state(trim(restart_file))
  else
    ! Otherwise, write out initial model state
    ierr = L96%nc_write_model_state()
  end if

  ! Run the model
  do t=start_t + 1, end_t
    ! Advance the model one time step
    call L96%adv_1step()
  end do

  ! Write out the final model state
  ierr = L96%nc_write_model_state()

end program
