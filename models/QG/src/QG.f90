program QG

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  use kind
  use QG_Config
  use QG_Model
  use gptl

  implicit none

  integer :: start_step = 0
  integer :: spinup_steps = 720
  integer :: run_steps = 1440
  integer :: output_interval_steps = 3
  logical :: readstart = .true.

  namelist /runtime/ start_step, spinup_steps, run_steps, output_interval_steps, readstart

  type(qg_config_type) :: config
  type(qg_model_type) :: model

  integer step, ret

  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('QG')

  ! Read namelist from stdin
  read(stdin,nml=runtime)
  rewind(stdin)

  config = qg_config_type(stdin)
  model = qg_model_type(config)

  write(*,*) 'Integrating transient days: ', spinup_steps

  do step = 1, spinup_steps
    call model%forward
  enddo

  write(*,*) 'Integrating trajectory of days: ', run_steps

  call model%diag(0)

  do step = 1, run_steps
    call model%forward
    call model%diag(step)
  enddo

  call model%writestate

  ret = gptlstop ('QG') 
  ret = gptlpr (0) 
  !ret = gptlpr_summary (MPI_COMM_WORLD) 
  ret = gptlfinalize ()

 
end program QG
