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

  type(qg_config_type) :: config
  type(qg_model_type) :: model

  integer istep, nstep, ret
  real(r8kind) :: dt

  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('QG')

  config = qg_config_type(stdin)
  model = qg_model_type(config)

  write(*,*) 'Integrating transient days: ', config%get_ndayskip()

  dt = 1d0 / real(config%get_nstepsperday())
  nstep = config%get_ndayskip() / dt

  do istep = 1, nstep
    call model%forward
  enddo

  write(*,*) 'Integrating trajectory of days: ', config%get_nday()

  istep = 0
  nstep = config%get_nday() / dt

  call model%diag(istep)

  do istep = 1, nstep
    call model%forward
    call model%diag(istep)
  enddo

  call model%writestate

  ret = gptlstop ('QG') 
  ret = gptlpr (0) 
  !ret = gptlpr_summary (MPI_COMM_WORLD) 
  ret = gptlfinalize ()

 
end program QG
