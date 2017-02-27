program runqgmodel

  use QG
  use gptl

  implicit none

  type(qg_model_type) :: qgmodel

  integer istep, nstep, ret

  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('QG')

  qgmodel = qg_model_type()

  write(*,*) 'Integrating transient days: ', ndayskip

  nstep = ndayskip / dt

  do istep = 1, nstep
    call qgmodel%forward
  enddo

  write(*,*) 'Integrating trajectory of days: ', nday

  istep = 0
  nstep = nday / dt

  call qgmodel%diag(istep)

  do istep = 1, nstep
    call qgmodel%forward
    call qgmodel%diag(istep)
  enddo

  call qgmodel%writestate

  ret = gptlstop ('QG') 
  ret = gptlpr (0) 
  !ret = gptlpr_summary (MPI_COMM_WORLD) 
  ret = gptlfinalize ()

 
end program runqgmodel
