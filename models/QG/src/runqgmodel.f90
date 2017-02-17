program runqgmodel

  use QG
  use gptl

  implicit none

  integer istep, nstep, ret

  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('QG')

  call initqg

  write(*,*) 'Integrating transient days: ', ndayskip

  nstep = ndayskip / dt

  do istep = 1, nstep
    call forward
  enddo

  write(*,*) 'Integrating trajectory of days: ', nday

  istep = 0
  nstep = nday / dt

  call diag(istep)

  do istep = 1, nstep
    call forward
    call diag(istep)
  enddo

  call writestate

  ret = gptlstop ('QG') 
  ret = gptlpr (0) 
  !ret = gptlpr_summary (MPI_COMM_WORLD) 
  ret = gptlfinalize ()

 
end program runqgmodel
