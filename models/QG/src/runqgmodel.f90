program runqgmodel

  use QG

  implicit none

  integer istep, nstep

  call initqg

  write(*,*) 'Integrating transient days: ', ndayskip

  nstep = ndayskip / dt

  do istep = 1, nstep
    call forward
  enddo

  write(*,*) 'Integrating trajectory of days: ', nday

  istep = 0
  nstep = nday / dt

! call diagsf(istep)
  call diag(istep)

  do istep = 1, nstep
    call forward
!   call diagsf(istep)
    call diag(istep)
  enddo

  call writestate
 
end program runqgmodel
