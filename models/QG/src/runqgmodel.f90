!23456789012345678901234567890123456789012345678901234567890123456789012
      program runqgmodel
!-----------------------------------------------------------------------
! *** integrates the qgmodel with parameters from inputdata/namelist
!-----------------------------------------------------------------------
      use QG

      implicit none
      
      integer istep,nstep

      call initqg

      write(*,*) 'Experiment ',expid
      write(*,*) 
      write(*,*) 'Integrating transient days: ',ndayskip
      
      nstep=ndayskip/dt
      
      do istep=1,nstep
        call forward
      enddo
      
      write(*,*) 'Integrating trajectory of days: ',nday
      
      istep=0
      nstep=nday/dt
      
!     call diagsf(istep)
      call diag(istep)
      
      do istep=1,nstep
        call forward
!       call diagsf(istep)
        call diag(istep)
      enddo
      
      call writestate
 
    end program runqgmodel
