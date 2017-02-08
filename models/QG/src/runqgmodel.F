c23456789012345678901234567890123456789012345678901234567890123456789012
      program runqgmodel
c-----------------------------------------------------------------------
c *** integrates the qgmodel with parameters from inputdata/namelist
c-----------------------------------------------------------------------
      implicit none
      
#include "truncation.h"
#include "comqg.h"
      
      integer istep,nstep
      
!      rootdir="${qgdir}"
      
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
      
c     call diagsf(istep)
      call diag(istep)
      
      do istep=1,nstep
        call forward
c       call diagsf(istep)
        call diag(istep)
      enddo
      
      call writestate
 
c     return
      end
