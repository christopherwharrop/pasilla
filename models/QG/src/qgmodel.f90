!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine initqg
!-----------------------------------------------------------------------
! *** initialise parameters and operators and read initial state
!-----------------------------------------------------------------------
      use ComQG
      implicit none



      integer :: resolution
      integer i,j,k1,k2,k,l,m,n,ifail,ii,jj,i1,j1,nn
      real*8  pigr4,dis,dif,rll
      real*8, allocatable :: ininag(:,:)
      real*8  r1,a,b,c,d,e,sqn,rsqn
      real*8  rnorm,rh0,dd,dlon
      real*8, allocatable :: agg(:,:), agg1(:,:), agg2(:,:) 
      real*8, allocatable :: fmu(:,:), wsx(:)
      
      namelist /param/ tdis,addisl,addish,trel,tdif,idif,h0, &
     &                 rrdef1,rrdef2
      namelist /control/resolution,nstepsperday,nstepsbetweenoutput, &
     &                  ndayskip,nday,obsfile,expid,inf,obsf,readstart
      
      rootdirl=index(rootdir,' ')-1
      
      resolution=21
      inf=.false.
      obsf=.false.
      readstart=.false.
      expid="0000"
      nstepsperday = 36
      nstepsbetweenoutput = 36
      ndayskip = 0
      nday = 10

! *** real parameters

      tdis=3.0
      addisl=0.5
      addish=0.5
      trel=25.
      tdif=3.0
      idif=4
      h0=3.
      rrdef1=.110
      rrdef2=.070

      
      OPEN(15,FILE='namelist.input',status='old',FORM='FORMATTED')

      read(15, NML = control)
      read(15, NML = param)
      
      close(15)
      
      OPEN(16,FILE='namelist.output')
      write(16, NML = control)
      write(16, NML = param)
      close(16)

      call init_comqg(resolution)

      allocate(ininag(nlat,nlon))
      allocate(agg(nlat,nlon), agg1(nlat,nlon), agg2(nlat,nlon))
      allocate(fmu(nlat,2),wsx(nsh2))

      OPEN(11,FILE='./qgcoefT'//trim(ft)//'.dat',FORM='FORMATTED')
      OPEN(12,FILE='./qgstartT'//trim(ft)//'.dat',FORM='FORMATTED')
      OPEN(13,FILE='./qgbergT'//trim(ft)//'.dat',FORM='FORMATTED')
      if (inf) then
        OPEN(14,FILE='./qgpvforT'//trim(ft)//'.dat',FORM='FORMATTED')
      endif

      do i=0,nm
        read(11,*) nshm(i)
      enddo
      do i=1,nsh
        read(11,*) ll(i)
      enddo
      
      pi=4d0*atan(1d0)
      radius=6.37e+6 
      om=4d0*pi/(24d0*3600d0)

      pigr4=4.d0*pi
      rl1=1.0d0/rrdef1**2
      rl2=1.0d0/rrdef2**2
      relt1=max(0.0d0,rl1/(trel*pigr4))
      relt2=max(0.0d0,rl2/(trel*pigr4))
      dis=max(0.0d0,1.0d0/(tdis*pigr4))
      rll=dble(ll(nsh))
      dif=max(0.0d0,1.0d0/(tdif*pigr4*(rll*(rll+1))**idif))
      
! *** time step of the model: 
! *** dt    : fraction of one day
! *** dtime : in seconds
! *** dtt   : dimensionless
 
      dt     = 1d0/real(nstepsperday)
      dtime  = dt*(24d0*3600d0)
      dtt    = dt*pi*4d0

! *** zonal derivative operator

      k2=0
      do m=0,nm
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          rm(k)=dble(m)
        enddo
      enddo

! *** laplace/helmholtz direct and inverse operators

      do j=0,5
        rinhel(1,j)=0.0d0
      enddo

      diss(1,1)=0.0d0
      diss(1,2)=0.0d0

      do k=2,nsh
        r1=dble(ll(k)*(ll(k)+1))
        a=-r1-3.0d0*rl1
        b=-r1-3.0d0*rl2
        c=-r1-rl1
        d=-r1-rl2
        e=a*d+b*c
        rinhel(k,0)=-r1
        rinhel(k,1)=-1.0d0/r1
        rinhel(k,2)= d/e
        rinhel(k,3)= b/e
        rinhel(k,4)=-c/e
        rinhel(k,5)= a/e
        diss(k,2)=dis*r1
        diss(k,1)=-dif*r1**idif
      enddo

      do j=0,5
        do k=1,nsh
          rinhel(k+nsh,j)=rinhel(k,j)
        enddo
      enddo
      
      do j=1,2
        do k=1,nsh
          diss(k+nsh,j)=diss(k,j)
        enddo
      enddo

! *** legendre associated functions and derivatives

      do k=1,nsh
        do j=1,nlat
          read(11,*) pp(j,k)
        enddo
      enddo
      do k=1,nsh
        do j=1,nlat
          read(11,*) pd(j,k)
        enddo
      enddo
      do k=1,nsh
        do j=1,nlat
          read(11,*) pw(j,k)
        enddo
      enddo

! *** compensation for normalization in nag fft routines

      sqn=sqrt(dble(nlon))
      rsqn=1d0/sqn
      do k=1,nsh
        do i=1,nlat
          pp(i,k)=pp(i,k)*sqn
          pd(i,k)=pd(i,k)*sqn
          pw(i,k)=pw(i,k)*rsqn
        enddo
      enddo

! *** initialization of coefficients for fft


      do j=1,nlon
        do i=1,nlat
          ininag(i,j)=1.0d0
        enddo
      enddo 
      
      ifail=0
      call c06fpf (nlat,nlon,ininag,'i',trigd,wgg,ifail)

      ifail=0
      call c06fqf (nlat,nlon,ininag,'i',trigi,wgg,ifail)

! *** orography and dissipation terms
      
! *** fmu(i,1): sin(phi(i))
! *** fmu(i,2): 1-sin**2(phi(i))      

      rnorm=1.0d0/sqrt(3.0d0*nlon)
      do i=1,nlat
        fmu(i,1)=rnorm*pp(i,2)
        fmu(i,2)=1.d0-fmu(i,1)**2
        sinfi(i)=fmu(i,1)
        phi(i)=asin(sinfi(i))
        cosfi(i)=cos(phi(i))
        phi(i)=180d0*phi(i)/pi
      enddo
      dlon=360d0/real(nlon)
      
! *** height of orography in meters
      
      do i=1,nlon
        do j=1,nlat
          read(13,*) agg1(J,I)
        enddo
      enddo
      
      rh0=max(0.0d0,0.001d0/h0)
      do j=1,nlon
        do i=1,nlat
          agg(i,j)=fmu(i,1)*agg1(i,j)*rh0
!          agg(i,j) = agg1(i,j)*rh0
        enddo
      enddo

              
! *** surface dependent friction

      lgdiss=((addisl.gt.0.0).or.(addish.gt.0.0))

      call ggtosp (agg,orog)
      call ddl (orog,ws)
      call sptogg (ws,dorodl,pp)
      call sptogg (orog,dorodm,pd)

      if (lgdiss) then

        do i=1,nlon
          do j=1,nlat
            read(13,*) agg2(j,i)
          enddo
        enddo
        
        do j=1,nlon
          do i=1,nlat
            agg(i,j)=1.0d0+addisl*agg2(i,j)+ &
     &                addish*(1.0d0-exp(-0.001d0*agg1(i,j)))
          enddo
        enddo

        call ggtosp (agg,ws)
        call ddl (ws,wsx)

        call sptogg (ws,rdiss,pp)
        call sptogg (wsx,ddisdx,pp)
        call sptogg (ws,ddisdy,pd)

        dd=0.5d0*diss(2,2)
        do j=1,nlon
          do i=1,nlat
            ddisdx(i,j)=dd*ddisdx(i,j)/fmu(i,2)
            ddisdy(i,j)=dd*ddisdy(i,j)*fmu(i,2)
          enddo
        enddo

      endif

! *** forcing term

      do l=1,3
        do k=1,nsh2
          for(k,l)=0d0
        enddo
      enddo
 
      if (inf) then
       
        read(14,'(1e12.5)') ((for(k,l),k=1,nsh2),l=1,3)
        
      endif
      
      if (obsf) then
              
        call artiforc
        
      endif

! *** input initial streamfunction

      if (readstart) then
        do l=1,3
          do k=1,nsh2
            read(12,*) psi(k,l)
          enddo
        enddo
      else
        do l=1,3
          do k=1,nsh2
            psi(k,l)=0d0
          enddo
        enddo
      endif
      
!
! *** Potential vorticity and streamfunction fields
!

      call psitoq
             
      close(11)
      close(12)
      close(13)
      close(14)
      
      OPEN(13,FILE='qgbergT'//trim(ft)//'.grads', &
     &        FORM='UNFORMATTED')
      write(13) ((real(agg1(j,i)),i=1,nlon),j=1,nlat)
      write(13) ((real(agg2(j,i)),i=1,nlon),j=1,nlat)
      close(13)
      open(50,file='qgbergT'//trim(ft)//'.ctl', &
     &          form='formatted')
        write(50,'(A)') 'dset ^qgbergT'//trim(ft)//'.grads'
        write(50,'(A)') 'undef 9.99e+10'
        write(50,'(A)') 'options sequential big_endian'
        write(50,'(A)') 'title three level QG model'
        write(50,'(A)') '*'
        write(50,'(A,i4,A,F19.14)') &
     &            'xdef ',nlon,' linear  0.000 ',dlon
        write(50,'(A)') '*'
        write(50,'(A,I4,A,1F19.14)') 'ydef ',nlat,' levels ',phi(1)
        write(50,'(F19.14)') (phi(j),j=2,nlat)
        write(50,'(A)') '*'
        write(50,'(A)') 'zdef  1 levels 1000'
        write(50,'(A)') '*'
        write(50,'(A)') 'tdef 1 linear 1jan0001 1dy'
        write(50,'(A)') '*'
        write(50,'(A)') 'vars  2'
        write(50,'(A)') 'oro    1  99 orography [m]'
        write(50,'(A)') 'friction    1  99 friction mask'
        write(50,'(A)') 'endvars'

      close(50)
      
      OPEN(14,FILE='qgpvforT'//trim(ft)//'.grads', &
     & FORM='UNFORMATTED')
      do l=1,nvl
        call sptogg(for(1,l),agg1,pp)
        write(14) ((real(agg1(j,i)),i=1,nlon),j=1,nlat)
      enddo
      close(14)
      
      open(50,file='qgpvforT'//trim(ft)//'.ctl', &
     &          form='formatted')
        write(50,'(A)') 'dset ^qgpvforT'//trim(ft)//'.grads'
        write(50,'(A)') 'undef 9.99e+10'
        write(50,'(A)') 'options sequential big_endian'
        write(50,'(A)') 'title three level QG model'
        write(50,'(A)') '*'
        write(50,'(A,i4,A,F19.14)') &
     &            'xdef ',nlon,' linear  0.000 ',dlon
        write(50,'(A)') '*'
        write(50,'(A,I4,A,1F19.14)') 'ydef ',nlat,' levels ',phi(1)
        write(50,'(F19.14)') (phi(j),j=2,nlat)
        write(50,'(A)') '*'
        write(50,'(A)') 'zdef  3 levels 800 500 200'
        write(50,'(A)') '*'
        write(50,'(A)') 'tdef 1 linear 1jan0001 1dy'
        write(50,'(A)') '*'
        write(50,'(A)') 'vars  1'
        write(50,'(A)') 'pvfor    3  99 pv forcing field [nondim]'
        write(50,'(A)') 'endvars'

      close(50)
      
      return
      end

!1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
      subroutine addperturb
      
      use ComQG
      implicit none
      


      
      integer ipert,i,j,l
      
      real*8 qpgg(nlat,nlon),ran1
      
      ipert=-1
      
      do l=1,nvl
        call sptogg(qprime(1,l), qpgg,pp)
        do j=1,nlon
          do i=1,nlat
            qpgg(i,j)=qpgg(i,j)*(1.025-0.05*ran1(ipert))
          enddo
        enddo
        call ggtosp(qpgg,qprime(1,l))
      enddo
      
      end
      
!1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
!  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      FUNCTION ran1(idum)
      use ComQG
      implicit none
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
     & NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ddt
!----------------------------------------------------------------------
! *** computation of time derivative of the potential vorticity fields
! *** input qprime, psi, psit
! *** output dqprdt
! *** NOTE psit is destroyed
!----------------------------------------------------------------------

      use ComQG
      implicit none



      integer k,l,i,j
      real*8  dum1,dum2
      
! *** advection of potential vorticity at upper level
 
      call jacob (psi(1,1),qprime(1,1),dqprdt(1,1))
 
! *** advection of potential vorticity at middle level
 
      call jacob (psi(1,2),qprime(1,2),dqprdt(1,2))
 
! *** advection of potential vorticity and dissipation at lower level

      call jacobd (psi(1,3),qprime(1,3),dqprdt(1,3))
 
! *** relaxation of temperature and forcing

      do k=1,nsh2
        dum1=relt1*psit(k,1)
        dum2=relt2*psit(k,2)
        dqprdt(k,1)=dqprdt(k,1)+dum1              +for(k,1)
        dqprdt(k,2)=dqprdt(k,2)-dum1+dum2         +for(k,2)
        dqprdt(k,3)=dqprdt(k,3)          -dum2    +for(k,3)
      enddo
 
! *** explicit horizontal diffusion
 
      do l=1,3
        do k=1,nsh2
          dqprdt(k,l)=dqprdt(k,l)+diss(k,1)*qprime(k,l)
        enddo
      enddo
                  
      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine jacob (psiloc,pvor,sjacob)
!----------------------------------------------------------------------
! *** advection of potential vorticity
! *** input psiloc, pvor
! *** output sjacob
!----------------------------------------------------------------------
      use ComQG
      implicit none



      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon), &
     &        dvordm(nlat,nlon), gjacob(nlat,nlon), dpsidls(nsh2)
 
! *** space derivatives of potential vorticity
 
      call ddl (pvor,vv)
      call sptogg (vv,dvordl,pp)
      call sptogg (pvor,dvordm,pd)
 
! *** space derivatives of streamfunction
 
      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)
 
! *** jacobian term
 
      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dvordl(i,j)-dpsidl(i,j)*dvordm(i,j)
        enddo
      enddo
 
      call ggtosp (gjacob,sjacob)
 
! *** planetary vorticity advection
 
      do k=1,nsh2
        sjacob(k)=sjacob(k)-dpsidls(k)
      enddo
 
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine jacobd (psiloc,pvor,sjacob)
!----------------------------------------------------------------------
! *** advection of potential vorticity and dissipation on gaussian grid
! *** input psiloc, pvor
! *** output sjacob
!----------------------------------------------------------------------
      use ComQG
      implicit none



      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon), &
     &        dvordm(nlat,nlon), gjacob(nlat,nlon), vv(nsh2), &
     &        azeta(nlat,nlon),dpsidls(nsh2)
 
! *** space derivatives of potential vorticity
 
      call ddl (pvor,vv)
      call sptogg (vv,dvordl,pp)
      call sptogg (pvor,dvordm,pd)
 
! *** space derivatives of streamfunction
 
      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)
 
! *** jacobian term + orographic forcing
 
      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*(dvordl(i,j)+sinfi(i)*dorodl(i,j))- &
     &                dpsidl(i,j)*(dvordm(i,j)+sinfi(i)*dorodm(i,j))
        enddo
      enddo

! *** dissipation 
 
 
      if (lgdiss) then

! ***   spatially varying dissipation 

        do k=1,nsh2
          vv(k)=diss(k,2)*psiloc(k)
        enddo
 
        call sptogg (vv,azeta,pp)
 
        do j=1,nlon
          do i=1,nlat
            gjacob(i,j)=gjacob(i,j) - dpsidm(i,j)*ddisdy(i,j) &
     &              -dpsidl(i,j)*ddisdx(i,j) &
     &              +rdiss(i,j)*azeta(i,j)       
          enddo
        enddo

        call ggtosp (gjacob,sjacob)

      else

! ***   uniform dissipation

        call ggtosp (gjacob,sjacob)

        do k=1,nsh2
          sjacob(k)=sjacob(k)+diss(k,2)*psi(k,3)
        enddo

      endif
  
! *** planetary vorticity advection
 
      do k=1,nsh2
        sjacob(k)=sjacob(k)-dpsidls(k)
      enddo
 
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ddl (as,dadl)
!-----------------------------------------------------------------------
! *** zonal derivative in spectral space
! *** input spectral field as
! *** output spectral field dadl which is as differentiated wrt lambda
!-----------------------------------------------------------------------

      use ComQG
      implicit none



      integer k
      real*8 as(nsh,2), dadl(nsh,2)
 
      do k=1,nsh
        dadl(k,1)=-rm(k)*as(k,2)
        dadl(k,2)= rm(k)*as(k,1)
      enddo
 
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sptogg (as,agg,pploc)
 
!-----------------------------------------------------------------------
! *** conversion from spectral coefficients to gaussian grid
! *** input  spectral field as, legendre polynomials pploc (pp or pd) 
! ***        where pp are legendre polynomials and pd derivatives with
! ***        respect to sin(fi)
! *** output gaussian grid agg
!-----------------------------------------------------------------------
 
      use ComQG
      implicit none


      
      integer i,ifail,j,k,k1,k2,m,mi,mr,nlon1
      real*8  as(nsh,2), agg(nlat,nlon), pploc(nlat,nsh)
 
! *** inverse legendre transform
 
      do j=1,nlon
        do i=1,nlat
          agg(i,j)=0.0d0
        enddo
      enddo
 
      nlon1=nlon+1
      k2=nshm(0)
 
      do k=1,k2
        do i=1,nlat
          agg(i,1)=agg(i,1)+as(k,1)*pploc(i,k)
        enddo
      enddo
 
      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            agg(i,mr)=agg(i,mr)+as(k,1)*pploc(i,k)
          enddo
          do i=1,nlat
            agg(i,mi)=agg(i,mi)-as(k,2)*pploc(i,k)
          enddo
        enddo
      enddo
 
! *** inverse fourier transform
 
      ifail=0
      call c06fqf (nlat,nlon,agg,'r',trigi,wgg,ifail)
 
      return
      end
 
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ggtosp (agg,as)
!-----------------------------------------------------------------------
! *** conversion from gaussian grid (agg) to spectral coefficients (as)
! *** input array agg is destroyed
! *** output as contains spectral coefficients
!-----------------------------------------------------------------------

      use ComQG
      implicit none



      integer ir,ifail,j,k,k1,k2,m,mi,mr,nlon1,i
      real*8  as(nsh,2), agg(nlat,nlon)
!
! *** fourier transform
!
      ifail=0
      call c06fpf (nlat,nlon,agg,'r',trigd,wgg,ifail)
!
! *** legendre transform
!
      do ir=1,2
        do k=1,nsh
          as(k,ir)=0.0d0
        enddo
      enddo
 
      nlon1=nlon+1

      k2=nshm(0)
 
      do k=1,k2
        do i=1,nlat
          as(k,1)=as(k,1)+agg(i,1)*pw(i,k)
        enddo
      enddo
 
      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            as(k,1)=as(k,1)+agg(i,mr)*pw(i,k)
            as(k,2)=as(k,2)+agg(i,mi)*pw(i,k)
          enddo
        enddo
      enddo
 
      return
      end
 

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine qtopsi
!-----------------------------------------------------------------------
! *** computation of streamfunction from potential vorticity
! *** input  qprime which is potential vorticity field
! *** output psi, the streamfunction and psit, the layer thicknesses
!-----------------------------------------------------------------------
 
      use ComQG
      implicit none



      integer k
      real*8  r3

      do k=1,nsh2
        ws(k)=qprime(k,1)+qprime(k,3)
        psi(k,1)=rinhel(k,1)*(ws(k)+qprime(k,2))
        psi(k,2)=ws(k)-2.d0*qprime(k,2)
        psi(k,3)=qprime(k,1)-qprime(k,3)
      enddo
 
      do k=1,nsh2
        psit(k,1)=rinhel(k,2)*psi(k,2)+rinhel(k,3)*psi(k,3)
        psit(k,2)=rinhel(k,4)*psi(k,2)+rinhel(k,5)*psi(k,3)
      enddo
 
      r3=1./3.
      do k=1,nsh2
        psi(k,2)=r3*(psi(k,1)-psit(k,1)+psit(k,2))
        psi(k,1)=psi(k,2)+psit(k,1)
        psi(k,3)=psi(k,2)-psit(k,2)
      enddo
 
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine psitoq 
!-----------------------------------------------------------------------
! *** computation of potential vorticity from stream function
! *** input psi streamfunction
! *** output qprime, the potential vorticity and psit, the layer thick.
!-----------------------------------------------------------------------

      use ComQG
      implicit none


      integer k
       
      do k=1,nsh2
        psit(k,1)=psi(k,1)-psi(k,2)
        psit(k,2)=psi(k,2)-psi(k,3)
        qprime(k,1)=rinhel(k,0)*psi(k,1)-rl1*psit(k,1)
        qprime(k,2)=rinhel(k,0)*psi(k,2)+rl1*psit(k,1)-rl2*psit(k,2)
        qprime(k,3)=rinhel(k,0)*psi(k,3)+rl2*psit(k,2)
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine psiq(sfin,qout)
!-----------------------------------------------------------------------
! ***  computation of potential vorticity qout from stream function sfin
!-----------------------------------------------------------------------

      use ComQG
      implicit none


      integer k
      real*8  sfin(nsh2,nvl),qout(nsh2,nvl),tus(nsh2)

      do k=1,nsh2
        tus(k)=rl1*sfin(k,1)-rl1*sfin(k,2)
      enddo

      do k=1,nsh2
        qout(k,1)=rinhel(k,0)*sfin(k,1)-tus(k)
        qout(k,2)=rinhel(k,0)*sfin(k,2)+tus(k)
      enddo

      do k=1,nsh2
        tus(k)=rl2*sfin(k,2)-rl2*sfin(k,3)
      enddo

      do k=1,nsh2
        qout(k,2)=qout(k,2)-tus(k)
        qout(k,3)=rinhel(k,0)*sfin(k,3)+tus(k)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine qpsi(qin,sfout)
!-----------------------------------------------------------------------
! *** computation of streamfunction bb from potential vorticity qin
!-----------------------------------------------------------------------

      use ComQG
      implicit none

     
      real*8  qin(nsh2,nvl),sfout(nsh2,nvl), tus(nsh2,ntl), r3
      integer k

      do k=1,nsh2
        ws(k)=qin(k,1)+qin(k,3)
        sfout(k,1)=rinhel(k,1)*(ws(k)+qin(k,2))
        sfout(k,2)=ws(k)-2.*qin(k,2)
        sfout(k,3)=qin(k,1)-qin(k,3)
      enddo

      do k=1,nsh2
        tus(k,1)=rinhel(k,2)*sfout(k,2)+rinhel(k,3)*sfout(k,3)
        tus(k,2)=rinhel(k,4)*sfout(k,2)+rinhel(k,5)*sfout(k,3)
      enddo

      r3=1./3
      do k=1,nsh2
        sfout(k,2)=r3*(sfout(k,1)-tus(k,1)+tus(k,2))
        sfout(k,1)=sfout(k,2)+tus(k,1)
        sfout(k,3)=sfout(k,2)-tus(k,2)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine qpsit(qin,tus)
!-----------------------------------------------------------------------
! *** computation of thickness tus from potential vorticity qin
!-----------------------------------------------------------------------

      use ComQG
      implicit none

     
      real*8  qin(nsh2,nvl),tus(nsh2,ntl), r3,sfout(nsh2,nvl)
      integer k

      do k=1,nsh2
        ws(k)=qin(k,1)+qin(k,3)
        sfout(k,1)=rinhel(k,1)*(ws(k)+qin(k,2))
        sfout(k,2)=ws(k)-2.*qin(k,2)
        sfout(k,3)=qin(k,1)-qin(k,3)
      enddo

      do k=1,nsh2
        tus(k,1)=rinhel(k,2)*sfout(k,2)+rinhel(k,3)*sfout(k,3)
        tus(k,2)=rinhel(k,4)*sfout(k,2)+rinhel(k,5)*sfout(k,3)
      enddo

      return
      end

 
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fmtofs (y,z)
!-----------------------------------------------------------------------
! *** transforms francos format to the french format for global fields
! *** input  y spectral coefficients in francos format
! *** output z spectral coefficients in french format
! *** fm format:
! *** k       m  n
! *** 1       0  0
! *** 2       0  1
! *** 3       0  2
! *** :       :  :
! *** nm+1    0  nm
! *** nm+2    1  1 --> real part
! *** nm+3    1  2 --> real part
! *** :       :  :
! *** nm+nm+1 1  nm --> real part
! *** :       :  :
! *** :       nm nm --> real part
! ***  repeat for imaginary part
! ***  disadvantage: 0 0 mode and imaginary parts of m=0 modes are obsolete
! *** fs format stores all m for every n first and has no obsolete indices
! *** 
! *** k       m  n
! *** 1       0  1
! *** 2       1  1 --> real part
! *** 3       1  1 --> imaginary part: k=1-3 is T1 truncation
! *** 4       0  2
! *** 5       1  2 --> real part
! *** 6       1  2 --> imaginary part
! *** 7       2  2 --> real part
! *** 8       2  2 --> imaginary part: k=1-8 is T2 truncation
! *** etcetera
!-----------------------------------------------------------------------

      use ComQG
      implicit none



      integer   m,n,k,indx,l
      real*8    y(nsh2,nvl),z(nsh2,nvl)

      do l=1,nvl
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if (m.eq.0) then
              indx=n**2
            else
              indx=n**2+2*m-1
            end if
            z(indx,l)=y(k,l)
            if (m.ne.0) z(indx+1,l)=y(k+nsh,l)
          enddo
        enddo
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine fstofm (y,z,ntr)
!-----------------------------------------------------------------------
! *** transforms the french format to francos format for global fields
! *** input  y spectral coef. in french format, ntr is truncation limit
! *** output z spectral coefficients in francos format
! *** fm format:
! *** k       m  n
! *** 1       0  0
! *** 2       0  1
! *** 3       0  2
! *** :       :  :
! *** nm+1    0  nm
! *** nm+2    1  1 --> real part
! *** nm+3    1  2 --> real part
! *** :       :  :
! *** nm+nm+1 1  nm --> real part
! *** :       :  :
! *** :       nm nm --> real part
! ***  repeat for imaginary part
! ***  disadvantage: 0 0 mode and imaginary parts of m=0 modes are obsolete
! *** fs format stores all m for every n first and has no obsolete indices
! *** 
! *** k       m  n
! *** 1       0  1
! *** 2       1  1 --> real part
! *** 3       1  1 --> imaginary part: k=1-3 is T1 truncation
! *** 4       0  2
! *** 5       1  2 --> real part
! *** 6       1  2 --> imaginary part
! *** 7       2  2 --> real part
! *** 8       2  2 --> imaginary part: k=1-8 is T2 truncation
! *** etcetera
!-----------------------------------------------------------------------

      use ComQG
      implicit none



      integer   m,n,k,indx,i,l,ntr
      real*8    y(nsh2,nvl),z(nsh2,nvl)

      do l=1,nvl
        do i=1,nsh2
          z(i,l)=0d0
        enddo
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if ((m.le.ntr).and.(n.le.ntr)) then
              if (m.eq.0) then
                indx=n**2
              else
                indx=n**2+2*m-1
              end if
              z(k,l)=y(indx,l)
              if (m.ne.0) z(k+nsh,l)=y(indx+1,l)
            endif
          enddo
        enddo
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine truncate(y,yt,ntr)
!-----------------------------------------------------------------------
! *** truncates y to ntr and writes to z which is formatted for
! *** lower resolution model Tntr.
!-----------------------------------------------------------------------

      use ComQG
      implicit none



      integer   m,n,k,indx,l,ntr,nshntr,i
      real*8    y(nsh2,nvl),z(nsh2,nvl),yt(nsh2,nvl)

      nshntr=(ntr+1)*(ntr+2)*0.5
      
      do l=1,nvl
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if (m.eq.0) then
              indx=n**2
            else
              indx=n**2+2*m-1
            end if
            z(indx,l)=y(k,l)
            if (m.ne.0) z(indx+1,l)=y(k+nsh,l)
          enddo
        enddo
      enddo
      
      do l=1,nvl
        do i=1,nsh2
          yt(i,l)=0d0
        enddo
        k=1
        do m=0,ntr
          do n=max(m,1),ntr
            k=k+1
              if (m.eq.0) then
                indx=n**2
              else
                indx=n**2+2*m-1
              end if
              yt(k,l)=z(indx,l)
              if (m.ne.0) yt(k+nshntr,l)=z(indx+1,l)
          enddo
        enddo
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine forward
!-----------------------------------------------------------------------
! *** performs a fourth order runge kutta time step at truncation nm
! *** with time step dt
! *** dqdt calculates the time derivative
! *** input  qprime at current time
! *** output qprime at current time plus dt
!-----------------------------------------------------------------------
      use ComQG
      implicit none


      integer  k,l,nvar
      real*8   dt2,dt6
      real*8   y(nsh2,nvl),dydt(nsh2,nvl),yt(nsh2,nvl)
      real*8   dyt(nsh2,nvl),dym(nsh2,nvl)

      nvar=(nm+2)*nm
      dt2=dtt*0.5d0
      dt6=dtt/6d0
      call fmtofs(qprime,y)
      call dqdt(y,dydt)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dt2*dydt(k,l)
        enddo
      enddo
      call dqdt(yt,dyt)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dt2*dyt(k,l)
        enddo
      enddo
      call dqdt(yt,dym)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dtt*dym(k,l)
          dym(k,l)=dyt(k,l)+dym(k,l)
        enddo
      enddo
      call dqdt(yt,dyt)
      do l=1,nvl
        do k=1,nvar
          y(k,l)=y(k,l)+dt6*(dydt(k,l)+dyt(k,l)+2.*dym(k,l))
        enddo
      enddo
      call fstofm(y,qprime,nm)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine dqdt(y,dydt)
!-----------------------------------------------------------------------
! *** computation of time derivative of the potential vorticity field
! *** input  y potential vorticity in french format
! *** output dydt time derivative of y in french format
! *** values of qprime, psi and psit are changed
!-----------------------------------------------------------------------

      use ComQG
      implicit none


      real*8  y(nsh2,nvl),dydt(nsh2,nvl)

      call fstofm(y,qprime,nm)
      call qtopsi
      call ddt
      call fmtofs(dqprdt,dydt)      
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine gridfields
!-----------------------------------------------------------------------
! *** computation of geostrophic winds at all levels
! *** computes geopotential height in [m2/s2[=[gz] from streamfunction 
! *** by solving the linear balance equation: 
! *** del phi = (1 - mu**2 ) d psi/dmu + mu del psi
! *** the global mean value is not determined and set to zero
!  
!-----------------------------------------------------------------------
      use ComQG
      implicit none



      integer i,j,k,l
      real*8  facwind,facsf,facgp,facpv
      real*8  dpsdl(nlat,nlon),dpsdm(nlat,nlon),psik(nsh2),vv(nsh2)
      real*8  fmu(nlat)
      real*8  delpsis(nsh2),delpsig(nlat,nlon)
      real*8  dmupsig(nlat,nlon),delgeog(nlat,nlon)
      real*8  delgeos(nsh2),geos(nsh2)


      call qtopsi      
      
! *** space derivatives of streamfunction
 
      facwind=radius*om
      facsf=om*(radius)**2
      facgp=(om**2)*(radius**2)
      facpv=om
      
      do i=1,nlat
        fmu(i)=1-sin(pi*phi(i)/180d0)**2
      enddo

      do l=1,nvl
        call sptogg(psi(1,l),psig(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            psig(i,j,l)=facsf*psig(i,j,l)
          enddo
        enddo

        call sptogg(qprime(1,l),qgpv(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            qgpv(i,j,l)=facpv*qgpv(i,j,l)
          enddo
        enddo

        do k=1,nsh2
          psik(k)=psi(k,l)
        enddo


        call ddl (psik,vv)
        call sptogg (vv,dpsdl,pp)
        call sptogg (psik,dpsdm,pd)

        do j=1,nlon
          do i=1,nlat
            ug(i,j,l)=-facwind*dpsdm(i,j)*cosfi(i)
            vg(i,j,l)=+facwind*dpsdl(i,j)/cosfi(i)
          enddo
        enddo
        
! *** solve linear balance equation


        call lap(psi(1,l),delpsis)
        call sptogg(delpsis,delpsig,pp)
        call sptogg(psi(1,l),dmupsig,pd)

        do j=1,nlon
          do i=1,nlat
            delgeog(i,j)=fmu(i)*dmupsig(i,j)+ &
     &                      sinfi(i)*delpsig(i,j)
          enddo
        enddo
        call ggtosp(delgeog,delgeos)
        call lapinv(delgeos,geos)
        geos(1)=0.d0
        call sptogg(geos,geopg(1,1,l),pp)
        
        
        do j=1,nlon
          do i=1,nlat
            geopg(i,j,l)=facgp*geopg(i,j,l)
          enddo
        enddo

      enddo

      return
      end

      
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lap(xs,xsl)
!-----------------------------------------------------------------------
! *** computation of laplace operator in spectral domain
! *** input  xs  field in spectral form
! *** output xsl laplace of xs in spectral form
!-----------------------------------------------------------------------
      use ComQG
      implicit none


      integer k
      real*8  xs(nsh2),xsl(nsh2)

      do k=1,nsh2
        xsl(k)=xs(k)*rinhel(k,0)
      enddo

      return
      end

      

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine lapinv(xsl,xs)
!-----------------------------------------------------------------------
! *** computation of laplace operator in spectral domain
! *** input  xsl field in spectral form
! *** output xs  inverse laplace of xs in spectral form
!-----------------------------------------------------------------------
      use ComQG
      implicit none


      integer k
      real*8  xs(nsh2),xsl(nsh2)

      do k=1,nsh2
        xs(k)=xsl(k)*rinhel(k,1)
      enddo

      return
      end

!-----------------------------------------------------------------------
      subroutine artiforc_iter
!-----------------------------------------------------------------------
! computation of artifical forcing according to roads(1987)
! the forcing is computed from daily ecmwf data for the winter season
! the forcing is composed of a climatological forcing and a contribution
! of the eddy terms.
!------------------------------------------------------------------------

      use ComQG
      implicit none
      


      
      integer i,j,k,l,iy,id,iday,nyb,nye,nd
      real*4 psi4(nsh2,3)
      real*8 psim(nsh2,3),sum(nsh2,3)
      real*8 eddf(nsh2,3),totf(nsh2,3),climf(nsh2,3),forg(nlat,nlon)
      
      
      open(unit=46,file='./qgmodelT42.T21',status='old', &
     &       form='formatted')
      open(unit=32,file='./qgpvforT42world.grads', &
     &       form='unformatted')
      
      nyb=83
      nye=92
      nday=91

      do l=1,3
        do k=1,nsh2
          sum(k,l)=0.
        enddo
      enddo
      
      iday=0

      do iy=nyb,nye

        nd=90
        if (iy .eq. 84 .or. iy .eq. 88 .or. iy .eq. 92) nd=91
         
        do id=1,nd
          read(46,*)i
          read(46,'(5e12.5)')((psi4(k,l),k=1,nsh2),l=1,3)
          iday=iday+1
          do l=1,3
            do k=1,nsh2
              sum(k,l)=sum(k,l)+psi4(k,l)
            enddo
          enddo
        enddo
      enddo
      
      
      do l=1,3
        do k=1,nsh2
          psim(k,l)=sum(k,l)/iday
        enddo
      enddo

! *** calculate the climatological forcing

      do l=1,3
        do k=1,nsh2
          psi(k,l)=psim(k,l)
        enddo
      enddo
      
      call psitoq

      call ddt 
      
      do l=1,3
        do k=1,nsh2
          climf(k,l)=-dqprdt(k,l)
        enddo
      enddo

! *** calculate the eddy forcing
      rewind(46)
      
      do l=1,3
        do k=1,nsh2
          sum(k,l)=0.
        enddo
      enddo
      iday=0

      do iy=nyb,nye

        nd=90
        if (iy .eq. 84 .or. iy .eq. 88 .or. iy .eq. 92) nd=91
         
        do id=1,nday
          read(46,*)i
          read(46,'(5e12.5)')((psi4(k,l),k=1,nsh2),l=1,3)
          iday=iday+1
          do l=1,3
            do k=1,nsh2
              psi(k,l)=psi4(k,l)-psim(k,l)
            enddo
          enddo
          call psitoq
          call eddforc
          do l=1,3
            do k=1,nsh2
              sum(k,l)=sum(k,l)+dqprdt(k,l)
            enddo
          enddo
        enddo
      enddo
      
      do l=1,3
        do k=1,nsh2
          eddf(k,l)=-sum(k,l)/iday
        enddo
      enddo
      
!***  compute the total forcing

      do l=1,3
        do k=1,nsh2
          for(k,l)=climf(k,l)+eddf(k,l)
        enddo
      enddo
         
      write(14,'(1E12.5)') ((for(k,l),k=1,nsh2),l=1,3)
      
      do l=1,3
        call sptogg(for(1,l),forg,pp)
        write(32) ((real(forg(i,j)),j=1,nlon),i=1,nlat)
      enddo

      close(46)
      close(32)
 
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine eddforc
!-----------------------------------------------------------------------
!***  computation of the eddy forcing
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      




      call jacobedd (psi(1,1),qprime(1,1),dqprdt(1,1))

      call jacobedd (psi(1,2),qprime(1,2),dqprdt(1,2))

      call jacobedd (psi(1,3),qprime(1,3),dqprdt(1,3))

      return
      end
      
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine jacobedd (psiloc,pvor,sjacob)
!-----------------------------------------------------------------------
! *** advection of potential vorticity
! *** input psiloc, pvor
! *** output sjacob
! *** the only difference with the routine jacob 
! *** is that the planetary vorticity advection is omitted.
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      



      integer i,j,k

      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2), vv(nsh2)
      real*8  dpsidls(nsh2)
      real*8  dpsidl(ngp),dpsidm(ngp),dvordl(ngp),dvordm(ngp)
      real*8  gjacob(ngp)


! *** space derivatives of potential vorticity

      call ddl (pvor,vv)
      call sptogg (vv,dvordl,pp)
      call sptogg (pvor,dvordm,pd)

! *** space derivatives of streamfunction

      call ddl (psiloc,dpsidls)
      call sptogg (dpsidls,dpsidl,pp)
      call sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,ngp
          gjacob(j)=dpsidm(j)*dvordl(j)-dpsidl(j)*dvordm(j)
      enddo

      call ggtosp (gjacob,sjacob)

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine artiforc
!-----------------------------------------------------------------------
! computation of artifical forcing according to roads(1987)
! the forcing is computed from file obsfile
!------------------------------------------------------------------------

      use ComQG
      implicit none
      


      
      integer i,j,k,l,iday,fl,nvar
      real*4 psi4(nsh2,3)
      real*8 sum(nsh2,3),forg(nlat,nlon),dlon,psifs(nsh2,3)
      
      nvar=(nm+2)*nm

      dlon=360d0/real(nlon)
      
      OPEN(14,FILE='qgpvforT'//trim(ft)//'.dat', &
     &     FORM='FORMATTED')
      
      open(unit=46,file='./'//obsfile, &
     &     status='old',form='unformatted')
      fl=index(obsfile," ")-1
      open(unit=32, &
     &     file='qgpvforT'//trim(ft)//'.grads', &
     &       form='unformatted')
     
      open(unit=99,file='./'//obsfile(1:fl)//trim(ft)//'.grads', &
     &  form='unformatted')
     
      write(*,'(A,A)') "Calculating forcing from ",obsfile
      
      do l=1,nvl
        do k=1,nsh2
          sum(k,l)=0.
          for(k,l)=0d0
        enddo
      enddo
      

! *** calculate the mean tendency

      iday=0

 10   continue
        do l=1,nvl
          read(46,end=20)(psi4(k,nvl-l+1),k=1,nvar)
        enddo
        
        iday=iday+1
        do l=nvl,1,-1
          do k=1,nvar
            psifs(k,l)=psi4(k,l)
          enddo
        enddo
        call fstofm(psifs,psi,nm)
        do l=nvl,1,-1
          call sptogg(psi(1,l),forg,pp)
          write(99)((real(forg(j,i)),i=1,nlon),j=1,nlat)
        enddo
        call psitoq
        call ddt
        do l=1,nvl
          do k=1,nsh2
            sum(k,l)=sum(k,l)+dqprdt(k,l)
          enddo
        enddo
        goto 10
 20   continue
      close(99)
      
      do l=1,nvl
        do k=1,nsh2
          for(k,l)=-sum(k,l)/real(iday)
        enddo
      enddo
               
      write(14,'(1E12.5)') ((for(k,l),k=1,nsh2),l=1,nvl)
      
      do l=nvl,1,-1
        call sptogg(for(1,l),forg,pp)
        write(32) ((real(forg(i,j)),j=1,nlon),i=1,nlat)
      enddo

      close(46)
      close(32)
      
      open(50,file='./'//obsfile(1:fl)//trim(ft)//'.ctl', &
     &          form='formatted')
        write(50,'(A)') 'dset ^'//obsfile(1:fl)//trim(ft)//'.grads'
        write(50,'(A)') 'undef 9.99e+10'
        write(50,'(A)') 'options sequential big_endian'
        write(50,'(A)') 'title three level QG model'
        write(50,'(A)') '*'
        write(50,'(A,i4,A,F19.14)') &
     &            'xdef ',nlon,' linear  0.000 ',dlon
        write(50,'(A)') '*'
        write(50,'(A,I4,A,1F19.14)') 'ydef ',nlat,' levels ',phi(1)
        write(50,'(F19.14)') (phi(j),j=2,nlat)
        write(50,'(A)') '*'
        write(50,'(A)') 'zdef  3 levels 800 500 200'
        write(50,'(A)') '*'
        write(50,'(A,I6,A)') 'tdef ',iday,' linear 1jan0001 1dy'
        write(50,'(A)') '*'
        write(50,'(A)') 'vars  1'
        write(50,'(A)') 'sf    3  99 streamfunction (nondim)'
        write(50,'(A)') 'endvars'

      close(50)
      
      write(*,'(A,I6)') &
     &  "Number of states used to calculate forcing: ",iday
      write(*,'(A)') "Forcing saved in files: "
      write(*,'(A)') 'qgpvforT'//trim(ft)//'.grads'
      write(*,'(A)') 'qgpvforT'//trim(ft)//'.dat'
 
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine analyses
!-----------------------------------------------------------------------
! *** convert ecmwf winter analyses files to asc file to be read by
! *** artiforc
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      



      integer i,j,k,nyb,nye,npw,iy,id,l,iday,irec,nsh2ntr

      real*8  psiloc(nsh2,3),psiT21(nsh2,3)
      real*4  psi4(nsh2,3),psig4(nlat,nlon,nvl),scalesf
      character*4 fy
      
      scalesf=1d0/(radius*radius*om)
      nsh2ntr=22*23
      
      open(98,file='./anwin79_10T42.dat',form='formatted')
      open(96,file='./anwin79_10T21.dat',form='formatted')
!      open(unit=97,file='./eracheck.grads',
!     *  form='unformatted')
      
      nyb=1979
      nye=2010
      npw=90
      
      write(fy,'(I4.4)') iy

      do iy=nyb,nye
        write(fy,'(I4.4)') iy
        open(unit=99,file='./era'//fy//'djf.grads', &
     &  form='unformatted',access='direct',recl=nlat*nlon*4)

        do id=1,npw
          do l=nvl,1,-1
            irec=(id-1)*3+nvl-l+1
            read(99,rec=irec) ((psig4(j,i,l),i=1,nlon),j=1,nlat)
          enddo
          do l=nvl,1,-1
            do i=1,nlon
              do j=1,nlat
                psig(j,i,l)=psig4(j,i,l)*scalesf
              enddo
            enddo
!            write(97) ((real(psig(j,i,l)),i=1,nlon),j=1,nlat)
          enddo
          do l=1,nvl
            call ggtosp(psig(1,1,l),psiloc(1,l))
          enddo
          write(98,*) id
          write(98,'(5e12.5)')((real(psiloc(k,l)),k=1,nsh2),l=1,nvl)
          call truncate(psiloc,psiT21,21)
          write(96,*) id
          write(96,'(5e12.5)')((real(psiT21(k,l)),k=1,nsh2ntr),l=1,nvl)
        enddo        
        close(99)
      enddo
      close(98)
!      close(97)
      
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine diagsf(istep)
!-----------------------------------------------------------------------
! *** output streamfunction data to outputfile
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      



      integer istep,i,j,k,nout

      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2),dlon
      real*4  psi4(nsh2,3)
      
      nout=nday*nstepsperday/nstepsbetweenoutput+1
      
      dlon=360d0/real(nlon)
      
      if (istep.eq.0) then
        open(52, &
     &  file='qgmodelsfT'//trim(ft)//'.ctl', &
     &          form='formatted')
        write(52,'(A)') 'dset ^qgmodelsfT'//trim(ft)//'.grads'
        write(52,'(A)') 'undef 9.99e+10'
        write(52,'(A)') 'options sequential big_endian'
        write(52,'(A)') 'title T'//trim(ft)//' QG model exp '//expid
        write(52,'(A)') '*'
        write(52,'(A,I6,A,F19.14)')  &
     &                       'xdef ',nlon,' linear  0.000 ',dlon
        write(52,'(A)') '*'
        write(52,'(A,I4,A,1F19.14)') 'ydef ',nlat,' levels ',phi(1)
        write(52,'(F19.14)') (phi(j),j=2,nlat)
        write(52,'(A)') '*'
        write(52,'(A)') 'zdef  3 levels 800 500 200'
        write(52,'(A)') '*'
        write(52,'(A,I6,A)') 'tdef ',nout,' linear 1jan0001 1dy'
        write(52,'(A)') '*'
        write(52,'(A)') 'vars  1'
        write(52,'(A)') 'psi    3  99 streamfunction [m2/s]'
        write(52,'(A)') 'endvars'

        close(52)
        open(52, &
     &       file='qgmodelsfT'//trim(ft)//'.grads', &
     &          form='unformatted')
      endif
      
      if (mod(istep,nstepsbetweenoutput).eq.0) then
        call gridfields
        do k=nvl,1,-1
          write(52) ((real(psig(j,i,k)),i=1,nlon),j=1,nlat)
        enddo
      endif
      
      end
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine diag(istep)
!-----------------------------------------------------------------------
! *** output model data to outputfile
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      



      integer istep,i,j,k,nout

      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2), dlon
      real*4  psi4(nsh2,3)
      
      nout=nday*nstepsperday/nstepsbetweenoutput+1
      dlon=360d0/real(nlon)
      
      if (istep.eq.0) then
        open(50,file='qgmodelT'//trim(ft)//'.ctl', &
     &          form='formatted')
        write(50,'(A)') 'dset ^qgmodelT'//trim(ft)//'.grads'
        write(50,'(A)') 'undef 9.99e+10'
        write(50,'(A)') 'options sequential big_endian'
        write(50,'(A)') 'title T'//trim(ft)//' QG model exp '//expid
        write(50,'(A)') '*'
        write(50,'(A,i4,A,F19.14)') &
     &            'xdef ',nlon,' linear  0.000 ',dlon
        write(50,'(A)') '*'
        write(50,'(A,I4,A,1F19.14)') 'ydef ',nlat,' levels ',phi(1)
        write(50,'(F19.14)') (phi(j),j=2,nlat)
        write(50,'(A)') '*'
        write(50,'(A)') 'zdef  3 levels 800 500 200'
        write(50,'(A)') '*'
        write(50,'(A,I6,A)') 'tdef ',nout,' linear 1jan0001 1dy'
        write(50,'(A)') '*'
        write(50,'(A)') 'vars   5'
        write(50,'(A)') 'geopg  3  99 geopotential [m2/s2]'
        write(50,'(A)') 'psi    3  99 streamfunction [m2/s]'
        write(50,'(A)') 'pv     3  99 pv [1/s]'
        write(50,'(A)') 'u      3  99 zonal velocity [m/s]'
        write(50,'(A)') 'v      3  99 meridional velocity [m/s]'
        write(50,'(A)') 'endvars'

        close(50)
        open(50, &
     &  file='qgmodelT'//trim(ft)//'.grads', &
     &          form='unformatted')
      endif
      
      if (mod(istep,nstepsbetweenoutput).eq.0) then
        call gridfields
        do k=nvl,1,-1
          write(50) ((real(geopg(j,i,k)),i=1,nlon),j=1,nlat)
        enddo
        do k=nvl,1,-1
          write(50) ((real(psig(j,i,k)),i=1,nlon),j=1,nlat)
        enddo
        do k=nvl,1,-1
          write(50) ((real(qgpv(j,i,k)),i=1,nlon),j=1,nlat)
        enddo
        do k=nvl,1,-1
          write(50) ((real(ug(j,i,k)),i=1,nlon),j=1,nlat)
        enddo
        do k=nvl,1,-1
          write(50) ((real(vg(j,i,k)),i=1,nlon),j=1,nlat)
        enddo
      endif
      
      end
      
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine outputT21(istep)
!-----------------------------------------------------------------------
! *** output T21 truncated data to an ascci outputfile that can be
! *** read by artiforc
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      



      integer istep,k,l,nsh2ntr
      real*8  psiT21(nsh2,nvl),y(nsh2,nvl)

      nsh2ntr=22*23
      if (istep.eq.0) then
        open(51, &
     &file='qgmodelT42.'//expid//'.T21.dat', &
     &form='formatted')
      endif
            
      if (mod(istep,nstepsbetweenoutput).eq.0) then
        call qtopsi
        call truncate(psi,psiT21,21)
        write(51,*)istep/nstepsbetweenoutput
        write(51,'(5e12.5)')((psiT21(k,l),k=1,nsh2ntr),l=1,nvl)
      endif
      
      end
      
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine writestate
!-----------------------------------------------------------------------
! *** output streamfunction state that can be read as initial state
!-----------------------------------------------------------------------

      use ComQG
      implicit none
      



      integer i,j,k,l

      OPEN(12,FILE='qgendT'//trim(ft)//'.dat', &
     &        FORM='FORMATTED')
      do l=1,3
        do k=1,nsh2
           write(12,*) psi(k,l)
        enddo
      enddo
      close(12)
      end
