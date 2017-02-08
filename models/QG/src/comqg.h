!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comqg.h
! *** Contents: Common declarations for three level qg model
!-----------------------------------------------------------------------
! *** COMMON  /ch/ ft,expid
!     ft :     character string with truncation
!     expid:   character string with experiment identifier used to create
!              the outputdata data directory
!     rootdir: path of qgmodel directory
!
!
! *** COMMON  /rvari/ pi,dp,om,rgas,grav,radius
!     pi :     value of pi
!     fzero:   value of f at 45 degrees
!     dp :     layer thicknesses [Pa]
!     om :     angular velocity of earth [rad/s]
!     rgas :   gas constant
!     grav :   gravity acceleration [m/s^2]
!     radius:  radius of earth [m]
!
!
! *** COMMON  /ggrid/ phi,cosfi,sinfi
!     phi :    Gauss points in radians
!     cosfi:   cosine of phi
!     sinfi:   sine of phi
!
! *** COMMON  /ctstep/ dt,dtime,dtt
!     dt :     timestep in fraction of one day
!     dtime :  timestep in seconds
!     dtt :    timestep in non-dimensional units
!
! *** COMMON  /intpar/ nshm, ll
!     nshm:   contains numbers 22 down to 1 for index 0 to 21
!     ll:     contains total wavenumber n of each spherical harmonic of
!             the corresponding index
!
! *** COMMON  /linopr/ rm, rinhel, diss 
!     rm:     contains zonal wavenumber m of each spherical harmonic of
!             the corresponding index for zonal derivative operator
!     rinhel: Laplace and Helmholtz operator for Q-PSI inversion
!     diss:   dissipation coefficients for each spherical harmonic
!             diss(k,1) : hyperviscosity at the three levels 
!                         (tdif sets timescale)
!             diss(k,2) : Ekman friction at lower level 
!                         (tdis sets timescale)
!
! *** COMMON  /logpar/ lgdiss,inf,readfor
!     lgdiss: if .true. then orography and land-sea mask dependent 
!             friction at the lower level plus Ekman friction, 
!             else only Ekman friction
!     inf:    if .true. then artificial PV forcing read from file
!     obsf:   if .true. PV forcing is calculated from observations in routine artiforc
!     readstart: if .true. initial state is read from inputfile
!
! *** COMMON  /metras/
!     pp:     Legendre polynomials defined at Gausian latitudes
!     pd:     mu derivative of Legendre polynomials
!     pw:     weights for Legendre integrals
!
! *** COMMON  /phypar/ rdiss, ddisdx, ddisdy
!     rdiss:  landsea-mask/orography dependent friction
!     ddisdx: landsea-mask/orography dependent friction
!     ddisdy: landsea-mask/orography dependent friction
!
! *** COMMON  /modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
!                      tdis,trel,tdif,addisl,addish,h0,idif
!     rrdef1: Rossby radius of deformation of 200-500 thickness
!     rrdef2: Rossby radius of deformation of 500-800 thickness
!     rl1:    one over Rossby rad. of def. squared of 200-500 thickness
!     rl2:    one over Rossby rad. of def. squared of 500-800 thickness
!     relt1:  nondimensional relaxation coefficient of 200-500 thickness
!     relt2:  nondimensional relaxation coefficient of 500-800 thickness
!     tdis:   Ekman dissipation timescale in days at lower level
!     trel:   relaxation time scale in days of the temperature
!     tdif:   dissipation timescale of scale-selective diffusion in
!             days for wavenumber nm
!     addisl: parameter used in the computation of the dissipation
!             timescale at the lower level over land
!     addish: parameter used in the computation of the dissipation
!             timescale at the lower level as a function of topography
!     h0:     scale factor for the topographically induced upward motion
!             at the lower level
!     idif:   determines scale-selectivity of hyperviscosity; power of
!             laplace operator
! 
! *** COMMON  /sfield/ psi, psit, qprime, dqprdt, for, ws
!     psi:    stream function at the nvl levels
!     psit:   thickness at the ntl levels
!     qprime: potential vorticity
!     dqprdt: time derivative of qprime
!     for:    constant potential vorticity forcing at the nvl levels
!     ws:     only used as portable workspace
!
! *** COMMON  /corog/ orog, dorodl, dorodm
!     orog:   orography in m. divided by h0
!     dorodl: derivative of orog wrt lambda
!     dorodm: derivative of orag wrt sin(fi)
!
! *** COMMON  /zotras/ trigd, trigi, wgg
!             arrays used by the nag version of the fft
!
! *** COMMON /cgpsi/  psig,qgpv,ug,vg
!     psig: grid values of dimensional streamfunction at the three levels
!     qgpv: grid values of dimensional pv at the three levels
!     ug:   grid values of zonal velocity at the three levels in m/s
!     vg:   grid values of meridional velocity at the three levels in m/s
!     geopg: geopotential on the grid
!-----------------------------------------------------------------------

      character*4 expid
      character*32 obsfile
      character*256 rootdir
      integer nstepsperday,nstepsbetweenoutput,ndayskip,nday,rootdirl
      real*8  pi,fzero,dp,om,rgas,grav,radius
      real*8  phi(nlat),cosfi(nlat),sinfi(nlat)
      real*8  dt,dtime,dtt

      common /ch/expid,obsfile,rootdir
      common /rvari/ pi,fzero,dp,om,rgas,grav,radius
      common /ggrid/ phi,cosfi,sinfi
      common /ctstep/ dt,dtime,dtt
     
      integer nshm(0:nm),ll(nsh)
      real*8  rm(nsh),rinhel(nsh2,0:5),diss(nsh2,2)
      logical lgdiss,inf,obsf,readstart
      real*8  pp(nlat,nsh),pd(nlat,nsh),pw(nlat,nsh)
      real*8  rdiss(nlat,nlon),ddisdx(nlat,nlon),ddisdy(nlat,nlon)
      real*8  rrdef1,rrdef2,rl1,rl2,relt1,relt2,tdis,trel,tdif
      real*8  addisl,addish,h0
      integer idif

      common /intpar/ nshm,ll,nstepsperday,nstepsbetweenoutput, &
     &                ndayskip,nday,rootdirl
      common /linopr/ rm, rinhel, diss 
      common /logpar/ lgdiss,inf,obsf,readstart
      common /metras/ pp, pd, pw
      common /phypar/ rdiss, ddisdx, ddisdy
      common /modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 , &
     &                tdis,trel,tdif,addisl,addish,h0,idif
			
			
      real*8  psi(nsh2,nvl),psit(nsh2,ntl)
      real*8  qprime(nsh2,nvl),dqprdt(nsh2,nvl),for(nsh2,nvl),ws(nsh2)
      real*8  orog(nsh2),dorodl(nlat,nlon),dorodm(nlat,nlon)
      real*8  trigd(nlon,2),trigi(nlon,2),wgg(nlat,nlon)
      real*8  psig(nlat,nlon,nvl),qgpv(nlat,nlon,nvl)
      real*8  ug(nlat,nlon,nvl),vg(nlat,nlon,nvl),geopg(nlat,nlon,nvl)


      common /sfield/ psi, psit, qprime, dqprdt, for, ws
      common /corog/  orog, dorodl , dorodm
      common /zotras/ trigd, trigi, wgg
      common /cgpsi/  psig,qgpv,ug,vg,geopg
