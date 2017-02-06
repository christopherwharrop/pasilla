c23456789012345678901234567890123456789012345678901234567890123456789012
c-----------------------------------------------------------------------
c *** File:     comqg.h                                                    
c *** Contents: Common declarations for three level qg model
c-----------------------------------------------------------------------
c *** COMMON  /ch/ ft,expid
c     ft :     character string with truncation
c     expid:   character string with experiment identifier used to create
c              the outputdata data directory
c     rootdir: path of qgmodel directory
c
c
c *** COMMON  /rvari/ pi,dp,om,rgas,grav,radius
c     pi :     value of pi
c     fzero:   value of f at 45 degrees
c     dp :     layer thicknesses [Pa]
c     om :     angular velocity of earth [rad/s]
c     rgas :   gas constant
c     grav :   gravity acceleration [m/s^2]
c     radius:  radius of earth [m]
c
c
c *** COMMON  /ggrid/ phi,cosfi,sinfi
c     phi :    Gauss points in radians
c     cosfi:   cosine of phi
c     sinfi:   sine of phi
c
c *** COMMON  /ctstep/ dt,dtime,dtt
c     dt :     timestep in fraction of one day
c     dtime :  timestep in seconds
c     dtt :    timestep in non-dimensional units
c
c *** COMMON  /intpar/ nshm, ll
c     nshm:   contains numbers 22 down to 1 for index 0 to 21
c     ll:     contains total wavenumber n of each spherical harmonic of
c             the corresponding index
c
c *** COMMON  /linopr/ rm, rinhel, diss 
c     rm:     contains zonal wavenumber m of each spherical harmonic of
c             the corresponding index for zonal derivative operator
c     rinhel: Laplace and Helmholtz operator for Q-PSI inversion
c     diss:   dissipation coefficients for each spherical harmonic
c             diss(k,1) : hyperviscosity at the three levels 
c                         (tdif sets timescale)
c             diss(k,2) : Ekman friction at lower level 
c                         (tdis sets timescale)
c
c *** COMMON  /logpar/ lgdiss,inf,readfor
c     lgdiss: if .true. then orography and land-sea mask dependent 
c             friction at the lower level plus Ekman friction, 
c             else only Ekman friction
c     inf:    if .true. then artificial PV forcing read from file
c     obsf:   if .true. PV forcing is calculated from observations in routine artiforc
c     readstart: if .true. initial state is read from inputfile
c
c *** COMMON  /metras/
c     pp:     Legendre polynomials defined at Gausian latitudes
c     pd:     mu derivative of Legendre polynomials
c     pw:     weights for Legendre integrals
c
c *** COMMON  /phypar/ rdiss, ddisdx, ddisdy
c     rdiss:  landsea-mask/orography dependent friction
c     ddisdx: landsea-mask/orography dependent friction
c     ddisdy: landsea-mask/orography dependent friction
c
c *** COMMON  /modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
c                      tdis,trel,tdif,addisl,addish,h0,idif
c     rrdef1: Rossby radius of deformation of 200-500 thickness
c     rrdef2: Rossby radius of deformation of 500-800 thickness
c     rl1:    one over Rossby rad. of def. squared of 200-500 thickness
c     rl2:    one over Rossby rad. of def. squared of 500-800 thickness
c     relt1:  nondimensional relaxation coefficient of 200-500 thickness
c     relt2:  nondimensional relaxation coefficient of 500-800 thickness
c     tdis:   Ekman dissipation timescale in days at lower level
c     trel:   relaxation time scale in days of the temperature
c     tdif:   dissipation timescale of scale-selective diffusion in
c             days for wavenumber nm
c     addisl: parameter used in the computation of the dissipation
c             timescale at the lower level over land
c     addish: parameter used in the computation of the dissipation
c             timescale at the lower level as a function of topography
c     h0:     scale factor for the topographically induced upward motion
c             at the lower level
c     idif:   determines scale-selectivity of hyperviscosity; power of
c             laplace operator
c 
c *** COMMON  /sfield/ psi, psit, qprime, dqprdt, for, ws
c     psi:    stream function at the nvl levels
c     psit:   thickness at the ntl levels
c     qprime: potential vorticity
c     dqprdt: time derivative of qprime
c     for:    constant potential vorticity forcing at the nvl levels
c     ws:     only used as portable workspace
c
c *** COMMON  /corog/ orog, dorodl, dorodm
c     orog:   orography in m. divided by h0
c     dorodl: derivative of orog wrt lambda
c     dorodm: derivative of orag wrt sin(fi)
c
c *** COMMON  /zotras/ trigd, trigi, wgg
c             arrays used by the nag version of the fft
c
c *** COMMON /cgpsi/  psig,qgpv,ug,vg
c     psig: grid values of dimensional streamfunction at the three levels
c     qgpv: grid values of dimensional pv at the three levels
c     ug:   grid values of zonal velocity at the three levels in m/s
c     vg:   grid values of meridional velocity at the three levels in m/s
c     geopg: geopotential on the grid
c-----------------------------------------------------------------------

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

      common /intpar/ nshm,ll,nstepsperday,nstepsbetweenoutput,
     *                ndayskip,nday,rootdirl
      common /linopr/ rm, rinhel, diss 
      common /logpar/ lgdiss,inf,obsf,readstart
      common /metras/ pp, pd, pw
      common /phypar/ rdiss, ddisdx, ddisdy
      common /modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
     *                tdis,trel,tdif,addisl,addish,h0,idif
			
			
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
