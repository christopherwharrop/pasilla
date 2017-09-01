module tapenade

    real*8,  parameter :: pi = 4d0 * atan(1d0)     ! value of pi

    integer, parameter :: nm = 21                         ! The truncation is of type T(riangular) nm
    integer, parameter :: nlon = 64                       ! Number of longitude points of the Gaussian grid
    integer, parameter :: nlat = 32                       ! Number of latitude  points of the Gaussian grid
    integer, parameter :: nvl = 3                         ! Number of vorticity levels in the vertical (should be set to 3)
    integer, parameter :: ntl = 2                         ! Number of temperature levels in the vertical (equal to nvl-1)
    integer, parameter :: nsh = ((nm + 1) * (nm + 2)) / 2 ! Half of nsh2
    integer, parameter :: nsh2 = 2 * nsh

    ! Model time step
    integer, parameter :: time_step = 1200
    real*8,  parameter :: nsteps_per_day = 24.0d0 * 3600.0d0 / real(time_step)
    real*8,  parameter :: dtt = (1d0 / nsteps_per_day) * pi * 4d0                      ! dimensionless time step

    ! Model state
    real*8 :: psi(nsh2,nvl)    ! Stream function at the nvl levels
    real*8 :: psit(nsh2,ntl)   ! Thickness at the ntl levels
    real*8 :: qprime(nsh2,nvl) ! Potential vorticity

    ! Model Forcing
    real*8, parameter :: for(nsh2,nvl)    ! Constant potential vorticity forcing at the nvl levels

    ! Spectral Coefficients
    integer, parameter :: nshm(0:nm)     ! Contains numbers 22 down to 1 for index 0 to 21
    integer, parameter :: ll(nsh)       ! Contains total wavenumber n of each spherical harmonic of the corresponding index
    real*8,  parameter :: pp(lat,nsh)     ! Legendre polynomials defined at Gausian latitudes
    real*8,  parameter :: pd(lat,nsh)     ! Mu derivative of Legendre polynomials
    real*8,  parameter :: pw(lat,nsh)     ! Weights for Legendre integrals
    real*8,  parameter :: rm(:)       ! contains zonal wavenumber m of each spherical harmonic of the corresponding index for zonal derivative operator

    ! Laplace/Helmholtz direct and inverse operators
    real*8, parameter :: rinhel(nsh2,0:5) ! Laplace and Helmholtz operator for Q-PSI inversion
    real*8, parameter :: diss(nsh2,2)   ! Dissipation coefficients for each spherical harmonic
                                             !   diss(k,1) : Hyperviscosity at the three levels (tdif sets timescale)
                                             !   diss(k,2) : Ekman friction at lower level (tdis sets timescale)
    real*8, parameter   :: rl1         ! One over Rossby rad. of def. squared of 200-500 thickness
    real*8, parameter   :: rl2         ! One over Rossby rad. of def. squared of 500-800 thickness
    real*8, parameter   :: relt1       ! Nondimensional relaxation coefficient of 200-500 thickness
    real*8, parameter   :: relt2       ! Nondimensional relaxation coefficient of 500-800 thickness

    ! Orography
    real*8, parameter  :: phi(nlat)      ! Gauss points in radians
    real*8, parameter  :: sinfi(nlat)    ! Sine of phi
    real*8, parameter  :: cosfi(nlat)    ! Cosine of phi
    logical, parameter         :: lgdiss      ! If .true. then orography and land-sea mask dependent friction at the lower level plus Ekman friction, else only Ekman friction
    real*8, parameter  :: dorodl(nlat,nlon) ! Derivative of orog wrt lambda
    real*8, parameter  :: dorodm(nlat,nlon) ! Derivative of orag wrt sin(fi)
    real*8, parameter  :: rdiss(nlat,nlon)  ! Landsea-mask/orography dependent friction
    real*8, parameter  :: ddisdx(nlat,nlon) ! Landsea-mask/orography dependent friction
    real*8, parameter  :: ddisdy(nlat,nlon) ! Landsea-mask/orography dependent friction

    real*8  :: agg_copy(nlat,nlon)  ! Copy of input gaussian grid field
    real*8  :: tmp(nlat,nlon)       ! Work space used by the nag version of the fft    

contains

  !-----------------------------------------------------------------------
  ! performs a fourth order runge kutta time step at truncation nm
  ! with time step dt
  ! dqdt calculates the time derivative
  ! input  qprime at current time
  ! output qprime at current time plus dt
  !-----------------------------------------------------------------------
  subroutine adv_nsteps(nsteps)

    integer              :: nsteps

    integer :: step, k, l, nvar
    real*8 :: dt2, dt6
    real*8 :: y(nsh2, nvl), dydt(nsh2, nvl), yt(nsh2, nvl)
    real*8 :: dyt(nsh2, nvl), dym(nsh2, nvl)

    nvar = (nm + 2) * nm
    dt2 = dtt * 0.5d0
    dt6 = dtt / 6d0

    ! Advance the model forward in time n steps
    y = fmtofs(qprime)
    call dqdt(y, dydt)
    do l = 1, nvl
      do k = 1, nvar
        yt(k, l) = y(k, l) + dt2 * dydt(k, l)
      enddo
    enddo
    call dqdt(yt, dyt)
    do l = 1, nvl
      do k = 1, nvar
        yt(k, l) = y(k, l) + dt2 * dyt(k, l)
      enddo
    enddo
    call dqdt(yt, dym)
    do l = 1, nvl
      do k = 1, nvar
        yt(k, l) = y(k, l) + dtt * dym(k, l)
        dym(k, l) = dyt(k, l) + dym(k, l)
      enddo
    enddo
    call dqdt(yt, dyt)
    do l = 1, nvl
      do k = 1, nvar
        y(k, l) = y(k, l) + dt6 * (dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
      enddo
    enddo
    qprime = fstofm(y, nm)

    ! Make stream function consistent with potential vorticity
    call qtopsi(qprime, psi, psit)

  end subroutine adv_nsteps


  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  subroutine dqdt(y, dydt)

    real*8,         intent(   in) :: y(:,:)
    real*8,         intent(  out) :: dydt(:,:)

    real*8 :: local_qprime(nsh2,nvl) ! qprime
    real*8 :: local_psi(nsh2,nvl)    ! psi
    real*8 :: local_psit(nsh2,ntl)   ! psit
    real*8 :: dqprdt(nsh2,nvl) ! time derivative of qprime

    local_qprime = fstofm(y, nm)
    call qtopsi(local_qprime, local_psi, local_psit)
    dqprdt = ddt(local_psi, local_psit, local_qprime, for) ! psi, psit, qprime, for, diss --> dqprdt
    dydt = fmtofs(dqprdt)

    return

  end subroutine dqdt


  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  function ddt(psi, psit, qprime, for) result(dqprdt)

    real*8,         intent(in) :: psi(nsh2,nvl)    ! stream function at the nvl levels
    real*8,         intent(in) :: psit(nsh2,ntl)   ! thickness at the ntl levels
    real*8,         intent(in) :: qprime(nsh2,nvl) ! potential vorticity
    real*8,         intent(in) :: for(nsh2,nvl)    ! constant potential vorticity forcing at the nvl levels
    real*8                     :: dqprdt(nsh2,nvl)

    integer :: k, l, i, j
    real*8 :: dum1, dum2

    ! advection of potential vorticity at upper level
    dqprdt(:, 1) = jacob (psi(:, 1), qprime(:, 1))

    ! advection of potential vorticity at middle level
    dqprdt(:, 2) = jacob (psi(:, 2), qprime(:, 2))

    ! advection of potential vorticity and dissipation at lower level
    dqprdt(:, 3) = jacobd (psi(:, 3), qprime(:, 3))

    ! relaxation of temperature and forcing
    do k = 1, nsh2
      dum1 = relt1 * psit(k, 1)
      dum2 = relt2 * psit(k, 2)
      dqprdt(k, 1) = dqprdt(k, 1) + dum1        + for(k, 1)
      dqprdt(k, 2) = dqprdt(k, 2) - dum1 + dum2 + for(k, 2)
      dqprdt(k, 3) = dqprdt(k, 3)        - dum2 + for(k, 3)
    enddo

    ! explicit horizontal diffusion
    do l = 1, 3
      do k = 1, nsh2
        dqprdt(k, l) = dqprdt(k, l) + diss(k, 1) * qprime(k, l)
      enddo
    enddo

    return

  end function ddt


  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacob (psiloc, pvor) result(sjacob)

    implicit none

    real*8, intent( in) :: psiloc(nsh2)
    real*8, intent( in) :: pvor(nsh2)
    real*8              :: sjacob(nsh2)

    integer      :: i, j, k
    real*8 :: vv(nsh2)
    real*8 :: dpsidl(nlat, nlon),  dpsidm(nlat, nlon),  dvordl(nlat, nlon)
    real*8 :: dvordm(nlat, nlon),  gjacob(nlat, nlon),  dpsidls(nsh2)

    ! space derivatives of potential vorticity
    vv = ddl (pvor)
    dvordl = sptogg_pp (vv)
    dvordm = sptogg_pd (pvor)

    ! space derivatives of streamfunction
    dpsidls = ddl (psiloc)
    dpsidl = sptogg_pp (dpsidls)
    dpsidm = sptogg_pd (psiloc)

    ! jacobian term
    do j = 1, nlon
      do i = 1, nlat
        gjacob(i, j) = dpsidm(i, j) * dvordl(i, j) - dpsidl(i, j) * dvordm(i, j)
      enddo
    enddo

    sjacob = ggtosp (gjacob)

    ! planetary vorticity advection
    do k = 1, nsh2
      sjacob(k) = sjacob(k) - dpsidls(k)
    enddo

    return

  end function jacob


  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  function jacobd (psiloc, pvor) result(sjacob)

    real*8,         intent(in) :: psiloc(nsh2)
    real*8,         intent(in) :: pvor(nsh2)
    real*8                     :: sjacob(nsh2)

    integer      :: i, j, k
    real*8 :: dpsidl(nlat, nlon),  dpsidm(nlat, nlon),  dvordl(nlat, nlon)
    real*8 :: dvordm(nlat, nlon),  gjacob(nlat, nlon),  vv(nsh2)
    real*8 :: azeta(nlat, nlon), dpsidls(nsh2)

    ! space derivatives of potential vorticity 
    vv = ddl (pvor)
    dvordl = sptogg_pp (vv)
    dvordm = sptogg_pd (pvor)

    ! space derivatives of streamfunction
    dpsidls = ddl (psiloc)
    dpsidl = sptogg_pp (dpsidls)
    dpsidm = sptogg_pd (psiloc)

    ! jacobian term + orographic forcing
    do j = 1, nlon
      do i = 1, nlat
        gjacob(i, j) = dpsidm(i, j) * (dvordl(i, j) + sinfi(i) * dorodl(i, j)) -  &
   &                   dpsidl(i, j) * (dvordm(i, j) + sinfi(i) * dorodm(i, j))
      enddo
    enddo

    ! dissipation 
    if (lgdiss) then

      !   spatially varying dissipation 
      do k = 1, nsh2
        vv(k) = diss(k, 2) * psiloc(k)
      enddo

      azeta = sptogg_pp (vv)

      do j = 1, nlon
        do i = 1, nlat
          gjacob(i, j) = gjacob(i, j) - dpsidm(i, j)     * ddisdy(i, j) &
   &                                  - dpsidl(i, j)     * ddisdx(i, j) &
   &                                  + rdiss(i, j) * azeta(i, j)
        enddo
      enddo

      sjacob = ggtosp (gjacob)

    else

      !   uniform dissipation
      sjacob = ggtosp (gjacob)

      do k = 1, nsh2
        sjacob(k) = sjacob(k) + diss(k, 2) * psiloc(k)
      enddo

    endif

    ! planetary vorticity advection
    do k = 1, nsh2
      sjacob(k) = sjacob(k) - dpsidls(k)
    enddo

    return

  end function jacobd


  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  subroutine qtopsi(qprime, psi, psit)

    real*8,         intent( in) :: qprime(:,:) ! potential vorticity
    real*8,         intent(out) :: psi(:,:)    ! stream function at the nvl levels
    real*8,         intent(out) :: psit(:,:)   ! thickness at the ntl levels

    integer :: k
    real*8 :: r3
    real*8 :: ws(nsh2)       ! only used as portable workspace

    do k = 1, size(psi,1)
      ws(k) = qprime(k, 1) + qprime(k, 3)
      psi(k, 1) = rinhel(k, 1) * (ws(k) + qprime(k, 2))
      psi(k, 2) = ws(k) - 2.d0 * qprime(k, 2)
      psi(k, 3) = qprime(k, 1) - qprime(k, 3)
    enddo

    do k = 1, size(psit,1)
      psit(k, 1) = rinhel(k, 2) * psi(k, 2) + rinhel(k, 3) * psi(k, 3)
      psit(k, 2) = rinhel(k, 4) * psi(k, 2) + rinhel(k, 5) * psi(k, 3)
    enddo

    r3 = 1. / 3.
    do k = 1, size(psi,1)
      psi(k, 2) = r3 * (psi(k, 1) - psit(k, 1) + psit(k, 2))
      psi(k, 1) = psi(k, 2) + psit(k, 1)
      psi(k, 3) = psi(k, 2) - psit(k, 2)
    enddo

    return

  end subroutine qtopsi


  !-----------------------------------------------------------------------
  ! computation of potential vorticity from stream function
  ! input psi streamfunction
  ! output qprime,  the potential vorticity and psit,  the layer thick.
  !-----------------------------------------------------------------------
  subroutine psitoq(psi, psit, qprime)
      
    real*8,         intent( in) :: psi(:,:)    ! stream function at the nvl levels
    real*8,         intent(out) :: psit(:,:)   ! thickness at the ntl levels
    real*8,         intent(out) :: qprime(:,:) ! potential vorticity

    integer :: k

    do k = 1, size(psit,1)
      psit(k, 1) = psi(k, 1) - psi(k, 2)
      psit(k, 2) = psi(k, 2) - psi(k, 3)
      qprime(k, 1) = rinhel(k, 0) * psi(k, 1) - rl1 * psit(k, 1)
      qprime(k, 2) = rinhel(k, 0) * psi(k, 2) + rl1 * psit(k, 1) - rl2 * psit(k, 2)
      qprime(k, 3) = rinhel(k, 0) * psi(k, 3) + rl2 * psit(k, 2)
    enddo

    return

  end subroutine psitoq


  !-----------------------------------------------------------------------
  ! transforms francos format to the french format for global fields
  ! input  y spectral coefficients in francos format
  ! output z spectral coefficients in french format
  ! fm format:
  ! k       m  n
  ! 1       0  0
  ! 2       0  1
  ! 3       0  2
  ! :       :  :
  ! nm+1    0  nm
  ! nm+2    1  1 --> real part
  ! nm+3    1  2 --> real part
  ! :       :  :
  ! nm+nm+1 1  nm --> real part
  ! :       :  :
  ! :       nm nm --> real part
  !  repeat for imaginary part
  !  disadvantage: 0 0 mode and imaginary parts of m = 0 modes are obsolete
  ! fs format stores all m for every n first and has no obsolete indices
  ! 
  ! k       m  n
  ! 1       0  1
  ! 2       1  1 --> real part
  ! 3       1  1 --> imaginary part: k = 1-3 is T1 truncation
  ! 4       0  2
  ! 5       1  2 --> real part
  ! 6       1  2 --> imaginary part
  ! 7       2  2 --> real part
  ! 8       2  2 --> imaginary part: k = 1-8 is T2 truncation
  ! etcetera
  !-----------------------------------------------------------------------
  function fmtofs (y) result(z)

    real*8, intent( in)        :: y(:,:)

    real*8, dimension(size(y,1),size(y,2)) :: z

    integer ::  m, n, k, indx, l

    do l = 1, size(y,2)
      k = 1
      do m = 0, nm
        do n = max(m, 1), nm
          k = k + 1
          if (m .eq. 0) then
            indx = n**2
          else
            indx = n**2 + 2 * m - 1
          end if
          z(indx, l) = y(k, l)
          if (m .ne. 0) z(indx + 1, l) = y(k + nsh, l)
        enddo
      enddo
    enddo

    return

  end function fmtofs


  !-----------------------------------------------------------------------
  ! transforms the french format to francos format for global fields
  ! input  y spectral coef. in french format,  ntr is truncation limit
  ! output z spectral coefficients in francos format
  ! fm format:
  ! k       m  n
  ! 1       0  0
  ! 2       0  1
  ! 3       0  2
  ! :       :  :
  ! nm+1    0  nm
  ! nm+2    1  1 --> real part
  ! nm+3    1  2 --> real part
  ! :       :  :
  ! nm+nm+1 1  nm --> real part
  ! :       :  :
  ! :       nm nm --> real part
  !  repeat for imaginary part
  !  disadvantage: 0 0 mode and imaginary parts of m = 0 modes are obsolete
  ! fs format stores all m for every n first and has no obsolete indices
  ! 
  ! k       m  n
  ! 1       0  1
  ! 2       1  1 --> real part
  ! 3       1  1 --> imaginary part: k = 1-3 is T1 truncation
  ! 4       0  2
  ! 5       1  2 --> real part
  ! 6       1  2 --> imaginary part
  ! 7       2  2 --> real part
  ! 8       2  2 --> imaginary part: k = 1-8 is T2 truncation
  ! etcetera
  !-----------------------------------------------------------------------
  function fstofm (y, ntr) result(z)

    real*8,         intent(in) :: y(:,:)
    integer,              intent(in) :: ntr

    real*8, dimension(size(y,1),size(y,2)) :: z

    integer :: m, n, k, indx, i, l

    do l = 1, size(y,2)
      do i = 1, size(y,1)
        z(i, l) = 0d0
      enddo
      k = 1
      do m = 0, nm
        do n = max(m, 1), nm
          k = k + 1
          if ((m .le. ntr).and.(n .le. ntr)) then
            if (m .eq. 0) then
              indx = n**2
            else
              indx = n**2 + 2 * m - 1
            end if
            z(k, l) = y(indx, l)
            if (m .ne. 0) z(k + nsh, l) = y(indx + 1, l)
          endif
        enddo
      enddo
    enddo

    return

  end function fstofm


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid using
  ! legendre polynomials
  !
  ! input  spectral field as
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg_pp (as) result(agg)

    ! Input
    real*8,        intent(   in) :: as(nsh, 2)

    ! Return value
    real*8 :: agg(nlat, nlon)

    agg = sptogg(as, pp)

  end function sptogg_pp


  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid using
  ! derivatives with respect to sin(fi) .
  !
  ! input  spectral field as
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg_pd (as) result(agg)

    ! Input
    real*8,        intent(   in) :: as(nsh, 2)

    ! Return value
    real*8 :: agg(nlat, nlon)

    agg = sptogg(as, pd)

  end function sptogg_pd

 
  !-----------------------------------------------------------------------
  ! conversion from spectral coefficients to gaussian grid
  ! input  spectral field as,  legendre polynomials pploc (pp or pd) 
  !        where pp are legendre polynomials and pd derivatives with
  !        respect to sin(fi)
  ! output gaussian grid agg
  !-----------------------------------------------------------------------
  function sptogg (as, pploc) result(agg)

    ! Input 
    real*8,        intent(   in) :: as(:,:)
    real*8,        intent(   in) :: pploc(:,:)

    ! Return value
    real*8 :: agg(nlat, nlon)

    ! Local data
!    real*8, allocatable :: tmp(:,:)   ! Work space used by the nag version of the fft    
    integer :: i, j, k, k1, k2, m, mi, mr, nlon1
    integer :: ifail

    ! inverse legendre transform
    do j = 1, nlon
      do i = 1, nlat
        agg(i, j) = 0.0d0
      enddo
    enddo

    nlon1 = nlon + 1
    k2 = nshm(0)

    do k = 1, k2
      do i = 1, nlat
        agg(i, 1) = agg(i, 1) + as(k, 1) * pploc(i, k)
      enddo
    enddo

    do m = 1, nm
      mr = m + 1
      mi = nlon1 - m
      k1 = k2 + 1
      k2 = k2 + nshm(m)
      do k = k1, k2
        do i = 1, nlat
          agg(i, mr) = agg(i, mr) + as(k, 1) * pploc(i, k)
        enddo
        do i = 1, nlat
          agg(i, mi) = agg(i, mi) - as(k, 2) * pploc(i, k)
        enddo
      enddo
    enddo

    ! inverse fourier transform
    ifail = 0
 !   allocate(tmp(nlat,nlon))
    call c06fqf (nlat, nlon, agg, 'r', trigi, tmp, ifail)

  end function sptogg
 

  !-----------------------------------------------------------------------
  ! conversion from gaussian grid (agg) to spectral coefficients (as)
  ! input gaussian grid field agg
  ! output as contains spectral coefficients
  !-----------------------------------------------------------------------
  function ggtosp (agg) result (as)

    ! Input
    real*8,        intent(   in) :: agg(:,:)

    ! Return value
    real*8 :: as(nsh, 2)

    ! Local data
    integer :: i, k, k1, k2, m, mi, mr, nlon1
    integer :: ifail

    ! Make a local copy of agg so it is not destroyed by c06fpf
!    allocate(agg_copy(nlat,nlon))
    agg_copy(:,:) = agg(:,:)

    ! fourier transform
    ifail = 0
!    allocate(tmp(nlat,nlon))
    call c06fpf (nlat, nlon, agg_copy, 'r', trigd, tmp, ifail)

    ! legendre transform
    do i = 1, 2
      do k = 1, nsh
        as(k, i) = 0.0d0
      enddo
    enddo

    nlon1 = nlon + 1

    k2 = nshm(0)

    do k = 1, k2
      do i = 1, nlat
        as(k, 1) = as(k, 1) + agg_copy(i, 1) * pw(i, k)
      enddo
    enddo

    do m = 1, nm
      mr = m + 1
      mi = nlon1 - m
      k1 = k2 + 1
      k2 = k2 + nshm(m)
      do k = k1, k2
        do i = 1, nlat
          as(k, 1) = as(k, 1) + agg_copy(i, mr) * pw(i, k)
          as(k, 2) = as(k, 2) + agg_copy(i, mi) * pw(i, k)
        enddo
      enddo
    enddo

  end function ggtosp


  !-----------------------------------------------------------------------
  ! zonal derivative in spectral space
  ! input spectral field as
  ! output spectral field dadl which is as differentiated wrt lambda
  !-----------------------------------------------------------------------
  function ddl (as) result(dadl)

    implicit none

    real*8,        intent(in) :: as(nsh, 2)
    real*8                    :: dadl(nsh, 2)

    integer :: k

    do k = 1, nsh
      dadl(k, 1) = -rm(k) * as(k, 2)
      dadl(k, 2) =  rm(k) * as(k, 1)
    enddo

    return

  end function ddl


end module tapenade
