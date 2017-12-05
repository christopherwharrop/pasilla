submodule(QG_Model) QG_Model_TL

!  use kind
!  use QG_Config
!  use QG_GGSP
!  use QG_Util

  implicit none

contains


  !-------------------------------------------------------------------------------
  ! constructor_qg_tl
  !-------------------------------------------------------------------------------
  module function constructor_qg_tl(config, state, state_vector, trajectory, trajectory_vector, for, step) result (qg_tl)

    type(qg_config_type),   intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:,:)
    real(r8kind), optional, intent(in) :: state_vector(:)
    real(r8kind), optional, intent(in) :: trajectory(:,:)
    real(r8kind), optional, intent(in) :: trajectory_vector(:)
    real(r8kind), optional, intent(in) :: for(:,:)
    integer,      optional, intent(in) :: step
    type(qg_tl_type)                   :: qg_tl

    real(r8kind), allocatable :: psig3d(:,:,:), psisp(:,:)
    real(r8kind), allocatable :: trajg3d(:,:,:)

    real(r8kind) :: nsteps_per_day

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    ! Set the model config
    qg_tl%config = config

    ! Get model grid dimensions
    qg_tl%ft = config%get_ft()
    qg_tl%nm = config%get_nm()
    qg_tl%nlon = config%get_nlon()
    qg_tl%nlat = config%get_nlat()
    qg_tl%nvl = config%get_nvl()
    qg_tl%ntl = config%get_ntl()
    qg_tl%nsh = config%get_nsh()
    qg_tl%nsh2 = config%get_nsh2()

    ! Get model grid converter
    qg_tl%ggsp = config%get_ggsp()

    ! Get model forcing from configuration
    allocate(qg_tl%for, source = config%get_for())

    ! Get Nondimensional relaxation coefficients
    qg_tl%relt1 = config%get_relt1()
    qg_tl%relt2 = config%get_relt2()

    ! Get dissipation coefficients for each spherical harmonic
    allocate(qg_tl%diss, source = config%get_diss())

    ! Get gauss points in radians and sine and cosine of it
    allocate(qg_tl%phi, source = config%get_phi())
    allocate(qg_tl%sinfi, source = config%get_sinfi())
    allocate(qg_tl%cosfi, source = config%get_cosfi())

    ! Get derivatives of orog
    allocate(qg_tl%dorodl, source = config%get_dorodl())
    allocate(qg_tl%dorodm, source = config%get_dorodm())

    ! Get orography and land-sea mask dependent friction option
    qg_tl%lgdiss = config%get_lgdiss()

    ! Get landsea-mask/orography dependent friction
    allocate(qg_tl%rdiss, source = config%get_rdiss())
    allocate(qg_tl%ddisdx, source = config%get_ddisdx())
    allocate(qg_tl%ddisdy, source = config%get_ddisdy())

    ! Get Laplace and Helmholtz operator for Q-PSI inversion
    allocate(qg_tl%rinhel(qg_tl%nsh2,0:5))
    qg_tl%rinhel = config%get_rinhel()

    ! Get one over Rossby rad. of def. squared of 200-500, and 500-800 thicknesses
    qg_tl%rl1 = config%get_rl1()
    qg_tl%rl2 = config%get_rl2()

    ! Initialize model step and clock
    if (present(step)) then
      qg_tl%step = step
    else
      qg_tl%step = 0
    end if
    qg_tl%clock = qg_tl%step * config%get_time_step()

    ! Initialize time step of the model:
    nsteps_per_day = 24.0d0 * 3600.0d0 / real(config%get_time_step())
    qg_tl%dtt = (1d0 / nsteps_per_day) * pi * 4d0

    ! Initialize streamfunction
    if (present(state)) then
      call qg_tl%init_state(psi=state)
    else if (present(state_vector)) then
      allocate(psig3d(qg_tl%nlon, qg_tl%nlat, qg_tl%nvl))
      allocate(psisp(qg_tl%nsh2, qg_tl%nvl))
      psig3d = reshape(state_vector,(/qg_tl%nlon, qg_tl%nlat, qg_tl%nvl/))
      psig3d(:,:,:) = psig3d(:,:,:) / facsf
      psisp(:,1) = reshape(qg_tl%ggsp%ggtosp(transpose(psig3d(:,:,1))), (/qg_tl%nsh2/))
      psisp(:,2) = reshape(qg_tl%ggsp%ggtosp(transpose(psig3d(:,:,2))), (/qg_tl%nsh2/))
      psisp(:,3) = reshape(qg_tl%ggsp%ggtosp(transpose(psig3d(:,:,3))), (/qg_tl%nsh2/))
      call qg_tl%init_state(psi=psisp)
    else
      call qg_tl%init_state()
    end if

    ! Initialize trajectory
    if (present(trajectory)) then
      allocate(qg_tl%trajectory, source = trajectory)
    else if (present(trajectory_vector)) then
      allocate(trajg3d(qg_tl%nlon, qg_tl%nlat, qg_tl%nvl))
      allocate(qg_tl%trajectory(qg_tl%nsh2,qg_tl%nvl))
      trajg3d = reshape(trajectory_vector,(/qg_tl%nlon, qg_tl%nlat, qg_tl%nvl/))
      trajg3d(:,:,:) = trajg3d(:,:,:) / facsf
      qg_tl%trajectory(:,1) = reshape(qg_tl%ggsp%ggtosp(transpose(trajg3d(:,:,1))), (/qg_tl%nsh2/))
      qg_tl%trajectory(:,2) = reshape(qg_tl%ggsp%ggtosp(transpose(trajg3d(:,:,2))), (/qg_tl%nsh2/))
      qg_tl%trajectory(:,3) = reshape(qg_tl%ggsp%ggtosp(transpose(trajg3d(:,:,3))), (/qg_tl%nsh2/))
    else
      allocate(qg_tl%trajectory, source = qg_tl%psi)
    end if
    allocate(qg_tl%psitd(qg_tl%nsh2,qg_tl%ntl))
    allocate(qg_tl%qprimed(qg_tl%nsh2,qg_tl%nvl))
    call qg_tl%psitoq(qg_tl%trajectory, qg_tl%psitd, qg_tl%qprimed)

  end function constructor_qg_tl


  !------------------------------------------------------------------
  ! destructor_qg_tl
  !
  ! Deallocates pointers used by a qg_tl_type object (none currently)
  !------------------------------------------------------------------
  module elemental subroutine destructor_qg_tl(this)

    type(qg_tl_type), intent(inout) :: this

    ! No pointers in qg_tl_type object so we do nothing

  end subroutine destructor_qg_tl


  !-----------------------------------------------------------------------
  ! performs a fourth order runge kutta time step at truncation nm
  ! with time step dt
  ! dqdt calculates the time derivative
  ! input  qprime at current time
  ! output qprime at current time plus dt
  !-----------------------------------------------------------------------
  module subroutine adv_nsteps_tl(this, nsteps)

    class(qg_tl_type) :: this
    integer           :: nsteps

    integer :: step, k, l, nvar
    real(r8kind) :: dt2, dt6
    real(r8kind) :: y(this%nsh2, this%nvl), dydt(this%nsh2, this%nvl), yt(this%nsh2, this%nvl)
    real(r8kind) :: dyt(this%nsh2, this%nvl), dym(this%nsh2, this%nvl)
    real(r8kind) :: yd(this%nsh2, this%nvl), dydtd(this%nsh2, this%nvl), ytd(this%nsh2, this%nvl)
    real(r8kind) :: dytd(this%nsh2, this%nvl), dymd(this%nsh2, this%nvl)

    if (nsteps > 0) then

      nvar = (this%nm + 2) * this%nm
      dt2 = this%dtt * 0.5d0
      dt6 = this%dtt / 6d0

      ! Advance the model forward in time n steps
      do step = 1, nsteps

        yd = this%fmtofs(this%qprimed)
        y = this%fmtofs(this%qprime)
        call this%dqdt_d(y, yd, dydt, dydtd)
        ytd = 0.0_8
        do l = 1, this%nvl
          do k = 1, nvar
            ytd(k, l) = yd(k, l) + dt2 * dydtd(k, l)
            yt(k, l) = y(k, l) + dt2 * dydt(k, l)
          end do
        end do
        call this%dqdt_d(yt, ytd, dyt, dytd)
        do l = 1, this%nvl
          do k = 1, nvar
            ytd(k, l) = yd(k, l) + dt2 * dytd(k, l)
            yt(k, l) = y(k, l) + dt2 * dyt(k, l)
          end do
        end do
        call this%dqdt_d(yt, ytd, dym, dymd)
        do l = 1, this%nvl
          do k = 1, nvar
            ytd(k, l) = yd(k, l) + this%dtt * dymd(k, l)
             yt(k, l) = y(k, l) + this%dtt * dym(k, l)
            dymd(k, l) = dytd(k, l) + dymd(k, l)
             dym(k, l) = dyt(k, l) + dym(k, l)
          end do
        end do
        call this%dqdt_d(yt, ytd, dyt, dytd)
        do l = 1, this%nvl
          do k = 1, nvar
            yd(k, l) = yd(k, l) + dt6 * (dydtd(k, l) + dytd(k, l) + 2. * dymd(k, l))
            y(k, l)  = y(k, l) + dt6 * ( dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
          end do
        end do

        ! Update model state with original trajectory
        this%qprime = this%qprime + this%qprimed

        ! Update trajectory
        this%qprimed = this%fstofm(yd, this%nm)

        ! Inrement the step count
        this%step = this%step + 1

      end do

      ! Make stream function consistent with potential vorticity
      call this%qtopsi(this%qprimed, this%trajectory, this%psitd)
      call this%qtopsi(this%qprime, this%psi, this%psit)

    end if

  end subroutine adv_nsteps_tl


  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  module subroutine dqdt_d(this, y, yd, dydt, dydtd)

    class(qg_tl_type) :: this
    real(r8kind), intent( in) :: y(:,:)
    real(r8kind), intent( in) :: yd(:,:)
    real(r8kind), intent(out) :: dydt(:,:)
    real(r8kind), intent(out) :: dydtd(:,:)

    real(r8kind) :: qprime(this%nsh2,this%nvl) ! qprime
    real(r8kind) :: qprimed(this%nsh2,this%nvl) ! qprime
    real(r8kind) :: psi(this%nsh2,this%nvl)    ! psi
    real(r8kind) :: psid(this%nsh2,this%nvl)    ! psi
    real(r8kind) :: psit(this%nsh2,this%ntl)   ! psit
    real(r8kind) :: psitd(this%nsh2,this%ntl)   ! psit
    real(r8kind) :: dqprdt(this%nsh2, this%nvl)  ! time derivative of qprime
    real(r8kind) :: dqprdtd(this%nsh2, this%nvl) ! time derivative of qprime

    qprimed = this%fstofm(yd, this%nm)
    qprime = this%fstofm(y, this%nm)

    call this%qtopsi(qprimed, psid, psitd) ! qprimed --> psid and psitd
    call this%qtopsi(qprime, psi, psit)    ! qprime --> psi and psit

    ! psi, psit, qprime, for, diss --> dqprdt
    call this%ddt_d(psi, psid, psit, psitd, qprime, qprimed, this%for, dqprdt, dqprdtd)

    dydtd = this%fmtofs(dqprdtd)
    dydt = this%fmtofs(dqprdt)


  end subroutine dqdt_d


  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  module subroutine ddt_d(this, psi, psid, psit, psitd, qprime, qprimed, for, dqprdt, dqprdtd)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),      intent( in) :: psi(:,:)    ! stream function at the nvl levels
    real(r8kind),      intent( in) :: psid(:,:)    ! stream function at the nvl levels
    real(r8kind),      intent( in) :: psit(:,:)   ! thickness at the ntl levels
    real(r8kind),      intent( in) :: psitd(:,:)   ! thickness at the ntl levels
    real(r8kind),      intent( in) :: qprime(:,:) ! potential vorticity
    real(r8kind),      intent( in) :: qprimed(:,:) ! potential vorticity
    real(r8kind),      intent( in) :: for(:,:)    ! constant potential vorticity forcing at the nvl levels
    real(r8kind),      intent(out) :: dqprdt(:,:)
    real(r8kind),      intent(out) :: dqprdtd(:,:)

    integer      :: k, l, i, j
    real(r8kind) :: dum1, dum2
    real(r8kind) :: dum1d, dum2d

    dqprdtd = 0.0_8

    ! advection of potential vorticity at upper level
    call this%jacob_d(psi(:, 1), psid(:, 1), qprime(:, 1), qprimed(:, 1), dqprdt(:, 1), dqprdtd(:, 1))

    ! advection of potential vorticity at middle level
    call this%jacob_d(psi(:, 2), psid(:, 2), qprime(:, 2), qprimed(:, 2), dqprdt(:, 2), dqprdtd(:, 2))

    ! advection of potential vorticity and dissipation at lower level
    call this%jacobd_d(psi(:, 3), psid(:, 3), qprime(:, 3), qprimed(:, 3), dqprdt(:, 3), dqprdtd(:, 3))

    ! relaxation of temperature and forcing
    do k = 1, this%nsh2
      dum1d = this%relt1 * psitd(k, 1)
      dum1 = this%relt1 * psit(k, 1)
      dum2d = this%relt2 * psitd(k, 2)
      dum2 = this%relt2 * psit(k, 2)
      dqprdtd(k, 1) = dqprdtd(k, 1) + dum1d
      dqprdt(k, 1) = dqprdt(k, 1) + dum1 + for(k, 1)
      dqprdtd(k, 2) = dqprdtd(k, 2) - dum1d + dum2d
      dqprdt(k, 2) = dqprdt(k, 2) - dum1 + dum2 + for(k, 2)
      dqprdtd(k, 3) = dqprdtd(k, 3) - dum2d
      dqprdt(k, 3) = dqprdt(k, 3) - dum2 + for(k, 3)
    end do

    ! explicit horizontal diffusion
    do l = 1, this%nvl
      do k = 1, this%nsh2
        dqprdtd(k, l) = dqprdtd(k, l) + this%diss(k, 1) * qprimed(k, l)
        dqprdt(k, l) = dqprdt(k, l) + this%diss(k, 1) * qprime(k, l)
      end do
    end do

  end subroutine ddt_d


  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  module subroutine jacob_d(this, psiloc, psilocd, pvor, pvord, sjacob, sjacobd)

    class(qg_tl_type), intent(in) :: this
    real(r8kind), intent( in) :: psiloc(:)
    real(r8kind), intent( in) :: psilocd(:)
    real(r8kind), intent( in) :: pvor(:)
    real(r8kind), intent( in) :: pvord(:)
    real(r8kind), intent(out) :: sjacob(:)
    real(r8kind), intent(out) :: sjacobd(:)

    integer      :: i, j, k
    real(r8kind) :: vv(this%nsh2)
    real(r8kind) :: vvd(this%nsh2)
    real(r8kind) :: dpsidl(this%nlat, this%nlon),  dpsidm(this%nlat, this%nlon),  dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidld(this%nlat, this%nlon),  dpsidmd(this%nlat, this%nlon),  dvordld(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon),  gjacob(this%nlat, this%nlon),  dpsidls(this%nsh2)
    real(r8kind) :: dvordmd(this%nlat, this%nlon),  gjacobd(this%nlat, this%nlon),  dpsidlsd(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vvd = reshape(ggsp%ddl(pvord), (/this%nsh2/))
    vv = reshape(ggsp%ddl(pvor), (/this%nsh2/))
    dvordld = ggsp%sptogg_pp(vvd)
    dvordl = ggsp%sptogg_pp(vv)
    dvordmd = ggsp%sptogg_pd(pvord)
    dvordm = ggsp%sptogg_pd(pvor)

    ! space derivatives of streamfunction
    dpsidlsd = reshape(ggsp%ddl(psilocd), (/this%nsh2/))
    dpsidls = reshape(ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidld = ggsp%sptogg_pp(dpsidlsd)
    dpsidl = ggsp%sptogg_pp(dpsidls)
    dpsidmd = ggsp%sptogg_pd(psilocd)
    dpsidm = ggsp%sptogg_pd(psiloc)

    gjacobd = 0.0_8

    ! jacobian term
    do j = 1, this%nlon
      do i = 1, this%nlat
        gjacobd(i, j) = dpsidmd(i, j) * dvordl(i, j) + dpsidm(i, j) * dvordld(i, j) - dpsidld(i, j) * dvordm(i, j) - dpsidl(i, j) * dvordmd(i, j)
        gjacob(i, j) = dpsidm(i, j) * dvordl(i, j) - dpsidl(i, j) * dvordm(i, j)
      end do
    end do
    sjacobd = reshape(ggsp%ggtosp(gjacobd), (/this%nsh2/))
    sjacob = reshape(ggsp%ggtosp(gjacob), (/this%nsh2/))

    ! planetary vorticity advection
    do k = 1, this%nsh2
      sjacobd(k) = sjacobd(k) - dpsidlsd(k)
      sjacob(k) = sjacob(k) - dpsidls(k)
    end do

  end subroutine jacob_d


  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  module subroutine jacobd_d(this, psiloc, psilocd, pvor, pvord, sjacob, sjacobd)

    class(qg_tl_type), intent( in) :: this
    real(r8kind),      intent( in) :: psiloc(:)
    real(r8kind),      intent( in) :: psilocd(:)
    real(r8kind),      intent( in) :: pvor(:)
    real(r8kind),      intent( in) :: pvord(:)
    real(r8kind),      intent(out) :: sjacob(:)
    real(r8kind),      intent(out) :: sjacobd(:)

    integer      :: i, j, k
    real(r8kind) :: dpsidl(this%nlat, this%nlon),  dpsidm(this%nlat, this%nlon),  dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidld(this%nlat, this%nlon),  dpsidmd(this%nlat, this%nlon),  dvordld(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon),  gjacob(this%nlat, this%nlon),  vv(this%nsh2)
    real(r8kind) :: dvordmd(this%nlat, this%nlon),  gjacobd(this%nlat, this%nlon),  vvd(this%nsh2)
    real(r8kind) :: azeta(this%nlat, this%nlon), dpsidls(this%nsh2)
    real(r8kind) :: azetad(this%nlat, this%nlon), dpsidlsd(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vvd = reshape(ggsp%ddl(pvord), (/this%nsh2/))
    vv = reshape(ggsp%ddl(pvor), (/this%nsh2/))
    dvordld = ggsp%sptogg_pp(vvd)
    dvordl = ggsp%sptogg_pp(vv)
    dvordmd = ggsp%sptogg_pd(pvord)
    dvordm = ggsp%sptogg_pd(pvor)

    ! space derivatives of streamfunction
    dpsidlsd = reshape(ggsp%ddl(psilocd), (/this%nsh2/))
    dpsidls = reshape(ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidld = ggsp%sptogg_pp(dpsidlsd)
    dpsidl = ggsp%sptogg_pp(dpsidls)
    dpsidmd = ggsp%sptogg_pd(psilocd)
    dpsidm = ggsp%sptogg_pd(psiloc)

    gjacobd = 0.0_8

    ! jacobian term + orographic forcing
    do j = 1, this%nlon
      do i = 1, this%nlat
        gjacobd(i, j) = dpsidmd(i, j) * (dvordl(i, j) + this%sinfi(i) * this%dorodl(i, j)) + dpsidm(i, j) * dvordld(i, j) - &
                     &  dpsidld(i, j) * (dvordm(i, j) + this%sinfi(i) * this%dorodm(i, j)) - dpsidl(i, j) * dvordmd(i, j)
        gjacob(i, j) = dpsidm(i, j) * (dvordl(i, j) + this%sinfi(i) * this%dorodl(i, j)) - dpsidl(i, j) * (dvordm(i, j) + &
                     & this%sinfi(i) * this%dorodm(i, j))
      end do
    end do

    ! dissipation 
    if (this%lgdiss) then

      !   spatially varying dissipation 
      do k = 1, this%nsh2
        vvd(k) = this%diss(k, 2) * psilocd(k)
        vv(k) = this%diss(k, 2) * psiloc(k)
      end do

      azetad = ggsp%sptogg_pp(vvd)
      azeta = ggsp%sptogg_pp(vv)

      do j = 1, this%nlon
        do i = 1, this%nlat
          gjacobd(i, j) = gjacobd(i, j) - this%ddisdy(i, j) * dpsidmd(i, j) - this%ddisdx(i, j) * dpsidld(i, j) + this%rdiss(i, j) * azetad(i, j)
          gjacob(i, j) = gjacob(i, j) - dpsidm(i, j) * this%ddisdy(i, j) - dpsidl(i, j) * this%ddisdx(i, j) + this%rdiss(i, j) * azeta(i, j)
        end do
      end do

      sjacobd = reshape(ggsp%ggtosp(gjacobd), (/this%nsh2/))
      sjacob = reshape(ggsp%ggtosp(gjacob), (/this%nsh2/))

    else

      !   uniform dissipation
      sjacobd = reshape(ggsp%ggtosp(gjacobd), (/this%nsh2/))
      sjacob = reshape(ggsp%ggtosp(gjacob), (/this%nsh2/))
      do k = 1, this%nsh2
        sjacobd(k) = sjacobd(k) + this%diss(k, 2) * psilocd(k)
        sjacob(k) = sjacob(k) + this%diss(k, 2) * psiloc(k)
      end do

    end if

    ! planetary vorticity advection
    do k = 1, this%nsh2
      sjacobd(k) = sjacobd(k) - dpsidlsd(k)
      sjacob(k) = sjacob(k) - dpsidls(k)
    end do

  end subroutine jacobd_d


  !-------------------------------------------------------------------------------
  ! get_trajectory
  !-------------------------------------------------------------------------------
  module function get_trajectory_tl(this) result(trajectory)

    class(qg_tl_type),            intent(in) :: this
!    real(r8kind), dimension(this%nsh2,this%nvl) :: trajectory
    real(r8kind), allocatable :: trajectory(:,:)

!    trajectory = this%trajectory
    allocate(trajectory, source = this%trajectory)

  end function get_trajectory_tl


  !-------------------------------------------------------------------------------
  ! get_qprimed
  !-------------------------------------------------------------------------------
  module function get_qprimed(this) result(qprimed)

    class(qg_tl_type),            intent(in) :: this
!    real(r8kind), dimension(this%nsh2,this%nvl) :: qprimed
    real(r8kind), allocatable :: qprimed(:,:)

!    qprimed = this%qprimed
    allocate(qprimed, source = this%qprimed)

  end function get_qprimed


end submodule QG_Model_TL
