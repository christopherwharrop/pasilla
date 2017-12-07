submodule(QG_Model) QG_Model_ADJ

  use kind
  use QG_Config
  use QG_GGSP
  use QG_Util

  implicit none

contains


  !-------------------------------------------------------------------------------
  ! constructor_qg_adj
  !-------------------------------------------------------------------------------
  module function constructor_qg_adj(config, state, state_vector, trajectory, trajectory_vector, for, step) result (qg_adj)

    type(qg_config_type),   intent(in) :: config
    real(r8kind), optional, intent(in) :: state(:,:)
    real(r8kind), optional, intent(in) :: state_vector(:)
    real(r8kind), optional, intent(in) :: trajectory(:,:)
    real(r8kind), optional, intent(in) :: trajectory_vector(:)
    real(r8kind), optional, intent(in) :: for(:,:)
    integer,      optional, intent(in) :: step
    type(qg_adj_type)                  :: qg_adj

    real(r8kind), allocatable          :: psig3d(:,:,:), psisp(:,:)
    real(r8kind), allocatable          :: trajg3d(:,:,:)

    real(r8kind) :: nsteps_per_day

    ! Physical constants
    real(r8kind), parameter :: radius = 6.37e+6
    real(r8kind), parameter :: om = 4d0 * pi / (24d0 * 3600d0)
    real(r8kind), parameter :: facsf = om * (radius)**2

    ! Set the model config
    qg_adj%config = config

    ! Get model grid dimensions
    qg_adj%ft = config%get_ft()
    qg_adj%nm = config%get_nm()
    qg_adj%nlon = config%get_nlon()
    qg_adj%nlat = config%get_nlat()
    qg_adj%nvl = config%get_nvl()
    qg_adj%ntl = config%get_ntl()
    qg_adj%nsh = config%get_nsh()
    qg_adj%nsh2 = config%get_nsh2()

    ! Get model grid converter
    qg_adj%ggsp = config%get_ggsp()

    ! Get model forcing from configuration
    allocate(qg_adj%for, source = config%get_for())

    ! Get Nondimensional relaxation coefficients
    qg_adj%relt1 = config%get_relt1()
    qg_adj%relt2 = config%get_relt2()

    ! Get dissipation coefficients for each spherical harmonic
    allocate(qg_adj%diss, source = config%get_diss())

    ! Get gauss points in radians and sine and cosine of it
    allocate(qg_adj%phi, source = config%get_phi())
    allocate(qg_adj%sinfi, source = config%get_sinfi())
    allocate(qg_adj%cosfi, source = config%get_cosfi())

    ! Get derivatives of orog
    allocate(qg_adj%dorodl, source = config%get_dorodl())
    allocate(qg_adj%dorodm, source = config%get_dorodm())

    ! Get orography and land-sea mask dependent friction option
    qg_adj%lgdiss = config%get_lgdiss()

    ! Get landsea-mask/orography dependent friction
    allocate(qg_adj%rdiss, source = config%get_rdiss())
    allocate(qg_adj%ddisdx, source = config%get_ddisdx())
    allocate(qg_adj%ddisdy, source = config%get_ddisdy())

    ! Get Laplace and Helmholtz operator for Q-PSI inversion
    allocate(qg_adj%rinhel(qg_adj%nsh2,0:5))
    qg_adj%rinhel = config%get_rinhel()

    ! Get one over Rossby rad. of def. squared of 200-500, and 500-800 thicknesses
    qg_adj%rl1 = config%get_rl1()
    qg_adj%rl2 = config%get_rl2()

    ! Initialize model step and clock
    if (present(step)) then
      qg_adj%step = step
    else
      qg_adj%step = 0
    end if
    qg_adj%clock = qg_adj%step * config%get_time_step()

    ! Initialize time step of the model:
    nsteps_per_day = 24.0d0 * 3600.0d0 / real(config%get_time_step())
    qg_adj%dtt = (1d0 / nsteps_per_day) * pi * 4d0

    ! Initialize streamfunction
    if (present(state)) then
      call qg_adj%init_state(psi=state)
    else if (present(state_vector)) then
      allocate(psig3d(qg_adj%nlon, qg_adj%nlat, qg_adj%nvl))
      allocate(psisp(qg_adj%nsh2, qg_adj%nvl))
      psig3d = reshape(state_vector,(/qg_adj%nlon, qg_adj%nlat, qg_adj%nvl/))
      psig3d(:,:,:) = psig3d(:,:,:) / facsf
      psisp(:,1) = reshape(qg_adj%ggsp%ggtosp(transpose(psig3d(:,:,1))), (/qg_adj%nsh2/))
      psisp(:,2) = reshape(qg_adj%ggsp%ggtosp(transpose(psig3d(:,:,2))), (/qg_adj%nsh2/))
      psisp(:,3) = reshape(qg_adj%ggsp%ggtosp(transpose(psig3d(:,:,3))), (/qg_adj%nsh2/))
      call qg_adj%init_state(psi=psisp)
    else
      call qg_adj%init_state()
    end if

    ! Initialize trajectory
    if (present(trajectory)) then
      allocate(qg_adj%trajectory, source = trajectory)
    else if (present(trajectory_vector)) then
      allocate(trajg3d(qg_adj%nlon, qg_adj%nlat, qg_adj%nvl))
      allocate(qg_adj%trajectory(qg_adj%nsh2,qg_adj%nvl))
      trajg3d = reshape(trajectory_vector,(/qg_adj%nlon, qg_adj%nlat, qg_adj%nvl/))
      trajg3d(:,:,:) = trajg3d(:,:,:) / facsf
      qg_adj%trajectory(:,1) = reshape(qg_adj%ggsp%ggtosp(transpose(trajg3d(:,:,1))), (/qg_adj%nsh2/))
      qg_adj%trajectory(:,2) = reshape(qg_adj%ggsp%ggtosp(transpose(trajg3d(:,:,2))), (/qg_adj%nsh2/))
      qg_adj%trajectory(:,3) = reshape(qg_adj%ggsp%ggtosp(transpose(trajg3d(:,:,3))), (/qg_adj%nsh2/))
    else
      allocate(qg_adj%trajectory, source = qg_adj%psi)
    end if
    allocate(qg_adj%psitb(qg_adj%nsh2,qg_adj%ntl))
    allocate(qg_adj%qprimeb(qg_adj%nsh2,qg_adj%nvl))
    call qg_adj%psitoq(qg_adj%trajectory, qg_adj%psitb, qg_adj%qprimeb)

  end function constructor_qg_adj


  !------------------------------------------------------------------
  ! destructor_qg_adj
  !
  ! Deallocates pointers used by a qg_adj_type object (none currently)
  !------------------------------------------------------------------
  module elemental subroutine destructor_qg_adj(this)

    type(qg_adj_type), intent(inout) :: this

    ! No pointers in qg_adj_type object so we do nothing

  end subroutine destructor_qg_adj


  !-----------------------------------------------------------------------
  ! performs a fourth order runge kutta time step at truncation nm
  ! with time step dt
  ! dqdt calculates the time derivative
  ! input  qprime at current time
  ! output qprime at current time plus dt
  !-----------------------------------------------------------------------
  module subroutine adv_nsteps_adj(this, nsteps)

    class(qg_adj_type), intent(inout) :: this
    integer           , intent(   in) :: nsteps

    integer :: step, k, l, nvar, i
    real(r8kind) :: dt2, dt6
    real(r8kind) :: y(this%nsh2, this%nvl), dydt(this%nsh2, this%nvl), yt(this%nsh2, this%nvl)
    real(r8kind) :: yt1(this%nsh2, this%nvl), yt2(this%nsh2, this%nvl), y1(this%nsh2, this%nvl)
    real(r8kind) :: dyt(this%nsh2, this%nvl), dym(this%nsh2, this%nvl)
    real(r8kind) :: yb(this%nsh2, this%nvl), dydtb(this%nsh2, this%nvl), ytb(this%nsh2, this%nvl)
    real(r8kind) :: dytb(this%nsh2, this%nvl), dymb(this%nsh2, this%nvl)
    real(r8kind) :: qprimeb_save(this%nsh2, this%nvl)
    real(r8kind) :: tempb
    type(qg_config_type) :: config
    type(qg_tl_type) :: qg_tl

    config = this%config

    if (nsteps > 0) then

      nvar = (this%nm + 2) * this%nm
      dt2 = this%dtt * 0.5d0
      dt6 = this%dtt / 6d0

      ! Advance the model forward in time n steps
      do step = 1, nsteps

        ! Initialize a QG TL model with the current state and trajectory
        qg_tl = qg_tl_type(this%config, state=this%psi, trajectory=this%trajectory)

        ! Advance the QG TL model one step
        call qg_tl%adv_nsteps(1)

        ! Save the current trajectory
        qprimeb_save = this%qprimeb

        ! Set the adjoint trajectory to the initial perturbation
        this%qprimeb = qg_tl%get_qprimed() - this%qprimeb

        ! Move the perturbation backward
        y = this%fmtofs(this%qprime)
        call this%dqdt(y, dydt)
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + dt2 * dydt(k, l)
          end do
        end do
        call this%dqdt(yt, dyt)
        yt1 = yt
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + dt2 * dyt(k, l)
          end do
        end do
        call this%dqdt(yt, dym)
        yt2 = yt
        do l = 1, this%nvl
          do k = 1, nvar
            yt(k, l) = y(k, l) + this%dtt * dym(k, l)
            dym(k, l) = dyt(k, l) + dym(k, l)
          end do
        end do
        call this%dqdt(yt, dyt)
        y1 = y
        do l = 1, this%nvl
          do k = 1, nvar
            y(k, l) = y(k, l) + dt6 * (dydt(k, l) + dyt(k, l) + 2. * dym(k, l))
          end do
        end do
        yb = 0.0_8

       call this%fstofm_b(y, yb, this%nm, this%qprimeb)
!        yb = this%fmtofs(this%qprimeb)
!        this%qprimeb = 0.0_r8kind

        dymb = 0.0_8
        dytb = 0.0_8
        dydtb = 0.0_8
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            tempb = dt6 * yb(k, l)
            dydtb(k, l) = dydtb(k, l) + tempb
            dytb(k, l) = dytb(k, l) + tempb
            dymb(k, l) = dymb(k, l) + 2. * tempb
          end do
        end do
        ytb = 0.0_8

        call this%dqdt_b(yt, ytb, dyt, dytb)
        dytb = 0.0_8
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            dytb(k, l) = dytb(k, l) + dymb(k, l)
            yb(k, l) = yb(k, l) + ytb(k, l)
            dymb(k, l) = dymb(k, l) + this%dtt * ytb(k, l)
            ytb(k, l) = 0.0_8
          end do
        end do
        call this%dqdt_b(yt2, ytb, dym, dymb)
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            yb(k, l) = yb(k, l) + ytb(k, l)
            dytb(k, l) = dytb(k, l) + dt2 * ytb(k, l)
            ytb(k, l) = 0.0_8
          end do
        end do
        call this%dqdt_b(yt1, ytb, dyt, dytb)
        do l = this%nvl, 1, -1
          do k = nvar, 1, -1
            yb(k, l) = yb(k, l) + ytb(k, l)
            dydtb(k, l) = dydtb(k, l) + dt2 * ytb(k, l)
            ytb(k, l) = 0.0_8
          end do
        end do
        call this%dqdt_b(y1, yb, dydt, dydtb)

        call this%fmtofs_b(this%qprime, this%qprimeb, yb)
!       this%qprimeb = this%fstofm(yb, this%nm)
!       yb = 0.0_r8kind

        ! Update the trajectory
        this%qprimeb = qprimeb_save - this%qprimeb

        ! Update model state with original trajectory plus its increment
        this%qprime = this%qprime - this%qprimeb

        ! Make stream function consistent with potential vorticity
        call this%qtopsi(this%qprimeb, this%trajectory, this%psitb)
        call this%qtopsi(this%qprime, this%psi, this%psit)

        ! Inrement the step count
        this%step = this%step - 1

      end do

    end if

  end subroutine adv_nsteps_adj


  !-----------------------------------------------------------------------
  !  Differentiation of dqdt in reverse (adjoint) mode:
  !   gradient     of useful results: tmp y dydt
  !   with respect to varying inputs: tmp y
  !-----------------------------------------------------------------------
  ! computation of time derivative of the potential vorticity field
  ! input  y potential vorticity in french format
  ! output dydt time derivative of y in french format
  ! values of qprime,  psi and psit are changed
  !-----------------------------------------------------------------------
  module subroutine dqdt_b(this, y, yb, dydt, dydtb)

    class(qg_adj_type), intent(inout) :: this
    real(r8kind),       intent(   in) :: y(:,:)
    real(r8kind)                      :: yb(:,:)
    real(r8kind)                      :: dydt(:,:)
    real(r8kind)                      :: dydtb(:,:)

    ! qprime
    real(r8kind) :: local_qprime(this%nsh2, this%nvl)
    real(r8kind) :: local_qprimeb(this%nsh2, this%nvl)
    ! psi
    real(r8kind) :: local_psi(this%nsh2, this%nvl)
    real(r8kind) :: local_psib(this%nsh2, this%nvl)
    ! psit
    real(r8kind) :: local_psit(this%nsh2, this%ntl)
    real(r8kind) :: local_psitb(this%nsh2, this%ntl)
    ! time derivative of qprime
    real(r8kind) :: dqprdt(this%nsh2, this%nvl)
    real(r8kind) :: dqprdtb(this%nsh2, this%nvl)

    local_psi = 0.0_r8kind
    local_psit = 0.0_r8kind
    local_qprime = this%fstofm(y, this%nm)

    call this%qtopsi(local_qprime, local_psi, local_psit)

    dqprdt = this%ddt(local_psi, local_psit, local_qprime, this%for)

    call this%fmtofs_b(dqprdt, dqprdtb, dydtb)

    call this%ddt_b(local_psi, local_psib, local_psit, local_psitb, local_qprime, local_qprimeb, this%for, dqprdtb)

    local_psi = 0.0_r8kind
    local_psit = 0.0_r8kind
    call this%QTOPSI_B(local_qprime, local_qprimeb, local_psi, local_psib, local_psit, local_psitb)

    call this%fstofm_b(y, yb, this%nm, local_qprimeb)

  end subroutine dqdt_b


  !----------------------------------------------------------------------
  !  Differentiation of ddt in reverse (adjoint) mode:
  !   gradient     of useful results: tmp dqprdt
  !   with respect to varying inputs: tmp psi qprime psit
  !----------------------------------------------------------------------
  ! ddt
  !
  ! computation of time derivative of the potential vorticity fields
  !
  ! input qprime,  psi,  psit
  ! output dqprdt
  !----------------------------------------------------------------------
  module subroutine ddt_b(this, psi, psib, psit, psitb, qprime, qprimeb, for, dqprdtb)
    class(qg_adj_type), intent(inout) :: this    
    real(r8kind),       intent(   in) :: psi(:,:)
    real(r8kind)                      :: psib(:,:)
    real(r8kind),       intent(   in) :: psit(:,:)
    real(r8kind)                      :: psitb(:,:)
    real(r8kind),       intent(   in) :: qprime(:,:)
    real(r8kind)                      :: qprimeb(:,:)
    real(r8kind),       intent(   in) :: for(:,:)
    real(r8kind)                      :: dqprdtb(:,:)

    integer :: k, l, i, j
    real(r8kind) :: dum1, dum2
    real(r8kind) :: dum1b, dum2b
    real(r8kind), dimension(this%nsh2) :: res
    real(r8kind), dimension(this%nsh2) :: resb
    real(r8kind), dimension(this%nsh2) :: res0
    real(r8kind), dimension(this%nsh2) :: resb0

    ! advection of potential vorticity at upper level
    res = this%jacob(psi(:, 1), qprime(:, 1))

    ! advection of potential vorticity at middle level
    res0 = this%jacob(psi(:, 2), qprime(:, 2))

    ! advection of potential vorticity and dissipation at lower level
    qprimeb = 0.0_8
    do l=3,1,-1
      do k=this%nsh2,1,-1
        qprimeb(k, l) = qprimeb(k, l) + this%diss(k, 1) * dqprdtb(k, l)
      end do
    end do
    psitb = 0.0_8
    do k=this%nsh2,1,-1
      dum2b = dqprdtb(k, 2) - dqprdtb(k, 3)
      dum1b = dqprdtb(k, 1) - dqprdtb(k, 2)
      psitb(k, 2) = psitb(k, 2) + this%relt2 * dum2b
      psitb(k, 1) = psitb(k, 1) + this%relt1 * dum1b
    end do
    call this%jacobd_b(psi(:, 3), psib(:, 3), qprime(:, 3), qprimeb(:, 3), dqprdtb(:, 3))
    dqprdtb(:, 3) = 0.0_8
    resb0 = dqprdtb(:, 2)
    dqprdtb(:, 2) = 0.0_8

    call this%jacob_b(psi(:, 2), psib(:, 2), qprime(:, 2), qprimeb(:, 2), resb0)
    resb = dqprdtb(:, 1)

    call this%jacob_b(psi(:, 1), psib(:, 1), qprime(:, 1), qprimeb(:, 1), resb)

  end subroutine ddt_b


  !----------------------------------------------------------------------
  !  Differentiation of jacob in reverse (adjoint) mode:
  !   gradient     of useful results: tmp sjacob psiloc pvor
  !   with respect to varying inputs: tmp psiloc pvor
  !----------------------------------------------------------------------
  ! advection of potential vorticity
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  module subroutine jacob_b(this, psiloc, psilocb, pvor, pvorb, sjacobb)

    class(qg_adj_type), intent(in) :: this
    real(r8kind),       intent(in) :: psiloc(:)
    real(r8kind)                   :: psilocb(:)
    real(r8kind),       intent(in) :: pvor(:)
    real(r8kind)                   :: pvorb(:)
    real(r8kind)                   :: sjacobb(:)

    integer :: i, j, k
    real(r8kind) :: vv(this%nsh2)
    real(r8kind) :: vvb(this%nsh2)
    real(r8kind) :: dpsidl(this%nlat, this%nlon), dpsidm(this%nlat, this%nlon), dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidlb(this%nlat, this%nlon), dpsidmb(this%nlat, this%nlon), dvordlb(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon), gjacob(this%nlat, this%nlon), dpsidls(this%nsh2)
    real(r8kind) :: dvordmb(this%nlat, this%nlon), gjacobb(this%nlat, this%nlon), dpsidlsb(this%nsh2)
    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity
    vv = reshape(ggsp%ddl(pvor), (/this%nsh2/))
    dvordl = ggsp%sptogg_pp(vv)
    dvordm = ggsp%sptogg_pd(pvor)

    ! space derivatives of streamfunction
    dpsidls = reshape(ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidl = ggsp%sptogg_pp(dpsidls)
    dpsidm = ggsp%sptogg_pd(psiloc)
    dpsidlsb = 0.0_8
    do k = this%nsh2, 1, -1
      dpsidlsb(k) = dpsidlsb(k) - sjacobb(k)
    end do

    gjacobb = ggsp%sptogg_pw(sjacobb)

    dpsidlb = 0.0_8
    dpsidmb = 0.0_8
    dvordlb = 0.0_8
    dvordmb = 0.0_8
    do j = this%nlon, 1, -1
      do i = this%nlat, 1, -1
        dpsidmb(i, j) = dpsidmb(i, j) + dvordl(i, j) * gjacobb(i, j)
        dvordlb(i, j) = dvordlb(i, j) + dpsidm(i, j) * gjacobb(i, j)
        dpsidlb(i, j) = dpsidlb(i, j) - dvordm(i, j) * gjacobb(i, j)
        dvordmb(i, j) = dvordmb(i, j) - dpsidl(i, j) * gjacobb(i, j)
        gjacobb(i, j) = 0.0_8
      end do
    end do

    psilocb = reshape(ggsp%ggtosp(dpsidmb), (/this%nsh2/))
    dpsidlsb = reshape(ggsp%ggtosp(dpsidlb), (/this%nsh2/))

    call ggsp%ddl_b(psiloc, psilocb, dpsidlsb)

    pvorb = reshape(ggsp%ggtosp(dvordmb), (/this%nsh2/))

    vvb = 0.0_8
    vvb = reshape(ggsp%ggtosp(dvordlb), (/this%nsh2/))

    call ggsp%ddl_b(pvor, pvorb, vvb)

  end subroutine jacob_b


  !----------------------------------------------------------------------
  !  Differentiation of jacobd in reverse (adjoint) mode:
  !   gradient     of useful results: tmp sjacob pvor
  !   with respect to varying inputs: tmp psiloc pvor
  !----------------------------------------------------------------------
  ! advection of potential vorticity and dissipation on gaussian grid
  ! input psiloc,  pvor
  ! output sjacob
  !----------------------------------------------------------------------
  module subroutine jacobd_b(this, psiloc, psilocb, pvor, pvorb, sjacobb)

    class(qg_adj_type), intent(inout) :: this
    real(r8kind),       intent(   in) :: psiloc(:)
    real(r8kind),       intent(  out) :: psilocb(:)
    real(r8kind),       intent(   in) :: pvor(:)
    real(r8kind),       intent(  out) :: pvorb(:)
    real(r8kind),       intent(   in) :: sjacobb(:)

    integer :: i, j, k
    real(r8kind) :: dpsidl(this%nlat, this%nlon), dpsidm(this%nlat, this%nlon), dvordl(this%nlat, this%nlon)
    real(r8kind) :: dpsidlb(this%nlat, this%nlon), dpsidmb(this%nlat, this%nlon), dvordlb(this%nlat, this%nlon)
    real(r8kind) :: dvordm(this%nlat, this%nlon), gjacob(this%nlat, this%nlon), vv(this%nsh2)
    real(r8kind) :: dvordmb(this%nlat, this%nlon), gjacobb(this%nlat, this%nlon), vvb(this%nsh2)
    real(r8kind) :: azeta(this%nlat, this%nlon), dpsidls(this%nsh2)
    real(r8kind) :: azetab(this%nlat, this%nlon), dpsidlsb(this%nsh2)
    real(r8kind) :: tempb, tempb0

    type(qg_ggsp_type) :: ggsp

    ! Get grid conversion object
    ggsp = this%ggsp

    ! space derivatives of potential vorticity 
    vv = reshape(ggsp%ddl(pvor), (/this%nsh2/))
    dvordl = ggsp%sptogg_pp(vv)
    dvordm = ggsp%sptogg_pd(pvor)

    ! space derivatives of streamfunction
    dpsidls = reshape(ggsp%ddl(psiloc), (/this%nsh2/))
    dpsidl = ggsp%sptogg_pp(dpsidls)
    dpsidm = ggsp%sptogg_pd(psiloc)

    ! dissipation 
    dpsidlsb = 0.0_r8kind
    do k = this%nsh2, 1, -1
      dpsidlsb(k) = dpsidlsb(k) - sjacobb(k)
    end do

    if (this%lgdiss) then

      gjacobb = ggsp%sptogg_pw(sjacobb)
      dpsidlb = 0.0_r8kind
      dpsidmb = 0.0_r8kind
      azetab = 0.0_r8kind
      do j = this%nlon, 1, -1
        do i = this%nlat, 1, -1
          dpsidmb(i, j) = dpsidmb(i, j) - this%ddisdy(i, j)*gjacobb(i, j)
          azetab(i, j) = azetab(i, j) + this%rdiss(i, j)*gjacobb(i, j)
          dpsidlb(i, j) = dpsidlb(i, j) - this%ddisdx(i, j)*gjacobb(i, j)
        end do
      end do
      vvb = 0.0_r8kind
      vvb = reshape(ggsp%ggtosp(azetab), (/this%nsh2/))
      psilocb = 0.0_r8kind
      do k = this%nsh2, 1, -1
        psilocb(k) = psilocb(k) + this%diss(k, 2) * vvb(k)
        vvb(k) = 0.0_r8kind
      end do

    else

      psilocb = 0.0_r8kind
      do k = this%nsh2, 1, -1
        psilocb(k) = psilocb(k) + this%diss(k, 2) * sjacobb(k)
      end do
      gjacobb = ggsp%sptogg_pw(sjacobb)
      dpsidlb = 0.0_r8kind
      dpsidmb = 0.0_r8kind
      vvb = 0.0_r8kind

    end if

    dvordlb = 0.0_r8kind
    dvordmb = 0.0_r8kind
    do j = this%nlon, 1, -1
      do i = this%nlat, 1, -1
        dpsidmb(i, j) = dpsidmb(i, j) + (this%sinfi(i) * this%dorodl(i, j) + dvordl(i, j)) * gjacobb(i, j)
        dvordlb(i, j) = dvordlb(i, j) + dpsidm(i, j) * gjacobb(i, j)
        dpsidlb(i, j) = dpsidlb(i, j) - (this%sinfi(i) * this%dorodm(i, j) + dvordm(i, j)) * gjacobb(i, j)
        dvordmb(i, j) = dvordmb(i, j) - dpsidl(i, j) * gjacobb(i, j)
        gjacobb(i, j) = 0.0_8
      end do
    end do

    psilocb = reshape(ggsp%ggtosp(dpsidmb), (/this%nsh2/))
    dpsidlsb = reshape(ggsp%ggtosp(dpsidlb), (/this%nsh2/))

    call ggsp%ddl_b(psiloc, psilocb, dpsidlsb)

    pvorb = reshape(ggsp%ggtosp(dvordmb), (/this%nsh2/))
    vvb = reshape(ggsp%ggtosp(dvordlb), (/this%nsh2/))

    call ggsp%ddl_b(pvor, pvorb, vvb)

  end subroutine jacobd_b


  !-----------------------------------------------------------------------
  !  Differentiation of qtopsi in reverse (adjoint) mode:
  !   gradient     of useful results: psi psit qprime
  !   with respect to varying inputs: qprime
  !-----------------------------------------------------------------------
  ! computation of streamfunction from potential vorticity
  ! input  qprime which is potential vorticity field
  ! output psi,  the streamfunction and psit,  the layer thicknesses
  !-----------------------------------------------------------------------
  module subroutine qtopsi_b(this, qprime, qprimeb, psi, psib, psit, psitb)

    class(qg_adj_type), intent(in) :: this
    real(r8kind),       intent(in) :: qprime(:,:)   ! potential vorticity
    real(r8kind)                   :: qprimeb(:,:)
    real(r8kind)                   :: psi(:,:)      ! stream function at the nvl levels
    real(r8kind)                   :: psib(:,:)
    real(r8kind)                   :: psit(:,:)     ! thickness at the ntl levels
    real(r8kind)                   :: psitb(:,:)

    integer :: k
    real(r8kind) :: r3
    real(r8kind) :: ws(this%nsh2)    ! only used as portable workspace
    real(r8kind) :: wsb(this%nsh2)
    real(r8kind) :: tempb
    real(r8kind) :: tempb0
    integer :: ad_to
    integer :: ad_to0
    integer :: ad_to1

    do k=1,size(psi, 1)
      ws(k) = qprime(k, 1) + qprime(k, 3)
      psi(k, 1) = this%rinhel(k, 1)*(ws(k)+qprime(k, 2))
      psi(k, 2) = ws(k) - 2.d0*qprime(k, 2)
      psi(k, 3) = qprime(k, 1) - qprime(k, 3)
    end do
    do k=1,size(psit, 1)
      psit(k, 1) = this%rinhel(k, 2)*psi(k, 2) + this%rinhel(k, 3)*psi(k, 3)
      psit(k, 2) = this%rinhel(k, 4)*psi(k, 2) + this%rinhel(k, 5)*psi(k, 3)
    end do
    r3 = 1./3.
    do k=1,size(psi, 1)
      psi(k, 2) = r3*(psi(k, 1)-psit(k, 1)+psit(k, 2))
      psi(k, 1) = psi(k, 2) + psit(k, 1)
      psi(k, 3) = psi(k, 2) - psit(k, 2)
    end do
    do k=size(psi, 1),1,-1
      psib(k, 2) = psib(k, 2) + psib(k, 3)
      psitb(k, 2) = psitb(k, 2) - psib(k, 3)
      psib(k, 3) = 0.0_r8kind
      psib(k, 2) = psib(k, 2) + psib(k, 1)
      psitb(k, 1) = psitb(k, 1) + psib(k, 1)
      psib(k, 1) = 0.0_r8kind
      tempb0 = r3*psib(k, 2)
      psib(k, 1) = psib(k, 1) + tempb0
      psitb(k, 1) = psitb(k, 1) - tempb0
      psitb(k, 2) = psitb(k, 2) + tempb0
      psib(k, 2) = 0.0_r8kind
    end do
    do k=size(psit, 1),1,-1
      psib(k, 2) = psib(k, 2) + this%rinhel(k, 4)*psitb(k, 2)
      psib(k, 3) = psib(k, 3) + this%rinhel(k, 5)*psitb(k, 2)
      psitb(k, 2) = 0.0_r8kind
      psib(k, 2) = psib(k, 2) + this%rinhel(k, 2)*psitb(k, 1)
      psib(k, 3) = psib(k, 3) + this%rinhel(k, 3)*psitb(k, 1)
      psitb(k, 1) = 0.0_r8kind
    end do
    wsb = 0.0_r8kind
    do k=size(psi, 1),1,-1
      qprimeb(k, 1) = qprimeb(k, 1) + psib(k, 3)
      qprimeb(k, 3) = qprimeb(k, 3) - psib(k, 3)
      psib(k, 3) = 0.0_r8kind
      wsb(k) = wsb(k) + psib(k, 2)
      qprimeb(k, 2) = qprimeb(k, 2) - 2.d0*psib(k, 2)
      psib(k, 2) = 0.0_r8kind
      tempb = this%rinhel(k, 1)*psib(k, 1)
      wsb(k) = wsb(k) + tempb
      qprimeb(k, 2) = qprimeb(k, 2) + tempb
      psib(k, 1) = 0.0_r8kind
      qprimeb(k, 1) = qprimeb(k, 1) + wsb(k)
      qprimeb(k, 3) = qprimeb(k, 3) + wsb(k)
      wsb(k) = 0.0_r8kind
    end do

  end subroutine qtopsi_b


  !-----------------------------------------------------------------------
  !  Differentiation of fmtofs in reverse (adjoint) mode:
  !   gradient     of useful results: z
  !   with respect to varying inputs: y
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
  module subroutine fmtofs_b(this, y, yb, zb)

    class(qg_adj_type)       :: this
    real(r8kind), intent(in) :: y(:,:)
    real(r8kind)             :: yb(:,:)
    real(r8kind)             :: zb(:,:)

    integer :: m, n, k, indx, l
    integer, dimension(size(y, 2)) :: k2

    integer :: max1
    integer :: branch
    integer :: ad_from
    integer :: ad_to

    k = 0
    do l=1,size(y, 2)
      k2(l) = k
      k = 1
      do m=0,this%nm
        if (m < 1) then
          max1 = 1
        else
          max1 = m
        end if
        ad_from = max1
        do n=ad_from, this%nm
          k = k + 1
          if (m .eq. 0) then
            indx = n**2
          else
            indx = n**2 + 2*m - 1
          end if
        end do
      end do
    end do
    yb = 0.0_r8kind
    do l=size(y, 2),1,-1
      do m= this%nm,0,-1
        if (m < 1) then
          max1 = 1
        else
          max1 = m
        end if
        ad_from = max1
        do n= this%nm,ad_from,-1
          if (m .eq. 0) then
            indx = n**2
          else
            indx = n**2 + 2*m - 1
          end if
          if (m .ne. 0) then
            yb(k+this%nsh, l) = yb(k+this%nsh, l) + zb(indx+1, l)
            zb(indx+1, l) = 0.0_r8kind
          end if
          yb(k, l) = yb(k, l) + zb(indx, l)
          zb(indx, l) = 0.0_r8kind
          k = k - 1
        end do
      end do
      k = k2(l)
    end do

  end subroutine fmtofs_b


  !-----------------------------------------------------------------------
  !  Differentiation of fstofm in reverse (adjoint) mode:
  !   gradient     of useful results: y z
  !   with respect to varying inputs: y
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
  module subroutine fstofm_b(this, y, yb, ntr, zb)

    class(qg_adj_type), intent(in) :: this
    real(r8kind),       intent(in) :: y(:,:)
    real(r8kind)                   :: yb(:,:)
    integer,            intent(in) :: ntr
    real(r8kind)                   :: zb(:,:)

    integer :: m, n, k, indx, i, l
    integer, dimension(size(y, 2)) :: k2

    integer :: max1
    integer :: ad_to
    integer :: branch
    integer :: ad_from
    integer :: ad_to0

    k = 0
    do l=1,size(y, 2)

      k2(l) = k
      k = 1
      do m=0,this%nm
        if (m .lt. 1) then
          max1 = 1
        else
          max1 = m
        end if
        ad_from = max1
        do n=ad_from,this%nm
          k = k + 1
          if (m .le. ntr .and. n .le. ntr) then
            if (m .eq. 0) then
              indx = n**2
            else
              indx = n**2 + 2*m - 1
            end if
          end if
        end do
      end do
    end do
    ad_to0 = l - 1
    do l=ad_to0,1,-1
      do m=this%nm,0,-1
        if (m .lt. 1) then
          max1 = 1
        else
          max1 = m
        end if
        ad_from = max1
        do n=this%nm,ad_from,-1
          if (m .le. ntr .and. n .le. ntr) then
            if (m .eq. 0) then
              indx = n**2
            else
              indx = n**2 + 2*m - 1
            end if
            if (m .ne. 0) then
              yb(indx+1, l) = yb(indx+1, l) + zb(k+this%nsh, l)
              zb(k+this%nsh, l) = 0.0_r8kind
            end if
            yb(indx, l) = yb(indx, l) + zb(k, l)
            zb(k, l) = 0.0_r8kind
          end if
          k = k - 1
        end do
      end do

      k = k2(l)

    end do

  end subroutine fstofm_b


  !-------------------------------------------------------------------------------
  ! get_trajectory
  !-------------------------------------------------------------------------------
  module function get_trajectory_adj(this) result(trajectory)

    class(qg_adj_type), intent(in) :: this
    real(r8kind), allocatable      :: trajectory(:,:)

    allocate(trajectory, source = this%trajectory)

  end function get_trajectory_adj


end submodule QG_Model_ADJ
