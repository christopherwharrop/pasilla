module solver

  use kind, only                   : r8kind
  use module_constants, only       : PI
  use config, only                 : Config_Type
  use background, only             : Background_Type
  use background_covariance, only  : Background_Covariance_Type
  use observation_covariance, only : Observation_Covariance_Type
  use innovation_vector, only      : Innovation_Vector_Type
  use observation_operator, only   : Observation_Operator_Type
  use lorenz96, only               : lorenz96_type, lorenz96_TL_type, lorenz96_ADJ_type

  implicit none

  private

  public :: Solver_Type

  type Solver_Type
    private
    ! instance variables
    integer, public           :: method
    real(r8kind)              :: alpha
    real(r8kind), allocatable :: hrh_cov(:,:,:)
    real(r8kind), allocatable :: brh_cov(:,:,:)
    real(r8kind), allocatable :: anl_vec(:,:)
    real(r8kind), allocatable :: bht_ino(:,:,:)
    real(r8kind), allocatable :: htr_ino(:,:,:)
    real(r8kind)              :: jvc_for(1,1)

  contains
    ! methods
    procedure, private :: pre_con_dif
    procedure, private :: pre_sol
    procedure          :: solve
    procedure          :: put_anl_vec
    final :: destructor
  end type Solver_Type

  interface Solver_Type
    procedure constructor
  end interface

contains

  !------------------------------------------------------------------
  ! constructor
  !
  ! Returns an initialized Solver
  !------------------------------------------------------------------
  type(Solver_Type) function constructor(cfg)

    class(Config_Type),     intent(in) :: cfg

    constructor%method = cfg%method
    constructor%alpha = cfg%alpha

  end function


  !------------------------------------------------------------------
  ! destructor
  !
  ! Deallocates pointers used by a Solver object (none currently)
  !------------------------------------------------------------------
  elemental subroutine destructor(this)

    type(Solver_Type), intent(inout) :: this

    ! No pointers in Solver object so we do nothing

  end subroutine

  ! BJE
  ! GENERATE PRE-CONDITIONED VECTOR V FROM X
  ! THIS REPLACED THE NEED TO EXPLICITLY INVERT B
  subroutine pre_con_dif(this, bkg_cov, non_vec)

    class(Solver_Type), intent(inout) :: this
    real(r8kind),       intent(   in) :: bkg_cov(:,:)
    real(r8kind),       intent(inout) :: non_vec(:,:)

    real(r8kind), allocatable :: con_vec(:,:)
    real(r8kind), allocatable :: bkg_cpy(:,:)
    integer                   :: i, j, jj, rad, info 
    integer,      allocatable :: ipiv(:)
    real(r8kind)              :: var
    integer                   :: bkg_len

    print *,"PRE_CON_DIF"

    bkg_len = size(bkg_cov, 2)

    allocate (ipiv(bkg_len))
    allocate (con_vec(bkg_len, 1))
    allocate (bkg_cpy(bkg_len, bkg_len))
    con_vec = non_vec
    bkg_cpy = bkg_cov

    call dgesv(bkg_len, 1, bkg_cpy, bkg_len, ipiv, con_vec, bkg_len, info)
    non_vec = con_vec

    print *,"PRE_CON_DIF COMPLETE"

  end subroutine pre_con_dif


  ! BJE
  ! PREPARES MATRICES FOR USE IN THE SOLVER
  ! THERE ARE A NUMBER OF MATRICES THAT ARE REUSED IN THE SOLVER
  ! SPECIFICLY:  HRH_COV=       (H(T)R(-1)H)
  !              HTR_INO=       (H(T)R(-1))*(Y-HXb)
  !              BRH_COV=B(1/2)*(H(T)R(-1)H)
  !              BHT_INO=B(1/2)*(H(T)R(-1))*(Y-HXb)

  subroutine pre_sol(this, obs_opr, obs_cov, bkg_cov, inno_vec)

    class(Solver_Type),                 intent(inout) :: this
    class(Observation_Operator_Type),   intent(   in) :: obs_opr
    class(Observation_Covariance_Type), intent(   in) :: obs_cov
    class(Background_Covariance_Type),  intent(   in) :: bkg_cov
    class(Innovation_Vector_Type),      intent(   in) :: inno_vec

    integer :: i, j, t
    integer :: tim_len, bkg_len, obs_len

    real(r8kind), allocatable :: tmp_mat(:,:)
    real(r8kind), allocatable :: tmp_vec(:,:)

    real(r8kind), allocatable :: tim_bkc(:,:)
    real(r8kind), allocatable :: tim_obc(:,:)
    real(r8kind), allocatable :: tim_hrh(:,:)
    real(r8kind), allocatable :: tim_opr(:,:)
    real(r8kind), allocatable :: tim_htr(:,:)

    real(r8kind), allocatable :: obs_vvc(:,:)
    real(r8kind), allocatable :: tmp_jfo(:,:)
    real(r8kind), allocatable :: obs_opt(:,:)
    real(r8kind), allocatable :: tmp_rhh(:,:)
    real(r8kind), allocatable :: tmp_hrr(:,:)

    print *, "PRE_SOLVER"

    tim_len = size(obs_cov%covariance, 1)
    obs_len = size(obs_cov%covariance, 2)
    bkg_len = size(bkg_cov%covariance, 2)

    allocate (tmp_mat(bkg_len, bkg_len))
    allocate (tmp_vec(bkg_len,       1))

    allocate (tim_bkc(bkg_len, bkg_len))
    allocate (tim_obc(obs_len, obs_len))
    allocate (tim_hrh(bkg_len, bkg_len))
    allocate (tim_opr(obs_len, bkg_len))
    allocate (tim_htr(bkg_len,       1))

    allocate (obs_vvc(obs_len,       1))
    allocate (tmp_jfo(obs_len,       1))
    allocate (obs_opt(bkg_len, obs_len))
    allocate (tmp_rhh(obs_len, bkg_len))
    allocate (tmp_hrr(bkg_len, obs_len))

    this%jvc_for(1,1) = 0.0
    do t = 1, tim_len
       ! ASSUME THAT OBS_OPR=H, OBS_COV=R(-1/2), BKG_COV=B(1/2), OBS_VEC=(Y-HXb)
       tim_opr(:,:) = obs_opr%operator(t,:,:)
       tim_obc(:,:) = obs_cov%covariance(t,:,:)
       tim_bkc(:,:) = bkg_cov%covariance(t,:,:)
       obs_vvc(:,1) = inno_vec%value(:)

       ! CREATE THE OBS BASED MATRICES, FOR USE IN CALCULATING THE COST FUNCTION
       ! tim_hrh=H(T)R(-1)H, tim_htr=R(-1)H
       ! R(-1/2)H
       call dgemm("N", "N", obs_len, bkg_len, obs_len, 1.d0, tim_obc, obs_len, tim_opr, obs_len, 0.d0, tmp_rhh, obs_len)
       obs_opt = transpose(tmp_rhh)            ! H(T)R(-1/2)

       ! H(T)R(-1)H
       call dgemm("N", "N", bkg_len, bkg_len, obs_len, 1.d0, obs_opt, bkg_len, tmp_rhh, obs_len, 0.d0, tim_hrh, bkg_len)
       this%hrh_cov(t,:,:) = tim_hrh(:,:)

      ! CREATE COST FUNCTION TERM (Y-HXb)R(-1)(Y-HXb), A CONSTANT
       call dgemv("N", obs_len, obs_len, 1.d0, tim_obc, obs_len, obs_vvc, 1, 0.d0, tmp_jfo, 1)
       do i = 1, obs_len                        ! CREATE (Y-HXb)(T)R(-1)(Y-HXb)
          this%jvc_for(1,1) = this%jvc_for(1,1) + tmp_jfo(i,1) * tmp_jfo(i,1)
       end do

       ! CREATE R(-1) from R(-1/2)
       call dgemm("N", "N", obs_len, obs_len, obs_len, 1.d0, tim_obc, obs_len, tim_obc, obs_len, 0.d0, tim_obc, obs_len)

       ! H(T)R(-1)
       call dgemm("N", "N", bkg_len, obs_len, obs_len, 1.d0, obs_opt, bkg_len, tim_obc, obs_len, 0.d0, tmp_hrr, bkg_len)

       ! H(T)R(-1)(Y-HXb)
       call dgemv("N", bkg_len, obs_len, 1.d0, tmp_hrr, bkg_len, obs_vvc, 1, 0.d0, tim_htr, 1)
       this%htr_ino(t,:,1) = tim_htr(:,1)

       ! CREATE THE UNPRECONDITIONED MATRICES, FOR THE GRADIENT OF J
       ! B(1/2)*H(T)R(-1)
       call dgemv("N", bkg_len, bkg_len, 1.d0, tim_bkc, bkg_len, tim_htr, 1, 0.d0, tmp_vec, 1)

       ! B(1/2)*H(T)R(-1)H
       call dgemm("N", "N", bkg_len, bkg_len, bkg_len, 1.d0, tim_bkc, bkg_len, tim_hrh, bkg_len, 0.d0, tmp_mat, bkg_len)

       this%bht_ino(t,:,1) = tmp_vec(:,1)
       this%brh_cov(t,:,:) = tmp_mat(:,:)
    end do

    print *,"PRE_SOLVER COMPLETE"

  end subroutine pre_sol


  ! BJE
  ! THE SOLVER FOR THE VARIATIONAL ASSIMILATION
  !  J    = (1/2)*(X- X1)(T)*(    D(-1) )*(X -X1)
  !       + (1/2)*(X- Xb)(T)*(H(T)R(-1)H)*(X -Xb)
  !       + (1/2)*(Y-HXb)(T)*(    R(-1)H)*(X- Xb)
  !       + (1/2)*(X- Xb)(T)*(H(T)R(-1) )*(Y-HXb)
  !       + (1/2)*(Y-HXb)(T)*(    R(-1) )*(Y-HXb)
  !
  ! GRADJ =      D(-1)  *(X- X1)
  !       + (H(T)R(-1)H)*(X- Xb)
  !       - (H(T)R(-1) )*(y-HXb)
  !
  ! WHERE D(1/2)V= (X- X1)
  ! D(1)=B, D(>1)=Q
  !
  ! THUS, FOR T=1, B(1/2)V=(X-Xb)
  !       FOR T>1, Q(1/2)V=(X-MX)    - FORWARD
  !                Q(1/2)V=(X-M(T)X) - BACKWARD
  !   B(1/2)Q(-1/2)Q(1/2)V=B(1/2)Q(-1/2)(X-MX)
  !                B(1/2)V=B(1/2)Q(-1/2)(X-MX)
  !
  !   ASSUME Q=B/C, WHERE C IS A CONSTANT GREATER THAN ONE
  !                B(1/2)V=C(X-MX)
  !   IF Q=0, THEN X=MX
  !
  !   THUS V=D(-1/2)(X2-X1), NO MATTER WHAT X1 AND X2 ARE
  !
  ! WHICH CAN BE REWRITTEN AS:
  !
  ! BRH=D(1/2)L(-T)(H(T)R(-1)H), BHT=(D(1/2)L(-T)H(T)R(-1))*(Y-HXb)
  ! WHERE D(1/2)V=(X-X1)
  !
  !  J    = (1/2)* V(T) * I *  V
  !       + (1/2)*(X-Xb)*HRH*(X-Xb)
  !       +              HTR*(X-Xb)
  !       + CONSTANT TERM
  !
  ! GRADJ = D(1/2)* [   V
  !                  + BRH*(X-Xb)
  !                  - BHT        ]
  subroutine solve(this, bkg, bkg_cov, obs_cov, obs_opr, inno_vec)

    class(Solver_Type),                 intent(inout) :: this
    class(Background_Type),             intent(   in) :: bkg
    class(Background_Covariance_Type),  intent(   in) :: bkg_cov
    class(Observation_Covariance_Type), intent(   in) :: obs_cov
    class(Observation_Operator_Type),   intent(   in) :: obs_opr
    class(Innovation_Vector_Type),      intent(   in) :: inno_vec

    real(r8kind), allocatable   :: tim_htr(:,:)
    real(r8kind), allocatable   :: tim_bkc(:,:)
    real(r8kind), allocatable   :: tim_hrh(:,:)
    real(r8kind), allocatable   :: tim_anv(:,:)
    real(r8kind), allocatable   :: tim_bkv(:,:)

    real(r8kind), allocatable   :: new_vec(:,:)
    real(r8kind), allocatable   :: tlm_vec(:,:)
    real(r8kind), allocatable   :: mdl_vec(:,:)
    real(r8kind), allocatable   :: ges_vec(:,:)
    real(r8kind), allocatable   :: dif_vec(:,:)
    real(r8kind), allocatable   :: dif_tra(:,:)

    real(r8kind), allocatable   :: pre_tra(:,:)
    real(r8kind), allocatable   :: pre_dif(:,:)
    real(r8kind), allocatable   :: tmp_mat(:,:)
    real(r8kind), allocatable   :: tmp_vec(:,:)
    real(r8kind), allocatable   :: tmp_vvc(:,:)

    real(r8kind), allocatable   :: grd_jvc(:,:)
    real(r8kind)                :: jvc_one(1,1)
    real(r8kind)                :: jvc_two(1,1)
    real(r8kind)                :: jvc_the(1,1)
    integer                     :: i, j, t, nitr, mxit
    real(r8kind)                :: jold, jnew, jthr, B, Q
    real(r8kind), allocatable   :: jtim(:) 
    integer                     :: nthreads, tid
    integer                     :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
    integer                     :: tim_len, bkg_len
    type(lorenz96_TL_type),  allocatable :: fwmod_vec(:)
    type(lorenz96_ADJ_type), allocatable :: bwmod_vec(:)
    type(lorenz96_type),     allocatable :: model(:)
    type(lorenz96_TL_type)               :: model_TL
    type(lorenz96_ADJ_type)              :: model_ADJ

    ! Allocate space for reused matrix products
    allocate (this%hrh_cov(bkg%ntimes, bkg%npoints, bkg%npoints))
    allocate (this%brh_cov(bkg%ntimes, bkg%npoints, bkg%npoints))
    allocate (this%anl_vec(bkg%ntimes, bkg%npoints))
    allocate (this%bht_ino(bkg%ntimes, bkg%npoints, 1))
    allocate (this%htr_ino(bkg%ntimes, bkg%npoints, 1))

    ! GET THE NEEDED REUSED MATRIX PRODUCTS:
    !        H(T)R(-1)(Y-HXb) = htr_ino
    !        H(T)R(-1)H       = hrh_cov
    ! B(-1/2)H(T)R(-1)(Y-HXb) = btr_ino
    ! B(-1/2)H(T)R(-1)H       = brh_cov
    ! INPUTS: Y, H, Xb, R(-1), B
    ! USE SAME VARIABLE NAME PRE AND POST
    call this%pre_sol(obs_opr, obs_cov, bkg_cov, inno_vec)

    ! Get the number of assimilation times
    tim_len = bkg%ntimes

    ! Get the size of the background
    bkg_len = bkg%npoints

    ! Allocate storage for intermediate results
    allocate (jtim(tim_len))
    allocate (tim_htr(bkg_len,       1))
    allocate (tim_bkc(bkg_len, bkg_len))
    allocate (tim_hrh(bkg_len, bkg_len))
    allocate (tim_bkv(bkg_len,       1))
    allocate (tim_anv(bkg_len,       1))

    allocate (new_vec(tim_len, bkg_len))
    allocate (tlm_vec(tim_len, bkg_len))
    allocate (mdl_vec(bkg_len,       1))
    allocate (ges_vec(bkg_len,       1))
    allocate (dif_vec(bkg_len,       1))
    allocate (dif_tra(      1, bkg_len))

    allocate (pre_dif(bkg_len,       1))
    allocate (pre_tra(      1, bkg_len))
    allocate (tmp_mat(      1, bkg_len))
    allocate (tmp_vec(bkg_len,       1))
    allocate (tmp_vvc(bkg_len,       1))
    allocate (grd_jvc(bkg_len,       1))


    ! Allocate/initialize fw and bw models if we are doing 4DVar
    if (this%method >= 3) then
      ! Allocate the arrays to hold fw and bw models
      allocate(model(tim_len))
      allocate (fwmod_vec(2:tim_len))
      allocate (bwmod_vec(1:tim_len - 1))

      ! Initialize a real model for use in calculating trajectories
      do t = 1, tim_len
        model(t) = lorenz96_type(t,'NETCDF')
      end do

      ! Load fw models
      do t = 2, tim_len
         fwmod_vec(t) = lorenz96_TL_type(t - 1, "NETCDF")
      end do

      ! Load bw models
      do t = 1, tim_len - 1
         bwmod_vec(t) = lorenz96_ADJ_type(t + 1, "NETCDF")
      end do
    end if

    !   PARAMETERS FOR VAR - SHOULD BE FROM NAMELIST
    nitr = 0
    mxit = 50
    jold = 100.0 
    jnew = 0.0
    jthr = 0.01 

    !   PARAMETERS FOR MODEL ERROR - ALSO SHOULD BE FROM NAMELIST
    !   B = RATIO OF B/B = 1.0
    !   Q = RATIO OF Q/B = 0.2 (OR, PERFECT TL/AD MODEL, Q=0.0)
    B    = 1.0
    Q    = 0.0
    print *, "SOLVER"

    ! IF 3DVAR - NO Q TERM
    if(this%method <= 2) Q = 0.0

    ! FIRST GUESS IS THE BACKGROUND
    this%anl_vec=bkg%state

    ! ITERATE TO SOLVE THE COST FUNCTION
    do while ( abs(jold - jnew) > jthr)
       if (nitr > mxit) exit
       if (jnew < 0.0) exit
       jold = jnew
       jnew = 0.0

       new_vec(:,:) = 0.0
       tlm_vec = this%anl_vec

!$OMP PARALLEL SHARED (bkg_cov, this, bkg, tlm_vec, jtim, new_vec, tim_len, B, Q, fwmod_vec, model) DEFAULT(PRIVATE)
!$OMP DO
       do t = 1, tim_len
          tid = OMP_GET_THREAD_NUM()
          tim_bkc(:,:) = bkg_cov%covariance(t,:,:)
          tim_hrh(:,:) = this%hrh_cov(t,:,:)
          tim_htr(:,:) = this%htr_ino(t,:,:)
          tim_bkv(:,1) = bkg%state(t,:)

          if(t == 1) then
             ! FOR THE FIRST TIME STEP, THERE IS NO PRIOR FIELD TO PROPAGATE
             mdl_vec(:,1) = tlm_vec(t,:)
             if (this%method == 4) new_vec(t,:) = mdl_vec(:,1)
          else
             ! RUN THE FORWARD MODEL FOR ALL STEPS AFTER THE FIRST
             mdl_vec(:,1) = tlm_vec(t-1,:)
             model_TL = fwmod_vec(t)
             model_TL%state(:) = mdl_vec(:,1)
             model(t)%state = model_TL%state
             call model(t)%adv_nsteps(1)
             model_TL%trajectory = model(t)%state - model_TL%state
             call model_TL%adv_nsteps(10)
             mdl_vec(:,1) = model_TL%state
             if (this%method == 4) new_vec(t,:) = mdl_vec(:,1)
             if (this%method /= 4) tlm_vec(t,:) = mdl_vec(:,1)
          end if

          ! CARRY ON WITH THE MINIMIZATION 
          tmp_vec(:,1) = (mdl_vec(:,1) - tim_bkv(:,1)) * B
          dif_vec(:,1) = (mdl_vec(:,1) - tim_bkv(:,1)) * B
          dif_tra = transpose(dif_vec)
          call this%pre_con_dif(tim_bkc, tmp_vec)
          pre_tra = transpose(tmp_vec)
          pre_dif(:,1) = tmp_vec(:,1)

          ! SOLVE FOR COST FUNCTION J
          ! FIRST TERM
          jvc_one = matmul(pre_tra, pre_dif)
          ! SECOND TERM
          tmp_mat = matmul(dif_tra, tim_hrh)
          jvc_two = matmul(tmp_mat, dif_vec)
          ! THIRD TERM
          jvc_the = matmul(dif_tra, tim_htr)
          ! COST FUNCTION
          jtim(t) = 0.5 * (jvc_one(1,1) + jvc_two(1,1) - 2.0 * jvc_the(1,1))

       end do
!$OMP END DO
!$OMP END PARALLEL

       if(this%method == 4) tlm_vec = new_vec
       do t = 1, tim_len
          jnew = jnew + jtim(t)
       end do
       jnew = jnew + this%jvc_for(1,1)
       new_vec(:,:) = 0.0

!$OMP PARALLEL SHARED (this, bkg_cov, bkg, tlm_vec, new_vec, tim_len, B, Q, bwmod_vec, model) DEFAULT(PRIVATE)
!$OMP DO
       !   CALCULATE GRAD-J IN REVERSE TEMPORAL ORDER 
       do t = tim_len, 1, -1
          tim_bkc(:,:) = bkg_cov%covariance(t,:,:)
          tim_hrh(:,:) = this%brh_cov(t,:,:)
          tim_htr(:,:) = this%bht_ino(t,:,:)
          tim_bkv(:,1) = bkg%state(t,:)

          if(t == tim_len) then
            ! FOR THE LAST (OR ONLY) TIME STEP - NO ADJOINT TO RUN
            mdl_vec(:,1) = tlm_vec(t,:)
            if(this%method == 4) mdl_vec(:,1) = 0.5 * (tlm_vec(t,:) + this%anl_vec(t,:))
          else
             ! THIS ONLY RUNS FOR 4DVAR
             ! FOR ALL OTHER TIME STEPS - ADJOINT NEEDED
             if (this%method == 3) mdl_vec(:,1) = tlm_vec(t+1,:)
             if (this%method == 4) mdl_vec(:,1) = this%anl_vec(t+1,:)
             model_ADJ = bwmod_vec(t)
             model_ADJ%state(:) = mdl_vec(:,1)
             model(t)%state = model_ADJ%state
             call model(t)%adv_nsteps(1)
             model_ADJ%trajectory(:) = -(model(t)%state - model_ADJ%state)
             call model_ADJ%adv_nsteps(10)
             mdl_vec(:,1) = model_ADJ%state
          end if

          ! CHOOSE THE FIRST GUESS FIELD
          if(this%method /= 4) ges_vec(:,1) = mdl_vec(:,1)
          if(this%method == 4) ges_vec(:,1) = 0.5 * (tlm_vec(t,:) + mdl_vec(:,1))

          ! CALCULATE THE GRADIENT OF THE COST FUNCTION
          ! FIRST - DIFFERENCE BETWEEN FIRST GUESS AND BACKGROUND
          tmp_vec(:,1) = (ges_vec(:,1) - tim_bkv(:,1)) * (B + Q)
          dif_vec(:,1) = (ges_vec(:,1) - tim_bkv(:,1)) * (B + Q)

          ! OBTAIN THE PRE-CONDITIONED DIFFERENCE BETWEEN THE BACKGROUND AND
          ! THE FIRST GUESS
          call this%pre_con_dif(tim_bkc, tmp_vec)
          pre_dif(:,1) = tmp_vec(:,1)
          tmp_vec = matmul(tim_hrh, dif_vec)

          tmp_vec(:,1) = pre_dif(:,1) + tmp_vec(:,1) - tim_htr(:,1)
          grd_jvc = matmul(tim_bkc, tmp_vec)
          new_vec(t,:) = ges_vec(:,1) - grd_jvc(:,1) * this%alpha

          if(this%method /= 4) tlm_vec(t,:) = new_vec(t,:)
       end do
!$OMP END DO
!$OMP END PARALLEL

       if(this%method == 4) tlm_vec = new_vec
       this%anl_vec = tlm_vec
       if (nitr > 0 .and. jnew > jold) jnew = jold
       if (nitr == 0) print *, 'initial cost = ', jnew
       nitr = nitr + 1 
       print *,"Cost at ",nitr,jnew
    end do

    print *,'final cost = ',jnew,' after ',nitr,' iterations'
    print *,"SOLVER COMPLETE"

  end subroutine solve


  ! BJE
  ! OUTPUT THE ANALYSIS VECTOR
  subroutine put_anl_vec(this)

    class(Solver_Type),     intent(in) :: this

    integer             :: i, t, ierr
    type(lorenz96_type) :: model

    print *,"PUT_ANL_VEC"

    ! Write new analysis to model output file
    model=lorenz96_type(2, "NETCDF")
    if (this%method <= 2) then
      model%state = this%anl_vec(1,:)
    else
      model%state = this%anl_vec(2,:)
    end if
    ierr = model%write_model_state("NETCDF")

40  FORMAT(A8, 2I5, F10.4)
    do t = 1, size(this%anl_vec, 1)
       do i = 1, size(this%anl_vec, 2)
          ! FOR 3DVAR
          if(this%method <=2) then
             write(*, 40) "FIN", 2, i, this%anl_vec(t,i)
          else
             write(*, 40) "FIN", t, i, this%anl_vec(t,i)
          end if
       end do
    end do

    print *,"PUT_ANL_VEC COMPLETE"

  end subroutine put_anl_vec


end module solver
