! BJE = 12 SEPT 2016
! PROGRAM TO DO VARIATIONAL ASSIMILATION

module module_varsolver

  use gptl
  use module_constants

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  implicit none

  ! Define namelists and default values
  integer          :: mthd  = 4
  integer          :: tim_len = 3
  integer          :: bkg_len = 4000
  integer          :: obs_len = 20
  namelist /params/  mthd, bkg_len, tim_len, obs_len

contains

  ! BJE
  ! GET THE METHOD TO USE FOR THE DATA ASSIMILATION
  ! 1 = 3DVAR, BKG AND OBS TIME MISMATCH
  ! 2 = 3DVAR, BKG AND OBS TIME MATCHED
  ! 3 = 4DVAR, NOT TIME-PARALLEL
  ! 4 = 4DVAR, TIME-PARALLEL
  subroutine get_method

    print *,"GET_METHOD"

    ! Read namelists from stdin
    read(stdin,nml=params)

    print *,"METHOD = ",mthd

    ! Force tim_len=1 for 3DVAR
    if(mthd.le.2) tim_len=1

    print *,"GET_METHOD COMPLETE"

  end subroutine get_method

  ! BJE
  ! GENERATE THE OBSERVATION ERROR COVARIANCE MATRIX, "R"
  subroutine get_obs_cov(obs_tim,obs_pos,obs_cov) 

    implicit none
    integer, intent(in)         :: obs_tim(:)
    integer, intent(in)         :: obs_pos(:)
    real(KIND=8), intent(inout) :: obs_cov(:,:,:)
    integer                     :: i,t 
    print *,"GET_OBS_COV_MAT"

    obs_cov(:,:,:)=0.0

  ! IN THIS CASE, R=I, THE IDENTITY MATRIX
    do i=1,obs_len
       obs_cov(obs_tim(i),i,i)=1.0
    end do

    print *,"GET_OBS_COV_MAT COMPLETE"

  end subroutine get_obs_cov

  ! BJE
  ! GENERATE THE BACKGROUND ERROR COVARIANCE MATRIX, "B"
  subroutine get_bkg_cov(bkg_cov)

    implicit none
    real(KIND=8), intent(inout) :: bkg_cov(:,:,:)
    real(KIND=8)                :: var
    integer                     :: t,i,j,jj,rad 

    print *,"GET_BKG_COV_MAT"

    var=3.61
    bkg_cov(:,:,:)=0.0
    rad=bkg_len/10

    do t=1,tim_len
       do i=1,bkg_len
          do j=-rad,+rad
             jj=i+j
             if(jj.gt.bkg_len) jj=jj-bkg_len
             if(jj.lt.1) jj=bkg_len+jj
             bkg_cov(t,i,jj)=var*exp(-((float(j)*0.005)**2))
          end do
       end do
    end do

    print *,"GET_BKG_COV_MAT COMPLETE"

  end subroutine get_bkg_cov


  ! BJE
  ! GENERATE PRE-CONDITIONED VECTOR V FROM X
  ! THIS REPLACED THE NEED TO EXPLICITLY INVERT B
  subroutine pre_con_dif(bkg_cov,non_vec)

    implicit none
    real(KIND=8), intent(inout) :: non_vec(:,:)
    real(KIND=8), intent(in)    :: bkg_cov(:,:)
    real(KIND=8), allocatable   :: con_vec(:,:)
    real(KIND=8), allocatable   :: bkg_cpy(:,:)
    integer                     :: i,j,jj,rad,info 
    integer, allocatable        :: ipiv(:)
    real(KIND=8)                :: var
    print *,"PRE_CON_DIF"

    allocate (ipiv(bkg_len))
    allocate (con_vec(bkg_len,1))
    allocate (bkg_cpy(bkg_len,bkg_len))
    con_vec=non_vec
    bkg_cpy=bkg_cov

    call dgesv(bkg_len,1,bkg_cpy,bkg_len,ipiv,con_vec,bkg_len,info)
    non_vec=con_vec

    print *,"PRE_CON_DIF COMPLETE"

  end subroutine pre_con_dif


  ! BJE
  ! GENERATE THE OBSERVATIONS "Y", AND THEIR LOCATIONS - FOR "H"
  subroutine get_obs_vec(obs_tim,obs_pos,obs_vec)

    implicit none
    integer, intent(inout), allocatable      :: obs_tim(:)
    integer, intent(inout), allocatable      :: obs_pos(:) 
    real(KIND=8), intent(inout), allocatable :: obs_vec(:) 
    integer                                  :: i,x 

    print *,"GET_OBS_VEC"

    allocate(obs_tim(obs_len))
    allocate(obs_pos(obs_len))
    allocate(obs_vec(obs_len))

    do i=1,obs_len
        x=modulo(i,8)
        obs_tim(i)=2
        obs_pos(i)=i*200 
        if(x.gt.3) then
            obs_tim(i)=1
            obs_pos(i)=i*200-35 
        end if
        if(x.gt.5) then
            obs_tim(i)=3
            obs_pos(i)=i*200-75 
        end if 
!       FOR MATCHED 3DVAR
        if(mthd.eq.2) obs_tim(i)=2
        obs_vec(i)=50.0+50.0*sin(((20.0*float(obs_tim(i)-1)+float(obs_pos(i)))/1000.0)*PI)
!       FOR 3DVAR
        if(mthd.le.2) obs_tim(i)=1
    end do

    print *,"GET_OBS_VEC COMPLETE" 

  end subroutine get_obs_vec


  ! BJE
  ! GET THE OBSERVATION OPERATOR, H, FROM THE INPUTS
  subroutine get_obs_opr(obs_tim,obs_pos,obs_opr)

    implicit none
    integer, intent(inout)      :: obs_tim(:)
    integer, intent(inout)      :: obs_pos(:)
    real(KIND=8), intent(inout) :: obs_opr(:,:,:)
    integer                     :: i,t 
    print *,"GET_OBS_VEC"

    obs_opr(:,:,:)=0.0

    do i=1,obs_len
       obs_opr(obs_tim(i),i,obs_pos(i))=1.0
    end do

    print *,"GET_OBS_OPR COMPLETE"

  end subroutine get_obs_opr

  ! BJE
  ! GENERATE THE FIRST GUESS "Xb" - SHOULD USE A REAL MODEL
  subroutine get_bkg_vec(bkg_tim,bkg_pos,bkg_vec)

    implicit none
    real(KIND=8), intent(inout), allocatable :: bkg_vec(:,:)
    integer,intent(inout), allocatable       :: bkg_pos(:,:)
    integer,intent(inout), allocatable       :: bkg_tim(:) 
    integer                                  :: i,t,tt
    print *,"GET_BKG_VEC"

    allocate (bkg_vec(tim_len,bkg_len))
    allocate (bkg_pos(tim_len,bkg_len))
    allocate (bkg_tim(tim_len))
    do t=1,tim_len
       bkg_tim(t)=t
       tt=t
       if(mthd.le.2) tt=2
       do i=1,bkg_len
          bkg_pos(t,i)=i
          bkg_vec(t,i)=50.0+50.0*sin(((20.0*float(tt-2)+float(i))/1000.0)*PI)
       end do
    end do

    print *,"GET_BKG_VEC COMPLETE"

  end subroutine get_bkg_vec

  ! BJE
  ! GENERATE THE INNOVATION VECTOR (Y-HXb)
  ! USE OBS_VEC TO STORE THE OUTPUT
  subroutine get_ino_vec(bkg_tim,bkg_pos,obs_tim,obs_pos,obs_vec,bkg_vec)

    implicit none
    integer, intent(in)         :: bkg_tim(:)
    integer, intent(in)         :: bkg_pos(:,:)
    integer, intent(in)         :: obs_tim(:)
    integer, intent(in)         :: obs_pos(:)
    real(KIND=8), intent(inout) :: obs_vec(:)
    real(KIND=8), intent(in)    :: bkg_vec(:,:)
    integer                     :: i 
    print *,"GET_INO_VEC"

40  FORMAT(A8,3I5,2F10.4)

    do i=1,obs_len
       obs_vec(i)=obs_vec(i)-bkg_vec(obs_tim(i),obs_pos(i))
       write(*,40) "INO ",i,obs_tim(i),obs_pos(i),obs_vec(i),bkg_vec(obs_tim(i),obs_pos(i))
    end do

    print *,"GET_INO_VEC COMPLETE"

  end subroutine get_ino_vec


  ! BJE
  ! PREPARES MATRICES FOR USE IN THE SOLVER
  ! THERE ARE A NUMBER OF MATRICES THAT ARE REUSED IN THE SOLVER
  ! SPECIFICLY:  HRH_COV=       (H(T)R(-1)H)
  !              HTR_INO=       (H(T)R(-1))*(Y-HXb)
  !              BRH_COV=B(1/2)*(H(T)R(-1)H)
  !              BHT_INO=B(1/2)*(H(T)R(-1))*(Y-HXb)

  subroutine pre_sol(obs_opr,obs_cov,bkg_cov,hrh_cov,brh_cov,obs_vec,htr_ino,bht_ino,jvc_for)

    implicit none
    real(KIND=8), intent(in)    :: obs_cov(:,:,:) 
    real(KIND=8), intent(in)    :: bkg_cov(:,:,:)
    real(KIND=8), intent(inout) :: hrh_cov(:,:,:) 
    real(KIND=8), intent(inout) :: brh_cov(:,:,:)
    real(KIND=8), intent(in)    :: obs_vec(:) 
    real(KIND=8), intent(in)    :: obs_opr(:,:,:) 
    real(KIND=8), intent(inout) :: htr_ino(:,:,:) 
    real(KIND=8), intent(inout) :: bht_ino(:,:,:)
    real(KIND=8), intent(inout) :: jvc_for(:,:) 
    integer                     :: i,j,t 

    real(KIND=8), allocatable   :: tmp_mat(:,:)
    real(KIND=8), allocatable   :: tmp_vec(:,:)

    real(KIND=8), allocatable   :: tim_bkc(:,:)
    real(KIND=8), allocatable   :: tim_obc(:,:)
    real(KIND=8), allocatable   :: tim_hrh(:,:)
    real(KIND=8), allocatable   :: tim_opr(:,:)
    real(KIND=8), allocatable   :: tim_htr(:,:)

    real(KIND=8), allocatable   :: obs_vvc(:,:)
    real(KIND=8), allocatable   :: tmp_jfo(:,:)
    real(KIND=8), allocatable   :: obs_opt(:,:)
    real(KIND=8), allocatable   :: tmp_rhh(:,:)
    real(KIND=8), allocatable   :: tmp_hrr(:,:)
    print *, "PRE_SOLVER"

    allocate (tmp_mat(bkg_len,bkg_len))
    allocate (tmp_vec(bkg_len,      1))

    allocate (tim_bkc(bkg_len,bkg_len))
    allocate (tim_obc(obs_len,obs_len))
    allocate (tim_hrh(bkg_len,bkg_len))
    allocate (tim_opr(obs_len,bkg_len))
    allocate (tim_htr(bkg_len,      1))

    allocate (obs_vvc(obs_len,      1))
    allocate (tmp_jfo(obs_len,      1))
    allocate (obs_opt(bkg_len,obs_len))
    allocate (tmp_rhh(obs_len,bkg_len))
    allocate (tmp_hrr(bkg_len,obs_len))

    jvc_for(1,1)=0.0
    do t=1,tim_len
       ! ASSUME THAT OBS_OPR=H, OBS_COV=R(-1/2), BKG_COV=B(1/2), OBS_VEC=(Y-HXb)
       tim_opr(:,:)=obs_opr(t,:,:)
       tim_obc(:,:)=obs_cov(t,:,:) 
       tim_bkc(:,:)=bkg_cov(t,:,:)
       obs_vvc(:,1)=obs_vec(:)

       ! CREATE THE OBS BASED MATRICES, FOR USE IN CALCULATING THE COST FUNCTION
       ! tim_hrh=H(T)R(-1)H, tim_htr=R(-1)H
       tmp_rhh=matmul(tim_obc,tim_opr)       !     R(-1/2)H 
       obs_opt=transpose(tmp_rhh)            ! H(T)R(-1/2)
       tim_hrh=matmul(obs_opt,tmp_rhh)       ! H(T)R(-1/2)H 
       hrh_cov(t,:,:)=tim_hrh(:,:)

       tmp_hrr=matmul(obs_opt,tim_obc)       ! H(T)R(-1) 
       tim_htr=matmul(tmp_hrr,obs_vvc)       ! H(T)R(-1)(Y-HXb)
       htr_ino(t,:,1)=tim_htr(:,1) 
 
       ! CREATE COST FUNCTION TERM (Y-HXb)R(-1)(Y-HXb), A CONSTANT
       tmp_jfo=matmul(tim_obc,obs_vvc)       ! R(-1/2)(Y-HXb)
       do i=1,obs_len
          jvc_for(1,1)=jvc_for(1,1)+tmp_jfo(i,1)*tmp_jfo(i,1)
       end do

       ! CREATE THE UNPRECONDITIONED MATRICES, FOR THE GRADIENT OF J
       tmp_vec=matmul(tim_bkc,tim_htr)
       tmp_mat=matmul(tim_bkc,tim_hrh)
       bht_ino(t,:,1)=tmp_vec(:,1)
       brh_cov(t,:,:)=tmp_mat(:,:)
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
  subroutine var_solver(bkg_cov,hrh_cov,brh_cov,htr_ino,bht_ino,jvc_for,bkg_vec,anl_vec)

    implicit none
    real(KIND=8), intent(in)    :: htr_ino(:,:,:)
    real(KIND=8), intent(in)    :: bht_ino(:,:,:)
    real(KIND=8), intent(in)    :: bkg_cov(:,:,:)
    real(KIND=8), intent(in)    :: hrh_cov(:,:,:)
    real(KIND=8), intent(in)    :: brh_cov(:,:,:)
    real(KIND=8), intent(in)    :: bkg_vec(:,:)
    real(KIND=8), intent(inout) :: anl_vec(:,:)
    real(KIND=8), intent(in)    :: jvc_for(:,:)

    real(KIND=8), allocatable   :: tim_htr(:,:)
    real(KIND=8), allocatable   :: tim_bkc(:,:)
    real(KIND=8), allocatable   :: tim_hrh(:,:)
    real(KIND=8), allocatable   :: tim_anv(:,:)
    real(KIND=8), allocatable   :: tim_bkv(:,:)

    real(KIND=8), allocatable   :: new_vec(:,:)
    real(KIND=8), allocatable   :: tlm_vec(:,:)
    real(KIND=8), allocatable   :: mdl_vec(:,:)
    real(KIND=8), allocatable   :: ges_vec(:,:)
    real(KIND=8), allocatable   :: dif_vec(:,:)
    real(KIND=8), allocatable   :: dif_tra(:,:)

    real(KIND=8), allocatable   :: pre_tra(:,:)
    real(KIND=8), allocatable   :: pre_dif(:,:)
    real(KIND=8), allocatable   :: tmp_mat(:,:)
    real(KIND=8), allocatable   :: tmp_vec(:,:)
    real(KIND=8), allocatable   :: tmp_vvc(:,:)

    real(KIND=8), allocatable   :: grd_jvc(:,:)
    real(KIND=8)                :: jvc_one(1,1)
    real(KIND=8)                :: jvc_two(1,1)
    real(KIND=8)                :: jvc_the(1,1)
    integer                     :: i,j,t,nitr,mxit
    real(KIND=8)                :: jold,jnew,jthr,alph,B,Q 
    real(KIND=8), allocatable   :: jtim(:) 
    integer                     :: nthreads, tid
    integer                     :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    allocate (jtim(tim_len)) 
    allocate (tim_htr(bkg_len,      1)) 
    allocate (tim_bkc(bkg_len,bkg_len))
    allocate (tim_hrh(bkg_len,bkg_len))
    allocate (tim_bkv(bkg_len,      1))
    allocate (tim_anv(bkg_len,      1))

    allocate (new_vec(tim_len,bkg_len))
    allocate (tlm_vec(tim_len,bkg_len))
    allocate (mdl_vec(bkg_len,      1))
    allocate (ges_vec(bkg_len,      1))
    allocate (dif_vec(bkg_len,      1))
    allocate (dif_tra(      1,bkg_len))

    allocate (pre_dif(bkg_len,      1))
    allocate (pre_tra(      1,bkg_len))
    allocate (tmp_mat(      1,bkg_len))
    allocate (tmp_vec(bkg_len,      1))
    allocate (tmp_vvc(bkg_len,      1))
    allocate (grd_jvc(bkg_len,      1))

!   PARAMETERS FOR VAR - SHOULD BE FROM NAMELIST
    nitr = 0
    mxit = 25 
    jold = 100.0 
    jnew = 0.0
    jthr = 0.01 
    alph = 0.00005
    if(mthd.eq.4) alph=alph*3.8
    if(mthd.eq.1) alph=alph*4.8

!   PARAMETERS FOR MODEL ERROR - ALSO SHOULD BE FROM NAMELIST
!   B = RATIO OF B/B = 1.0  
!   Q = RATIO OF Q/B = 0.2 (OR, PERFECT TL/AD MODEL, Q=0.0) 
    B    = 1.0 
    Q    = 0.0
    print *,"SOLVER"

! IF 3DVAR - NO Q TERM
    if(mthd.le.2) Q = 0.0

! FIRST GUESS IS THE BACKGROUND
    anl_vec=bkg_vec

! ITERATE TO SOLVE THE COST FUNCTION
    do while ( abs(jold-jnew) > jthr)
       if (nitr.gt.mxit) exit
       if (jnew.lt.0.0) exit
       jold=jnew
       jnew=0.0
 
       new_vec(:,:)=0.0
       tlm_vec=anl_vec
!$OMP PARALLEL SHARED (bkg_cov,htr_ino,hrh_cov,bkg_vec,tlm_vec,jtim,new_vec,tim_len,mthd,B,Q) DEFAULT(PRIVATE)
!$OMP DO
       do t=1,tim_len
          tid=OMP_GET_THREAD_NUM()
          tim_bkc(:,:)=bkg_cov(t,:,:)
          tim_hrh(:,:)=hrh_cov(t,:,:)
          tim_htr(:,:)=htr_ino(t,:,:)
	  tim_bkv(:,1)=bkg_vec(t,:) 

          if(t.eq.1) then
!            FOR THE FIRST TIME STEP, THERE IS NO PRIOR FIELD TO PROPAGATE
	     mdl_vec(:,1)=tlm_vec(t,:)
             if (mthd.eq.4) new_vec(t,:)=mdl_vec(:,1)
          else
!            RUN THE FORWARD MODEL FOR ALL STEPS AFTER THE FIRST
	     mdl_vec(:,1)=tlm_vec(t-1,:)
             call forward_model(mdl_vec,t-1,1)
             if (mthd.eq.4) new_vec(t,:)=mdl_vec(:,1)
             if (mthd.ne.4) tlm_vec(t,:)=mdl_vec(:,1)
          end if

          !   CARRY ON WITH THE MINIMIZATION 
          tmp_vec(:,1)=(mdl_vec(:,1)-tim_bkv(:,1))*B 
          dif_vec(:,1)=(mdl_vec(:,1)-tim_bkv(:,1))*B
	  dif_tra=transpose(dif_vec)
          call pre_con_dif(tim_bkc,tmp_vec)
          pre_tra=transpose(tmp_vec)
          pre_dif(:,1)=tmp_vec(:,1)
 
          !   SOLVE FOR COST FUNCTION J
          !   FIRST TERM
          jvc_one=matmul(pre_tra,pre_dif)
          !   SECOND TERM
          tmp_mat=matmul(dif_tra,tim_hrh)
          jvc_two=matmul(tmp_mat,dif_vec)
          !   THIRD TERM
          jvc_the=matmul(dif_tra,tim_htr)

          !   COST FUNCTION
          jtim(t) = 0.5*(jvc_one(1,1)+jvc_two(1,1)-2.0*jvc_the(1,1)) 

       end do
!$OMP END DO
!$OMP END PARALLEL

       if(mthd.eq.4) tlm_vec=new_vec 
       do t=1,tim_len
          jnew=jnew+jtim(t)
       end do
       jnew=jnew+jvc_for(1,1)
       new_vec(:,:)=0.0

!$OMP PARALLEL SHARED (bht_ino,bkg_cov,brh_cov,bkg_vec,anl_vec,tlm_vec,new_vec,tim_len,alph,mthd,B,Q) DEFAULT(PRIVATE)
!$OMP DO
       !   CALCULATE GRAD-J IN REVERSE TEMPORAL ORDER 
       do t=tim_len,1,-1
          tim_bkc(:,:)=bkg_cov(t,:,:)
          tim_hrh(:,:)=brh_cov(t,:,:)
          tim_htr(:,:)=bht_ino(t,:,:)
          tim_bkv(:,1)=bkg_vec(t,:)
         
          if(t.eq.tim_len) then 
!            FOR THE LAST (OR ONLY) TIME STEP - NO ADJOINT TO RUN
	     mdl_vec(:,1)=tlm_vec(t,:)
             if(mthd.eq.4) mdl_vec(:,1)=0.5*(tlm_vec(t,:)+anl_vec(t,:))
          else 			
!            THIS ONLY RUNS FOR 4DVAR
!	     FOR ALL OTHER TIME STEPS - ADJOINT NEEDED
             if (mthd.eq.3) mdl_vec(:,1)=tlm_vec(t+1,:)
             if (mthd.eq.4) mdl_vec(:,1)=anl_vec(t+1,:)
             call bakward_model(mdl_vec,t+1,1)
          end if

!         CHOOSE THE FIRST GUESS FIELD
          if(mthd.ne.4) ges_vec(:,1)=mdl_vec(:,1)
          if(mthd.eq.4) ges_vec(:,1)=0.5*(tlm_vec(t,:)+mdl_vec(:,1))

!         CALCULATE THE GRADIENT OF THE COST FUNCTION
!	  FIRST - DIFFERENCE BETWEEN FIRST GUESS AND BACKGROUND
	  tmp_vec(:,1)=(ges_vec(:,1)-tim_bkv(:,1))*(B+Q)
	  dif_vec(:,1)=(ges_vec(:,1)-tim_bkv(:,1))*(B+Q)

!         OBTAIN THE PRE-CONDITIONED DIFFERENCE BETWEEN THE BACKGROUND AND 
!         THE FIRST GUESS
          call pre_con_dif(tim_bkc,tmp_vec)
          pre_dif(:,1)=tmp_vec(:,1)
          tmp_vec=matmul(tim_hrh,dif_vec)

          tmp_vec(:,1)=pre_dif(:,1)+tmp_vec(:,1)-tim_htr(:,1)
	  grd_jvc=matmul(tim_bkc,tmp_vec)
          new_vec(t,:)=ges_vec(:,1)-grd_jvc(:,1)*alph

          if(mthd.ne.4) tlm_vec(t,:)=new_vec(t,:)
       end do
!$OMP END DO
!$OMP END PARALLEL
       
       if(mthd.eq.4) tlm_vec=new_vec
       anl_vec=tlm_vec
       if (nitr.gt.0 .and. jnew.gt.jold) jnew=jold
       if (nitr.eq.0) print *,'initial cost = ',jnew
       nitr = nitr + 1 
       print *,"Cost at ",nitr,jnew
    end do

    print *,'final cost = ',jnew,' after ',nitr,' iterations'
    print *,"SOLVER COMPLETE"

  end subroutine var_solver


  ! BJE
  ! OUTPUT THE ANALYSIS VECTOR
  subroutine put_anl_vec(anl_vec,bkg_vec,bkg_tim)

    real(KIND=8), intent(in)    :: anl_vec(:,:) 
    real(KIND=8), intent(in)    :: bkg_vec(:,:) 
    integer, intent(in)         :: bkg_tim(:)
    integer                     :: i,t
    print *,"PUT_ANL_VEC"


40  FORMAT(A8,2I5,3F10.4)
    do t=1,tim_len
       do i=1,bkg_len
!         FOR 3DVAR
          if(mthd.le.2) then
	     write(*,40) "FIN",2,i,anl_vec(t,i),bkg_vec(t,i),50.0+50.0*sin(((20.0*float(2-1)+float(i))/1000.0)*PI)
	  else
             write(*,40) "FIN",t,i,anl_vec(t,i),bkg_vec(t,i),50.0+50.0*sin(((20.0*float(t-1)+float(i))/1000.0)*PI)
	  end if
       end do
    end do

    print *,"PUT_ANL_VEC COMPLETE"

  end subroutine put_anl_vec


  ! THE FORWARD MODEL FOR MY SINE WAVE
  ! MOVE IT 20 POINTS EVERY TIME STEP
  subroutine forward_model(mod_vec,t,steps)

    real(KIND=8), intent(inout)     :: mod_vec(:,:)
    integer, intent(in)             :: t,steps
    integer                         :: i
    print *,"FORWARD MODEL"

35  FORMAT (A4,3I4,2F10.5)
    do i=1,bkg_len
       mod_vec(i,1)=mod_vec(i,1)+PI*cos(((20.0*(float(t)-1.5)+float(i))/1000.0)*PI)*float(steps)
    end do

    print *,"END FORWARD_MODEL"

  end subroutine forward_model


  ! THE BACKWARD MODEL FOR MY SINE WAVE
  ! MOVE IT 20 POINTS EVERY TIME STEP
  subroutine bakward_model(mod_vec,t,steps)

    real(KIND=8), intent(inout)     :: mod_vec(:,:)
    integer, intent(in)             :: t,steps
    integer                         :: i

    print *,"BAKWARD_MODEL"

35  FORMAT (A4,3I4,2F10.5)
    do i=1,bkg_len
       mod_vec(i,1)=mod_vec(i,1)-PI*cos(((20.0*(float(t)-2.5)+float(i))/1000.0)*PI)*float(steps)
    end do

    print *,"END BAKWARD_MODEL"

  end subroutine bakward_model


  ! THAT IS ALL FOLKS!
  !
end module module_varsolver
