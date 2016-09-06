! BJE = 01 AUG 2016
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
  subroutine get_obs_vec(bkg_tim,bkg_pos,obs_tim,obs_pos,obs_vec)

    implicit none
    integer, intent(in)                      :: bkg_tim(:)
    integer, intent(in)                      :: bkg_pos(:,:)
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
  ! SPECIFICLY:  HRH_COV=(H(T)R(-1)H), HTR_INO=(H(T)R(-1))*(Y-HXb)
  !              BHT_INO=B(1/2)*(H(T)R(-1))*(Y-HXb)
  subroutine pre_sol(obs_opr,obs_cov,bkg_cov,hrh_cov,obs_vec,bht_ino,jvc_for)

    implicit none
    real(KIND=8), intent(in)    :: obs_cov(:,:,:) 
    real(KIND=8), intent(in)    :: bkg_cov(:,:,:) 
    real(KIND=8), intent(inout) :: hrh_cov(:,:,:) 
    real(KIND=8), intent(in)    :: obs_vec(:) 
    real(KIND=8), intent(in)    :: obs_opr(:,:,:) 
    real(KIND=8), intent(inout) :: bht_ino(:,:,:) 
    real(KIND=8), intent(inout) :: jvc_for(:,:) 
    integer                     :: i,j,t 

    real(KIND=8), allocatable   :: tim_obc(:,:)
    real(KIND=8), allocatable   :: tim_bkc(:,:)
    real(KIND=8), allocatable   :: tim_hrh(:,:)
    real(KIND=8), allocatable   :: tim_opr(:,:)
    real(KIND=8), allocatable   :: tim_bht(:,:)

    real(KIND=8), allocatable   :: obs_vvc(:,:)
    real(KIND=8), allocatable   :: tmp_jfo(:,:)
    real(KIND=8), allocatable   :: obs_opt(:,:)
    real(KIND=8), allocatable   :: tmp_htr(:,:)
    real(KIND=8), allocatable   :: tmp_bhh(:,:)
    real(KIND=8), allocatable   :: tmp_hhb(:,:)
    print *, "PRE_SOLVER"

    allocate (tim_obc(obs_len,obs_len))
    allocate (tim_bkc(bkg_len,bkg_len))
    allocate (tim_hrh(bkg_len,bkg_len))
    allocate (tim_opr(obs_len,bkg_len))
    allocate (tim_bht(bkg_len,      1))

    allocate (obs_vvc(obs_len,      1))
    allocate (tmp_jfo(obs_len,      1))
    allocate (obs_opt(bkg_len,obs_len))
    allocate (tmp_htr(bkg_len,obs_len))
    allocate (tmp_bhh(bkg_len,obs_len))
    allocate (tmp_hhb(obs_len,bkg_len))

    jvc_for(1,1)=0.0
    do t=1,tim_len
       tim_opr=obs_opr(t,:,:)
       tim_bkc=bkg_cov(t,:,:)
       tim_obc=obs_cov(t,:,:) 

       obs_opt=transpose(tim_opr)
       tmp_bhh=matmul(tim_bkc,obs_opt)
       tmp_hhb=matmul(tim_opr,tim_bkc)

       tmp_htr=matmul(tmp_bhh,tim_obc)
       tim_hrh=matmul(tmp_htr,tmp_hhb)

       obs_vvc(:,1)=obs_vec(:)
       tim_bht=matmul(tmp_htr,obs_vvc)

       bht_ino(t,:,:)=tim_bht(:,:) 
       hrh_cov(t,:,:)=tim_hrh(:,:)
       ! CREATE COST FUNCTION TERM (Y-HXb)R(-1)(Y-HXb), A CONSTANT
       tmp_jfo=matmul(tim_obc,obs_vvc)

       do i=1,obs_len
          jvc_for(1,1)=jvc_for(1,1)+obs_vvc(i,1)*tmp_jfo(i,1)
       end do
    end do

    print *,"PRE_SOLVER COMPLETE"

  end subroutine pre_sol


  ! BJE
  ! THE SOLVER FOR THE VARIATIONAL ASSIMILATION
  !  J    = (1/2)*(  V  )(T)*(    I     )*(  V  ) 
  !       + (1/2)*(X- Xb)(T)*(H(T)R(-1)H)*(X -Xb)
  !       + (1/2)*(Y-HXb)(T)*(    R(-1)H)*(X- Xb)
  !       + (1/2)*(X- Xb)(T)*(H(T)R(-1) )*(Y-HXb)
  !       + (1/2)*(Y-HXb)(T)*(    R(-1) )*(Y-HXb)
  !
  ! WHERE B(1/2)V= X-Xb
  !
  ! GRADJ =  V
  !       + (H(T)R(-1)H)*(X- Xb)
  !       - (H(T)R(-1) )*(y-HXb)
  !
  ! WHICH CAN BE REWRITTEN AS:
  !
  ! HRH=(B(1/2)H(T)R(-1)HB(1/2)), HTR=(B(1/2)H(T)R(-1))*(Y-HXb), VB(1/2)=(X- Xb) 
  !
  !  J    = (1/2)*V(T)*    V
  !       + (1/2)*V(T)*HRH*V
  !       +            HTR*V  
  !       + CONSTANT TERM 
  !
  ! GRADJ =            V
  !       +        HRH*V
  !       - B(1/2)*HTR
  subroutine var_solver(bkg_cov,hrh_cov,bht_ino,jvc_for,bkg_vec,anl_vec)

    implicit none
    real(KIND=8), intent(in)    :: bht_ino(:,:,:)
    real(KIND=8), intent(in)    :: bkg_cov(:,:,:)
    real(KIND=8), intent(in)    :: hrh_cov(:,:,:)
    real(KIND=8), intent(in)    :: bkg_vec(:,:)
    real(KIND=8), intent(inout) :: anl_vec(:,:)
    real(KIND=8), intent(in)    :: jvc_for(:,:)

    real(KIND=8), allocatable   :: tim_bht(:,:)
    real(KIND=8), allocatable   :: tim_bkc(:,:)
    real(KIND=8), allocatable   :: tim_hrh(:,:)
    real(KIND=8), allocatable   :: tim_anv(:,:)
    real(KIND=8), allocatable   :: tim_bkv(:,:)

    real(KIND=8), allocatable   :: pre_tra(:,:)
    real(KIND=8), allocatable   :: pre_bkg(:,:,:)
    real(KIND=8), allocatable   :: pre_anl(:,:,:)
    real(KIND=8), allocatable   :: pre_lam(:,:,:)
    real(KIND=8), allocatable   :: tmp_pra(:,:,:)
    real(KIND=8), allocatable   :: pre_dif(:,:)
    real(KIND=8), allocatable   :: tmp_mat(:,:)
    real(KIND=8), allocatable   :: tmp_vec(:,:)
    real(KIND=8), allocatable   :: grd_jvc(:,:)
    real(KIND=8)                :: jvc_one(1,1)
    real(KIND=8)                :: jvc_two(1,1)
    real(KIND=8)                :: jvc_the(1,1)
    real(KIND=8)                :: jvc_fiv(1,1)
    integer                     :: i,j,t,nitr,mxit
    real(KIND=8)                :: jold,jnew,jthr,alph
    real(KIND=8), allocatable   :: jtim(:) 
    integer                     :: nthreads, tid
    integer                     :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    allocate (jtim(tim_len)) 
    allocate (tim_bht(bkg_len,      1))
    allocate (tim_bkc(bkg_len,bkg_len))
    allocate (tim_hrh(bkg_len,bkg_len))
    allocate (tim_bkv(bkg_len,      1))
    allocate (tim_anv(bkg_len,      1))

    allocate (pre_bkg(tim_len,bkg_len,      1))
    allocate (pre_anl(tim_len,bkg_len,      1))
    allocate (pre_lam(tim_len,bkg_len,      1))
    allocate (pre_dif(        bkg_len,      1))
    allocate (tmp_pra(tim_len,bkg_len,      1))

    allocate (pre_tra(      1,bkg_len))
    allocate (tmp_mat(      1,bkg_len))
    allocate (tmp_vec(bkg_len,      1))
    allocate (grd_jvc(bkg_len,      1))

    nitr = 0
    mxit = 25 
    jold = 100.0 
    jnew = 0.0
    jthr = 0.01 
    alph = 0.0001
    print *,"SOLVER"

40  FORMAT(A8,I4,3F10.4)
50  FORMAT(2I4,3F10.4)

    do t=1,tim_len
       ! DO THE PRE-CONDITIONING OF X and Xb
       tim_bkc=bkg_cov(t,:,:)
       do i=1,bkg_len
          tim_bkv(i,1)=bkg_vec(t,i) 
          tim_anv(i,1)=bkg_vec(t,i)
       end do
       call pre_con_dif(tim_bkc,tim_bkv)
       call pre_con_dif(tim_bkc,tim_anv)
       do i=1,bkg_len
          pre_bkg(t,i,1)=tim_bkv(i,1)
          pre_anl(t,i,1)=tim_anv(i,1)
       end do
    end do

33  FORMAT(A10,3F10.5)
! ITERATE TO SOLVE THE COST FUNCTION
    do while ( abs(jold-jnew) > jthr)
       if (nitr.gt.mxit) exit
       if (jnew.lt.0.0) exit
       jold=jnew
       jnew=0.0

       write(*,33) "AFTER IT: ",pre_anl(:,3000,1)
       do t=2,tim_len
          tmp_pra(t,:,1)=pre_anl(t-1,:,1)
       end do
!$OMP PARALLEL SHARED (bht_ino,hrh_cov,pre_bkg,pre_anl,pre_lam,tmp_pra,bkg_cov,jtim,tim_len,bkg_len,nitr,mthd) DEFAULT(PRIVATE)
!$OMP DO
!   CALCULATE THE COST FUNCTION - J - FORWARD IN TIME
       do t=1,tim_len
          tid=OMP_GET_THREAD_NUM()
          tim_hrh(:,:)=hrh_cov(t,:,:)
          tim_bht(:,:)=bht_ino(t,:,:)
          tim_bkv(:,1)=pre_bkg(t,:,1)
          tim_anv(:,1)=pre_anl(t,:,1)

          !   RUN THE FORWARD MODEL FOR ALL STEPS AFTER THE FIRST
          if(t.gt.1) then 
             tim_bkc(:,:)=bkg_cov(t,:,:)
             tim_anv(:,1)=tmp_pra(t,:,1)
!            FOR NOT-TIME PARALLEL
             if (mthd.eq.3) print *,"NOT TIME PARALLEL",mthd 
             if (mthd.eq.3) tim_anv(:,1)=pre_anl(t-1,:,1)
             tmp_vec=matmul(tim_bkc,tim_anv)
             call forward_model(tmp_vec,t-1,1)
             call pre_con_dif(tim_bkc,tmp_vec)
             tim_anv=tmp_vec
          end if

          !   CARRY ON WITH THE MINIMIZATION 
          do i=1,bkg_len
             pre_dif(i,1)=tim_anv(i,1)-tim_bkv(i,1)
          end do
          pre_tra=transpose(pre_dif)   

          !   SOLVE FOR COST FUNCTION J
          !   FIRST TERM
          jvc_one=matmul(pre_tra,pre_dif)
          !   SECOND TERM
          tmp_mat=matmul(pre_tra,tim_hrh)
          jvc_two=matmul(tmp_mat,pre_dif)
          !   THIRD TERM
          jvc_the=matmul(pre_tra,tim_bht)

          !   FIFTH TERM (FOR TIME-PARALLEL)
          do i=1,bkg_len
             pre_dif(i,1)=tim_anv(i,1)-pre_anl(t,i,1)
             pre_lam(t,i,1)=pre_dif(i,1)
          end do
          pre_tra=transpose(pre_dif)
          jvc_fiv=matmul(pre_tra,pre_dif)

          !   COST FUNCTION
          jtim(t) = 0.5*(jvc_one(1,1)+jvc_two(1,1)-2.0*jvc_the(1,1)+jvc_fiv(1,1))

          tid = OMP_GET_THREAD_NUM()
44        FORMAT(A10,3I4,5F8.3)
          write(*,44) "JTIM =",nitr,t,tid,jtim(t),jvc_one(1,1),jvc_two(1,1),jvc_the(1,1),jvc_fiv(1,1)
          pre_anl(t,:,1)=tim_anv(:,1)
       end do
!$OMP END DO
!$OMP END PARALLEL
  
       write(*,33) "AFTER TL: ",pre_anl(:,3000,1)
       do t=1,tim_len
           jnew=jnew+jtim(t)
           if(t.lt.tim_len) tmp_pra(t,:,1)=pre_anl(t+1,:,1)
       end do
       jnew=jnew+jvc_for(1,1)

!$OMP PARALLEL SHARED (bht_ino,bkg_cov,hrh_cov,pre_bkg,pre_anl,pre_lam,tmp_pra,tim_len,bkg_len,alph,mthd) DEFAULT(PRIVATE)
!$OMP DO
       !   CALCULATE GRAD-J IN REVERSE TEMPORAL ORDER 
       do t=tim_len,1,-1
          tim_hrh(:,:)=hrh_cov(t,:,:)
          tim_bht(:,:)=bht_ino(t,:,:)
          tim_anv(:,1)=pre_anl(t,:,1)
          tim_bkv(:,1)=pre_bkg(t,:,1)
 
          if(t.lt.tim_len) then 
             tim_bkc=bkg_cov(t,:,:)
             tim_anv(:,1)=tmp_pra(t,:,1)
!            FOR NON TIME-PARALLEL 4DVAR
             if(mthd.eq.3) tim_anv(:,1)=pre_anl(t+1,:,1)
             tmp_vec=matmul(tim_bkc,tim_anv)
             call bakward_model(tmp_vec,t+1,1)
             call pre_con_dif(tim_bkc,tmp_vec)
             tim_anv(:,1)=tmp_vec(:,1)
          end if

          do i=1,bkg_len
             pre_dif(i,1)=tim_anv(i,1)-tim_bkv(i,1)
          end do

          tmp_vec=matmul(tim_hrh,pre_dif)

          do i=1,bkg_len
             grd_jvc(  i,1)=pre_dif(i,1)+tmp_vec(i,1)-tim_bht(i,1)
          end do

          do i=1,bkg_len
             pre_anl(t,i,1)=tim_anv(i,1)-grd_jvc(i,1)*alph
          end do

       end do
!$OMP END DO
!$OMP END PARALLEL

       write(*,33) "AFTER AD: ",pre_anl(:,3000,1)
       if (nitr.eq.0) print *,'initial cost = ',jnew
       nitr = nitr + 1 
       print *,"Cost at ",nitr,jnew
    end do

    print *,'final cost = ',jnew,' after ',nitr,' iterations'

    ! UNDO THE PRE-CONDITIONING X=B(1/2)V
    do t=1,tim_len
       tim_bkc(:,:)=bkg_cov(t,:,:)
       tim_anv(:,:)=pre_anl(t,:,:) 
       tmp_vec=matmul(tim_bkc,tim_anv)

       do i=1,bkg_len
          anl_vec(t,i)=tmp_vec(i,1)
       end do
    end do

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
          write(*,40) "FIN",t,i,anl_vec(t,i),bkg_vec(t,i),50.0+50.0*sin(((20.0*float(t-1)+float(i))/1000.0)*PI)
!         FOR 3DVAR
          if(mthd.le.2) write(*,40) "FIN",2,i,anl_vec(t,i),bkg_vec(t,i),50.0+50.0*sin(((20.0*float(2-1)+float(i))/1000.0)*PI)
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
