! BJE = 01 AUG 2016
! PROGRAM TO DO VARIATIONAL ASSIMILATION

program adept
use gptl

implicit none
integer                      ::      obs_len, bkg_len
integer, allocatable         ::      obs_pos(:) 
real(KIND=8)                 ::      jvc_for(1,1)
real(KIND=8), allocatable    ::      obs_opr(:,:)
real(KIND=8), allocatable    ::      obs_cov(:,:)
real(KIND=8), allocatable    ::      bkg_cov(:,:)
real(KIND=8), allocatable    ::      hrh_cov(:,:)
real(KIND=8), allocatable    ::      obs_vec(:)
real(KIND=8), allocatable    ::      bkg_vec(:)
real(KIND=8), allocatable    ::      anl_vec(:)
real(KIND=8), allocatable    ::      bht_ino(:,:)
real                         ::      ret 

! BJE
! INITIALIZE GPTL AND START A TIMER
ret = gptlinitialize ()                    
ret = gptlstart ('adept')                 

! BJE
! FIRST - NEED TO KNOW HOW MANY OBSERVATIONS AND STATE VECTOR
! OBTAIN THE OBSERATIONS, Y, AND THE BACKGROUND, Xb 
call get_bkg_vec(bkg_len,bkg_vec)
call get_obs_vec(bkg_len,obs_len,obs_pos,obs_vec)

! BJE
! KNOWING THE NUMBERS, ALLOCATE VECTORS/MATRICIES (ARRAYS) ACCORTINGLY
allocate (obs_opr(obs_len,bkg_len))
allocate (obs_cov(obs_len,obs_len))
allocate (bkg_cov(bkg_len,bkg_len))
allocate (hrh_cov(bkg_len,bkg_len))
allocate (anl_vec(bkg_len))
allocate (bht_ino(bkg_len,1))

! BJE
! GET THE INNOVATION VECTOR - (Y-HXb) - OVERWRITE OBS_VEC
call get_ino_vec(obs_pos,obs_vec,bkg_vec)

! BJE
! KNOWING THE LOCATION OF THE OBS, CREATE OBS OPERATOR, H 
call get_obs_opr(obs_pos,obs_opr)

! BJE   
! OBTAIN THE COVARIANCE MATRIX R - OBS ERROR
! SOMEDAY WILL COME TO US FROM THE UFO
call get_obs_cov(obs_pos,obs_cov)

! BJE
! OBTAIN THE B(1/2) MATRIX - FOR PRE CONDITIONING
call get_bkg_cov(bkg_cov)

! BJE
! GET THE NEEDED REUSED MATRIX PRODUCTS:
! B(1/2)(Y-HXb)(T)R(-1)H        = bht_ino
! B(1/2)      H(T)R(-1)HB(1/2)  = hrh_cov
! INPUTS: Y, H, Xb, R(-1)
! USE SAME VARIABLE NAME PRE AND POST
call pre_sol(obs_opr,obs_cov,bkg_cov,hrh_cov,obs_vec,bht_ino,jvc_for)

! BJE
! THE MAIN EVENT - THE SOLVER
call var_solver(bkg_cov,hrh_cov,bht_ino,jvc_for,bkg_vec,anl_vec)

! BJE
! OUTPUT THE NEW ANALYSIS
call put_anl_vec(anl_vec,bkg_vec)

! BJE
! END THE TIMER AND OUTPUT THE GPTL RESULTS
ret = gptlstop ('adept') 
ret = gptlpr (0) 
!ret = gptlpr_summary (MPI_COMM_WORLD) 
ret = gptlfinalize ()

! THAT IS THE END OF THE MAIN PROGRAM
! NO, REALLY, THAT WAS THE END OF THE MAIN PROGRAM
! SERIOUSLY< THERE IS NOTHING MORE IN THE MAIN PROGRAM
! ON TO THE SUBROUTINES!
!
!--------------------------------------------------------------------
contains


! BJE
! GENERATE THE OBSERVATION ERROR COVARIANCE MATRIX, "R"
subroutine get_obs_cov(obs_pos,obs_cov) 
implicit none
integer, intent(in)         :: obs_pos(:)
real(KIND=8), intent(inout) :: obs_cov(:,:)
integer                     :: obs_len
integer                     :: i
print *,"GET_OBS_COV_MAT"

obs_len=size(obs_cov,1)

do i=1,obs_len
    obs_cov(i,i)=1.0
end do

print *,"GET_OBS_COV_MAT COMPLETE"
end subroutine get_obs_cov

! BJE
! GENERATE THE BACKGROUND ERROR COVARIANCE MATRIX, "B"
subroutine get_bkg_cov(bkg_cov)
implicit none
real(KIND=8), intent(inout) :: bkg_cov(:,:)
real(KIND=8)                :: var
integer                     :: i,j,jj,rad 
integer                     :: bkg_len
print *,"GET_BKG_COV_MAT"

bkg_len=size(bkg_cov,1)
var=3.61
bkg_cov(:,:)=0.0
rad=bkg_len/10

do i=1,bkg_len
    do j=-rad,+rad
        jj=i+j
        if(jj.gt.bkg_len) jj=jj-bkg_len
        if(jj.lt.1) jj=bkg_len+jj
        bkg_cov(i,jj)=var*exp(-((float(j)*0.5)**2))
        if(i.eq.rad) print *,i,j,jj,bkg_cov(i,jj)
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
integer                     :: bkg_len
integer, allocatable        :: ipiv(:)
real(KIND=8)                :: var
print *,"PRE_CON_DIF"

bkg_len=size(bkg_cov,1)
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
subroutine get_obs_vec(bkg_len,obs_len,obs_pos,obs_vec)
implicit none
integer, intent(in)                      :: bkg_len
integer, intent(inout)                   :: obs_len
integer, intent(inout), allocatable      :: obs_pos(:) 
real(KIND=8), intent(inout), allocatable :: obs_vec(:) 
integer                                  :: i
print *,"GET_OBS_VEC"

obs_len=8
allocate(obs_pos(obs_len))
allocate(obs_vec(obs_len))

do i=1,obs_len
    if((2**(i-1)).lt.bkg_len) then 
        obs_pos(i)=2**(i-1)
    else
        obs_pos(i)=bkg_len-(i*2)
    end if
end do

do i=1,obs_len
  obs_vec(i)=50.0+50.0*sin(((0.2+float(obs_pos(i)))/10.0)*(3.1416))
end do

print *,"GET_OBS_VEC COMPLETE" 
end subroutine get_obs_vec


! BJE
! GET THE OBSERVATION OPERATOR, H, FROM THE INPUTS
subroutine get_obs_opr(obs_pos,obs_opr)
implicit none
integer, intent(inout)      :: obs_pos(:)
real(KIND=8), intent(inout) :: obs_opr(:,:)
integer                     :: bkg_len
integer                     :: obs_len
integer                     :: i
print *,"GET_OBS_VEC"

obs_len=size(obs_opr,1)
bkg_len=size(obs_opr,2)
obs_opr(:,:)=0.0

do i=1,obs_len
    obs_opr(i,obs_pos(i))=1.0
end do

print *,"GET_OBS_OPR COMPLETE"
end subroutine get_obs_opr

! BJE
! GENERATE THE FIRST GUESS "Xb" - SHOULD USE A REAL MODEL
subroutine get_bkg_vec(bkg_len,bkg_vec)
implicit none
real(KIND=8), intent(inout),allocatable :: bkg_vec(:)
integer,intent(inout)                   :: bkg_len
integer                                 :: i
print *,"GET_BKG_VEC"

bkg_len=40
allocate (bkg_vec(bkg_len))
do i=1,bkg_len
    bkg_vec(i)=50.0+50.0*sin((float(i)/10.0)*(3.1416))
end do

print *,"GET_BKG_VEC COMPLETE"
end subroutine get_bkg_vec

! BJE
! GENERATE THE INNOVATION VECTOR (Y-HXb)
! USE OBS_VEC TO STORE THE OUTPUT
subroutine get_ino_vec(obs_pos,obs_vec,bkg_vec)
implicit none
integer, intent(in)         :: obs_pos(:)
real(KIND=8), intent(inout) :: obs_vec(:)
real(KIND=8), intent(in)    :: bkg_vec(:)
integer                     :: obs_len, bkg_len
integer                     :: i
print *,"GET_INO_VEC"
obs_len=size(obs_vec)
bkg_len=size(bkg_vec)

40 FORMAT(A8,2I4,2F10.4)

do i=1,obs_len
    obs_vec(i)=obs_vec(i)-bkg_vec(obs_pos(i))
    write(*,40) "INO ",i,obs_pos(i),obs_vec(i),bkg_vec(obs_pos(i))
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
real(KIND=8), dimension(:,:), intent(in)    :: obs_cov
real(KIND=8), dimension(:,:), intent(in)    :: bkg_cov
real(KIND=8), dimension(:,:), intent(inout) :: hrh_cov
real(KIND=8), dimension(:),   intent(in)    :: obs_vec
real(KIND=8), dimension(:,:), intent(in)    :: obs_opr
real(KIND=8), dimension(:,:), intent(inout) :: bht_ino
real(KIND=8), dimension(:,:), intent(inout) :: jvc_for
integer                                     :: i,j
integer                                     :: obs_len, bkg_len
real(KIND=8), allocatable                   :: obs_vvc(:,:)
real(KIND=8), allocatable                   :: tmp_jfo(:,:)
real(KIND=8), allocatable                   :: obs_opt(:,:)
real(KIND=8), allocatable                   :: tmp_htr(:,:)
real(KIND=8), allocatable                   :: tmp_bhh(:,:)
real(KIND=8), allocatable                   :: tmp_hhb(:,:)
print *, "PRE_SOLVER"

obs_len=size(obs_opr,1)
bkg_len=size(obs_opr,2)
allocate (obs_vvc(obs_len,      1))
allocate (tmp_jfo(obs_len,      1))
allocate (obs_opt(bkg_len,obs_len))
allocate (tmp_htr(bkg_len,obs_len))
allocate (tmp_bhh(bkg_len,obs_len))
allocate (tmp_hhb(obs_len,bkg_len))

obs_opt=transpose(obs_opr)
tmp_bhh=matmul(bkg_cov,obs_opt)
tmp_hhb=matmul(obs_opr,bkg_cov)

tmp_htr=matmul(tmp_bhh,obs_cov)
hrh_cov=matmul(tmp_htr,tmp_hhb)

obs_vvc(:,1)=obs_vec(:)
bht_ino=matmul(tmp_htr,obs_vvc)

! CREATE COST FUNCTION TERM (Y-HXb)R(-1)(Y-HXb), A CONSTANT
tmp_jfo=matmul(obs_cov,obs_vvc)
jvc_for(1,1)=0.0
do i=1,obs_len
    jvc_for(1,1)=jvc_for(1,1)+obs_vvc(i,1)*tmp_jfo(i,1)
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
real(KIND=8), intent(in)    :: bht_ino(:,:)
real(KIND=8), intent(in)    :: bkg_cov(:,:)
real(KIND=8), intent(in)    :: hrh_cov(:,:)
real(KIND=8), intent(in)    :: bkg_vec(:)
real(KIND=8), intent(inout) :: anl_vec(:)
real(KIND=8), intent(in)    :: jvc_for(:,:)

real(KIND=8), allocatable   :: pre_tra(:,:)
real(KIND=8), allocatable   :: pre_bkg(:,:)
real(KIND=8), allocatable   :: pre_anl(:,:)
real(KIND=8), allocatable   :: pre_dif(:,:)
real(KIND=8), allocatable   :: tmp_mat(:,:)
real(KIND=8), allocatable   :: tmp_vec(:,:)
real(KIND=8), allocatable   :: grd_jvc(:,:)
real(KIND=8)                :: jvc_one(1,1)
real(KIND=8)                :: jvc_two(1,1)
real(KIND=8)                :: jvc_the(1,1)
real(KIND=8)                :: bvc_hlf(1,1)
integer                     :: bkg_len
integer                     :: i,j,nitr,mxit
real(KIND=8)                :: jold,jnew,jthr,alph 

bkg_len=size(bkg_cov,1)
allocate (pre_tra(      1,bkg_len))
allocate (pre_bkg(bkg_len,      1))
allocate (pre_anl(bkg_len,      1))
allocate (pre_dif(bkg_len,      1))
allocate (tmp_mat(      1,bkg_len))
allocate (tmp_vec(bkg_len,      1))
allocate (grd_jvc(bkg_len,      1))

nitr = 0
mxit = 500
jold = 100.0 
jnew = 0.0
jthr = 0.0001 
alph = 0.01       
print *,"SOLVER"

40 FORMAT(A8,I4,3F10.4)

! DO THE PRE-CONDITIONING OF X and Xb
anl_vec=bkg_vec
do i=1,bkg_len
    pre_bkg(i,1)=bkg_vec(i) 
    pre_anl(i,1)=anl_vec(i) 
end do
call pre_con_dif(bkg_cov,pre_bkg)
call pre_con_dif(bkg_cov,pre_anl)

! ITERATE TO SOLVE THE COST FUNCTION
do while ( abs(jold-jnew) > jthr)
    if (nitr.gt.mxit) exit
    if (jnew.lt.0.0) exit
    jold=jnew
    jnew=0.0

    do i=1,bkg_len
        pre_dif(i,1)=pre_anl(i,1)-pre_bkg(i,1)
    end do
    pre_tra=transpose(pre_dif)   
   
!   SOLVE FOR COST FUNCTION J
!   FIRST TERM
    jvc_one=matmul(pre_tra,pre_dif)
!   SECOND TERM
    tmp_mat=matmul(pre_tra,hrh_cov)
    jvc_two=matmul(tmp_mat,pre_dif)
!   THIRD TERM
    jvc_the=matmul(pre_tra,bht_ino)
!   COST FUNCTION
    jnew = 0.5*(jvc_one(1,1)+jvc_two(1,1)-2.0*jvc_the(1,1)+jvc_for(1,1))

!   SOLVE FOR GRADIENT OF COST FUNCTION
    tmp_vec=matmul(hrh_cov,pre_dif)

    do i=1,bkg_len
        grd_jvc(i,1)= pre_dif(i,1)+tmp_vec(i,1)-bht_ino(i,1)
    end do

    do i=1,bkg_len
        pre_anl(i,1)= pre_anl(i,1)-grd_jvc(i,1)*alph
    end do

    if (nitr.eq.0) print *,'initial cost = ',jnew
    nitr = nitr + 1 

end do
print *,'final cost = ',jnew,' after ',nitr,' iterations'

! UNDO THE PRE-CONDITIONING X=B(1/2)V

tmp_vec=matmul(bkg_cov,pre_anl)
    
do i=1,bkg_len
    anl_vec(i)=tmp_vec(i,1)
end do

print *,"SOLVER COMPLETE"
end subroutine var_solver

! BJE
! OUTPUT THE ANALYSIS VECTOR
subroutine put_anl_vec(anl_vec,bkg_vec)
real(KIND=8), intent(in)    :: anl_vec(:) 
real(KIND=8), intent(in)    :: bkg_vec(:) 
integer                     :: bkg_len
integer                     :: i
print *,"PUT_ANL_VEC"

40 FORMAT(A8,I4,3F10.4)
bkg_len=size((bkg_vec))
do i=1,bkg_len
    write(*,40) "FIN",i,anl_vec(i),bkg_vec(i),50.0+50.0*sin(((0.2+float(i))/10.0)*(3.1416))
end do

print *,"PUT_ANL_VEC COMPLETE"
end subroutine put_anl_vec
end program adept
