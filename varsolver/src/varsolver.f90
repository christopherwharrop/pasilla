! BJE = 01 AUG 2016
! PROGRAM TO DO VARIATIONAL ASSIMILATION

program adept

  use gptl
  use module_varsolver

  implicit none

  integer, allocatable         ::      bkg_tim(:)
  integer, allocatable         ::      bkg_pos(:,:)
  integer, allocatable         ::      obs_tim(:)
  integer, allocatable         ::      obs_pos(:) 
  real(KIND=8), allocatable    ::      obs_opr(:,:,:)
  real(KIND=8), allocatable    ::      obs_cov(:,:,:)
  real(KIND=8), allocatable    ::      bkg_cov(:,:,:)
  real(KIND=8), allocatable    ::      hrh_cov(:,:,:)
  real(KIND=8), allocatable    ::      brh_cov(:,:,:)
  real(KIND=8), allocatable    ::      obs_vec(:)
  real(KIND=8), allocatable    ::      bkg_vec(:,:)
  real(KIND=8), allocatable    ::      anl_vec(:,:)
  real(KIND=8), allocatable    ::      htr_ino(:,:,:)
  real(KIND=8), allocatable    ::      bht_ino(:,:,:)
  real(KIND=8)                 ::      jvc_for(1,1)
  real                         ::      ret 

  ! BJE
  ! INITIALIZE GPTL AND START A TIMER
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()                    
  ret = gptlstart ('adept')                 

  ! BJE
  ! GET THE METHOD TO USE
  call get_method

  ! BJE
  ! FIRST - NEED TO KNOW HOW MANY OBSERVATIONS AND STATE VECTOR
  ! OBTAIN THE OBSERATIONS, Y, AND THE BACKGROUND, Xb 
  call get_bkg_vec(bkg_tim,bkg_pos,bkg_vec)
  call get_obs_vec(obs_tim,obs_pos,obs_vec)

  ! BJE
  ! KNOWING THE NUMBERS, ALLOCATE VECTORS/MATRICIES (ARRAYS) ACCORTINGLY
  allocate (obs_opr(tim_len,obs_len,bkg_len))
  allocate (obs_cov(tim_len,obs_len,obs_len))
  allocate (bkg_cov(tim_len,bkg_len,bkg_len))
  allocate (hrh_cov(tim_len,bkg_len,bkg_len))
  allocate (brh_cov(tim_len,bkg_len,bkg_len))
  allocate (anl_vec(tim_len,bkg_len))
  allocate (htr_ino(tim_len,bkg_len,1))
  allocate (bht_ino(tim_len,bkg_len,1))

  ! BJE
  ! GET THE INNOVATION VECTOR - (Y-HXb) - OVERWRITE OBS_VEC
  call get_ino_vec(bkg_tim,bkg_pos,obs_tim,obs_pos,obs_vec,bkg_vec)

  ! BJE
  ! KNOWING THE LOCATION OF THE OBS, CREATE OBS OPERATOR, H 
  call get_obs_opr(obs_tim,obs_pos,obs_opr)

  ! BJE   
  ! OBTAIN THE COVARIANCE MATRIX R - OBS ERROR
  ! SOMEDAY WILL COME TO US FROM THE UFO
  call get_obs_cov(obs_tim,obs_pos,obs_cov)

  ! BJE
  ! OBTAIN THE B(1/2) MATRIX - FOR PRE CONDITIONING 
  call get_bkg_cov(bkg_cov)

  ! BJE
  ! GET THE NEEDED REUSED MATRIX PRODUCTS:
  !        H(T)R(-1)(Y-HXb) = htr_ino
  !        H(T)R(-1)H       = hrh_cov
  ! B(-1/2)H(T)R(-1)(Y-HXb) = btr_ino
  ! B(-1/2)H(T)R(-1)H       = brh_cov 
  ! INPUTS: Y, H, Xb, R(-1), B 
  ! USE SAME VARIABLE NAME PRE AND POST
  call pre_sol(obs_opr,obs_cov,bkg_cov,hrh_cov,brh_cov,obs_vec,htr_ino,bht_ino,jvc_for)

  ! BJE
  ! THE MAIN EVENT - THE SOLVER
  call var_solver(bkg_cov,hrh_cov,brh_cov,htr_ino,bht_ino,jvc_for,bkg_vec,anl_vec)

  ! BJE
  ! OUTPUT THE NEW ANALYSIS
  call put_anl_vec(anl_vec,bkg_vec,bkg_tim)

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

end program adept
