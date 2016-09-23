! BJE = 01 AUG 2016
! PROGRAM TO DO VARIATIONAL ASSIMILATION

program adept

  use gptl
  use varsolver
  use background, only   : Background_Type
  use observations, only : Observations_Type
  use background_covariance, only : Background_Covariance_Type
  use observation_covariance, only : Observation_Covariance_Type

  implicit none

  type(Background_Type)             :: bkg
  type(Observations_Type)           :: obs
  type(Background_Covariance_Type)  :: bkg_cov
  type(Observation_Covariance_Type) :: obs_cov
  real(KIND=8), allocatable    ::      obs_opr(:,:,:)
  real(KIND=8), allocatable    ::      hrh_cov(:,:,:)
  real(KIND=8), allocatable    ::      anl_vec(:,:)
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

  ! Initialize a background object, Xb
  bkg = Background_Type(bkg_len, tim_len, mthd)

  ! OBTAIN THE B(1/2) MATRIX - FOR PRE CONDITIONING
  bkg_cov = Background_Covariance_Type(bkg, tim_len)

  ! Initialize observations object, Y
  obs = Observations_Type(obs_len, mthd)

  ! OBTAIN THE COVARIANCE MATRIX R - OBS ERROR
  ! SOMEDAY WILL COME TO US FROM THE UFO
  obs_cov = Observation_Covariance_Type(obs, tim_len)

  ! BJE
  ! KNOWING THE NUMBERS, ALLOCATE VECTORS/MATRICIES (ARRAYS) ACCORTINGLY
  allocate (obs_opr(tim_len,obs_len,bkg_len))
  allocate (hrh_cov(tim_len,bkg_len,bkg_len))
  allocate (anl_vec(tim_len,bkg_len))
  allocate (bht_ino(tim_len,bkg_len,1))

  ! BJE
  ! GET THE INNOVATION VECTOR - (Y-HXb) - OVERWRITE OBS_VEC
  call get_ino_vec(bkg, obs)

  ! BJE
  ! KNOWING THE LOCATION OF THE OBS, CREATE OBS OPERATOR, H 
  call get_obs_opr(obs, obs_opr)

  ! BJE
  ! GET THE NEEDED REUSED MATRIX PRODUCTS:
  ! B(1/2)(Y-HXb)(T)R(-1)H        = bht_ino
  ! B(1/2)      H(T)R(-1)HB(1/2)  = hrh_cov
  ! INPUTS: Y, H, Xb, R(-1)
  ! USE SAME VARIABLE NAME PRE AND POST
  call pre_sol(obs_opr, obs_cov, bkg_cov, hrh_cov, obs, bht_ino, jvc_for)

  ! BJE
  ! THE MAIN EVENT - THE SOLVER
  call var_solver(bkg_cov, hrh_cov, bht_ino, jvc_for, bkg, anl_vec)

  ! BJE
  ! OUTPUT THE NEW ANALYSIS
  call put_anl_vec(anl_vec, bkg)

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
