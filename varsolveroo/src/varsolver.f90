! BJE = 01 AUG 2016
! PROGRAM TO DO VARIATIONAL ASSIMILATION

program adept

  use gptl
  use varsolver
  use background, only             : Background_Type
  use observations, only           : Observations_Type
  use background_covariance, only  : Background_Covariance_Type
  use observation_covariance, only : Observation_Covariance_Type
  use innovation_vector, only      : Innovation_Vector_Type
  use observation_operator, only   : Observation_Operator_Type

  implicit none

  type(Background_Type)             :: bkg
  type(Observations_Type)           :: obs
  type(Background_Covariance_Type)  :: bkg_cov
  type(Observation_Covariance_Type) :: obs_cov
  type(Innovation_Vector_Type)      :: inno_vec
  type(Observation_Operator_Type)   :: obs_opr
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

  ! KNOWING THE LOCATION OF THE OBS, CREATE OBS OPERATOR, H
  obs_opr = Observation_Operator_Type(bkg, obs)

  ! GET THE INNOVATION VECTOR - (Y-HXb) - DO NOT OVERWRITE OBS
  inno_vec = Innovation_Vector_Type(bkg, obs)

  ! BJE
  ! KNOWING THE NUMBERS, ALLOCATE VECTORS/MATRICIES (ARRAYS) ACCORTINGLY
  allocate (hrh_cov(tim_len,bkg_len,bkg_len))
  allocate (anl_vec(tim_len,bkg_len))
  allocate (bht_ino(tim_len,bkg_len,1))

  ! BJE
  ! GET THE NEEDED REUSED MATRIX PRODUCTS:
  ! B(1/2)(Y-HXb)(T)R(-1)H        = bht_ino
  ! B(1/2)      H(T)R(-1)HB(1/2)  = hrh_cov
  ! INPUTS: Y, H, Xb, R(-1)
  ! USE SAME VARIABLE NAME PRE AND POST
  call pre_sol(obs_opr, obs_cov, bkg_cov, hrh_cov, inno_vec, bht_ino, jvc_for)

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
