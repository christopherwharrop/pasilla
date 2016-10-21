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
  use solver, only                 : Solver_Type
  use model, only                  : model_type
  use lorenz96, only               : lorenz96_type, lorenz96_TL_type, lorenz96_ADJ_type

  implicit none

  type(Background_Type)             :: bkg
  type(Observations_Type)           :: obs
  type(Background_Covariance_Type)  :: bkg_cov
  type(Observation_Covariance_Type) :: obs_cov
  type(Innovation_Vector_Type)      :: inno_vec
  type(Observation_Operator_Type)   :: obs_opr
  type(Solver_Type)                 :: solve
  real                              ::      ret 

  integer :: t,i

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
  bkg = Background_Type(tim_len, mthd)

  ! OBTAIN THE B(1/2) MATRIX - FOR PRE CONDITIONING
  bkg_cov = Background_Covariance_Type(bkg, tim_len, sigma)

  ! Initialize observations object, Y
  obs = Observations_Type(mthd)

  ! OBTAIN THE COVARIANCE MATRIX R - OBS ERROR
  ! SOMEDAY WILL COME TO US FROM THE UFO
  obs_cov = Observation_Covariance_Type(obs, tim_len)

  ! KNOWING THE LOCATION OF THE OBS, CREATE OBS OPERATOR, H
  obs_opr = Observation_Operator_Type(bkg, obs)

  ! GET THE INNOVATION VECTOR - (Y-HXb) - DO NOT OVERWRITE OBS
  inno_vec = Innovation_Vector_Type(bkg, obs)

  ! Initialize a solver
  solve = Solver_Type(mthd, bkg)

  ! THE MAIN EVENT - THE SOLVER
  call solve%var_solver(bkg, bkg_cov, obs_cov, obs_opr, inno_vec, mthd, alph)

  ! OUTPUT THE NEW ANALYSIS
  call solve%put_anl_vec(bkg, mthd)

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
