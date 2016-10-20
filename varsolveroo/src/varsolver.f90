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
  use model, only                  : model_type
  use lorenz96, only               : lorenz96_type, lorenz96_TL_type, lorenz96_ADJ_type

  implicit none

  type(Background_Type)             :: bkg
  type(Observations_Type)           :: obs
  type(Background_Covariance_Type)  :: bkg_cov
  type(Observation_Covariance_Type) :: obs_cov
  type(Innovation_Vector_Type)      :: inno_vec
  type(Observation_Operator_Type)   :: obs_opr
  real(KIND=8), allocatable    ::      hrh_cov(:,:,:)
  real(KIND=8), allocatable    ::      brh_cov(:,:,:)
  real(KIND=8), allocatable    ::      anl_vec(:,:)
  real(KIND=8), allocatable    ::      bht_ino(:,:,:)
  real(KIND=8), allocatable    ::      htr_ino(:,:,:)
  real(KIND=8)                 ::      jvc_for(1,1)
  real                         ::      ret 

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

  ! BJE
  ! KNOWING THE NUMBERS, ALLOCATE VECTORS/MATRICIES (ARRAYS) ACCORTINGLY
  allocate (hrh_cov(tim_len, bkg%npoints, bkg%npoints))
  allocate (brh_cov(tim_len, bkg%npoints, bkg%npoints))
  allocate (anl_vec(tim_len, bkg%npoints))
  allocate (bht_ino(tim_len, bkg%npoints, 1))
  allocate (htr_ino(tim_len, bkg%npoints, 1))

  ! BJE
  ! GET THE NEEDED REUSED MATRIX PRODUCTS:
  !        H(T)R(-1)(Y-HXb) = htr_ino
  !        H(T)R(-1)H       = hrh_cov
  ! B(-1/2)H(T)R(-1)(Y-HXb) = btr_ino
  ! B(-1/2)H(T)R(-1)H       = brh_cov
  ! INPUTS: Y, H, Xb, R(-1), B
  ! USE SAME VARIABLE NAME PRE AND POST
  call pre_sol(obs_opr, obs_cov, bkg_cov, hrh_cov, brh_cov, inno_vec, htr_ino, bht_ino, jvc_for)

  ! BJE
  ! THE MAIN EVENT - THE SOLVER
  call var_solver(bkg_cov, hrh_cov, brh_cov, htr_ino, bht_ino, jvc_for, bkg,     anl_vec, mthd)
!  call var_solver(bkg_cov, hrh_cov, brh_cov, htr_ino, bht_ino, jvc_for, bkg_vec, anl_vec, fwmod_vec, bwmod_vec)


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
