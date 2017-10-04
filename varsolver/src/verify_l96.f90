! BJE = 15 MAY 2017
! PROGRAM TO DO VERIFICATION OF ARIATIONAL ASSIMILATION

program adept

  use gptl
  use module_verify_l96
  use Abstract_Model, only : abstract_model_type

  implicit none
  
  real(KIND=8), allocatable    ::      tru_pos(:)
  real(KIND=8), allocatable    ::      mod_pos(:)
  real(KIND=8), allocatable    ::      tru_interp(:)
  real(KIND=8), allocatable    ::      mod_vec(:)
  real(KIND=8), allocatable    ::      tru_vec(:)
  real                         ::      ret 

  ! BJE
  ! INITIALIZE GPTL AND START A TIMER
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('adept')

  ! BJE
  ! FIRST - NEED TO KNOW HOW MANY OBSERVATIONS AND STATE VECTOR
  ! OBTAIN THE OBSERATIONS, Y, AND THE BACKGROUND, Xb 
  call get_mod_vec(mod_pos,mod_vec)
  print *,"HAVE THE MODEL"

  call get_tru_vec(tru_pos,tru_vec)
  print *,"HAVE THE TRUTH"

  ! BJE
  ! KNOWING THE NUMBERS, ALLOCATE VECTORS/MATRICIES (ARRAYS) ACCORDINGLY
  mod_len=size(mod_vec)
 
  ! BJE
  ! GET THE INNOVATION VECTOR - (Y-HXb) - OVERWRITE OBS_VEC
  print *,"MOD = ",mod_vec(1:3)
  print *,"TRU = ",tru_vec(1:3)
  mod_vec = mod_vec - tru_vec

  ! BJE
  ! PRINT OUT THE MEAN ABSOLUTE ERROR AND MEAN SQUARED ERROR
  print *,"DIF = ",mod_vec(1:3)

  open(40,file="verify.txt",form="formatted", status="new")
  write(40,*) "MAE ALL = ",sum(mod_vec)/float(mod_len)
!  write(40,*) "MSE ALL = ",dot_product(mod_vec,mod_vec)/float(mod_len)
  write(40,*) "RMSE ALL = ",sqrt(dot_product(mod_vec,mod_vec)/float(mod_len))
  close (40)

  ! BJE
  ! END THE TIMER AND OUTPUT THE GPTL RESULTS
  ret = gptlstop ('adept') 
  ret = gptlpr (0) 
  ret = gptlfinalize ()

  ! THAT IS THE END OF THE MAIN PROGRAM
  ! NO, REALLY, THAT WAS THE END OF THE MAIN PROGRAM
  ! SERIOUSLY, THERE IS NOTHING MORE IN THE MAIN PROGRAM
  ! ON TO THE SUBROUTINES!
  !
  !--------------------------------------------------------------------

end program adept
