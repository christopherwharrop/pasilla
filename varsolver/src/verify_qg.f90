! BJE = 15 MAY 2017
! PROGRAM TO DO VERIFICATION OF ARIATIONAL ASSIMILATION

program adept

  use gptl
  use module_verify_qg
  use Abstract_Model, only : abstract_model_type

  implicit none
  
  real(KIND=8), allocatable    ::      tru_pos(:,:)
  real(KIND=8), allocatable    ::      mod_pos(:,:)
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
  allocate (tru_interp(mod_len))
 
  ! CWH
  ! INITIALIZE THE ALLOCATED VECTORS/MATRICES
  tru_interp(:)  = 0.0

  ! BJE
  ! KNOWING THE LOCATION OF THE OBS, CREATE OBS OPERATOR, H 
  call get_mod_opr(tru_vec,mod_pos,tru_pos,tru_interp)

  ! BJE
  ! GET THE INNOVATION VECTOR - (Y-HXb) - OVERWRITE OBS_VEC
  print *,"MOD = ",mod_vec(1:3)
  print *,"TRU = ",tru_interp(1:3)
  mod_vec=mod_vec-tru_interp
  mod_vec(:)=mod_vec(:)/100000.0d0
  ! BJE
  ! PRINT OUT THE MEAN ABSOLUTE ERROR AND MEAN SQUARED ERROR
  print *,"DIF = ",mod_vec(1:3)

  open(4,file="verify.txt",form="formatted")
  write(4,*) "MAE ALL = ",sum(mod_vec)/float(mod_len)
  write(4,*) "MSE ALL = ",dot_product(mod_vec,mod_vec)/float(mod_len)
  write(4,*) "RMSE ALL = ",sqrt(dot_product(mod_vec,mod_vec)/float(mod_len))

  write(4,*) "MAE 200 = ",sum(mod_vec(1:mod_len/3))/float(mod_len/3)
  write(4,*) "MSE 200 = ",dot_product(mod_vec,mod_vec)/float(mod_len/3)
  write(4,*) "RMSE 200 = ",sqrt(dot_product(mod_vec,mod_vec)/float(mod_len/3))

  write(4,*) "MAE 500 = ",sum(mod_vec(mod_len/3+1:2*mod_len/3))/float(mod_len/3)
  write(4,*) "MSE 500 = ",dot_product(mod_vec(mod_len/3+1:2*mod_len/3),mod_vec(mod_len/3+1:2*mod_len/3))/float(mod_len/3)
  write(4,*) "RMSE 500 = ",sqrt(dot_product(mod_vec(mod_len/3+1:2*mod_len/3),mod_vec(mod_len/3+1:2*mod_len/3))/float(mod_len/3))

  write(4,*) "MAE 800 = ",sum(mod_vec(2*mod_len/3+1:mod_len))/float(mod_len/3)
  write(4,*) "MSE 800 = ",dot_product(mod_vec(2*mod_len/3+1:mod_len),mod_vec(2*mod_len/3+1:mod_len))/float(mod_len/3)
  write(4,*) "RMSE 800 = ",sqrt(dot_product(mod_vec(2*mod_len/3+1:mod_len),mod_vec(2*mod_len/3+1:mod_len))/float(mod_len/3))

  close (4)

  ! BJE
  ! END THE TIMER AND OUTPUT THE GPTL RESULTS
  ret = gptlstop ('adept') 
  ret = gptlpr (0) 
  !ret = gptlpr_summary (MPI_COMM_WORLD) 
  ret = gptlfinalize ()

  ! THAT IS THE END OF THE MAIN PROGRAM
  ! NO, REALLY, THAT WAS THE END OF THE MAIN PROGRAM
  ! SERIOUSLY, THERE IS NOTHING MORE IN THE MAIN PROGRAM
  ! ON TO THE SUBROUTINES!
  !
  !--------------------------------------------------------------------

end program adept
