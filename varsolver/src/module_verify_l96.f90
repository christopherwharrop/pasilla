! BJE = 12 MAY 2017
! PROGRAM TO DO VERIFICATION OF VARIATIONAL ASSIMILATION

module module_verify_l96

  use module_constants
  use Abstract_Model, only : abstract_model_type
  use L96_Config,     only : l96_config_type
  use L96_Model,      only : l96_model_type
  use L96_Reader,     only : l96_reader_type
  use L96_Writer,     only : l96_writer_type

  ! Get unit numbers for stdin, stdout, stderr in a portable way
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  implicit none

  ! Define namelists and default values
  integer          :: mthd  = 9
  integer          :: tim_len = 1
  type(l96_config_type) :: tru_config
  integer          :: tru_len
  real(KIND=8)     :: tru_nx, tru_ny, tru_nz

  type(l96_config_type) :: mod_config
  integer          :: mod_len
  real(KIND=8)     :: mod_nx, mod_ny, mod_nz

  real(KIND=8)     :: alph
  real(KIND=8)     :: sigma
  namelist /control/ mthd, tim_len
  namelist /method1/  alph, sigma
  namelist /method2/  alph, sigma
  namelist /method3/  alph, sigma
  namelist /method4/  alph, sigma

contains

  ! BJE
  ! GENERATE THE MODERVATIONS "Y", AND THEIR LOCATIONS - FOR "H"
  subroutine get_mod_vec(mod_pos,mod_vec)

    implicit none
    real(KIND=8), intent(inout), allocatable :: mod_vec(:)
    real(KIND=8), intent(inout), allocatable :: mod_pos(:)

    integer                                 :: i,t,tt
    type(l96_config_type)                    :: config
    type(l96_model_type)                     :: model
    type(l96_reader_type)                    :: reader
    character(len=256)                      :: filename

    print *,"GET_MOD_VEC"
    reader = l96_reader_type("NETCDF")
    filename="model"
    call reader%read(model, filename)
    mod_config = model%get_config()
    mod_nx = mod_config%get_nx()
    mod_len = mod_nx
    allocate (mod_vec(mod_len))
    allocate (mod_pos(mod_len))
    mod_pos(:) = model%get_location()
    mod_vec(:) = model%get_state()

    print *,"GET_MOD_VEC COMPLETE"

  end subroutine get_mod_vec

  ! BJE
  ! GET THE MODERVATION OPERATOR, H, FROM THE INPUTS
  subroutine get_mod_opr(tru_vec,mod_pos,tru_pos,tru_interp)

    implicit none

    real(KIND=8), intent( in) :: tru_vec(:)
    real(KIND=8), intent( in) :: mod_pos(:)
    real(KIND=8), intent( in) :: tru_pos(:)
    real(KIND=8), intent(out) :: tru_interp(:)

    integer :: i
    integer :: lower_index, upper_index
    real(KIND=8) :: lower_weight, upper_weight
    type(l96_model_type) :: model

    print *,"GET_MOD_OPR"

    model = l96_model_type(tru_config)

    do i=1,mod_len
       call model%get_interpolation_weights(mod_pos(i), lower_index, upper_index, lower_weight, upper_weight)
!      call model%get_interpolation_weights(45.0d0,45.0d0,200.0d0,NW_index, NE_index, &
!        SW_index, SE_index,  NW_weight, NE_weight, SW_weight, SE_weight)
       if(i.le.4) then
         print *,"I!",i
         print *,mod_pos(i)
         print *,lower_index, upper_index
         print *,tru_pos(lower_index), tru_pos(upper_index)
       endif
       tru_interp(i) = tru_vec(lower_index) * lower_weight + &
                     & tru_vec(upper_index) * upper_weight
       if(i.le.4) then
         print *,tru_vec(lower_index),tru_vec(upper_index)
         print *,"INTERP = ",tru_interp(i)
         print *," "
       endif
    end do

    print *,"GET_MOD_OPR COMPLETE"

  end subroutine get_mod_opr

  ! BJE
  ! GENERATE THE FIRST GUESS "Xb" - SHOULD USE A REAL MODEL
  subroutine get_tru_vec(tru_pos,tru_vec)

    implicit none
    real(KIND=8), intent(inout), allocatable :: tru_vec(:)
    real(KIND=8), intent(inout), allocatable :: tru_pos(:)

    integer                                 :: i,t,tt
    type(l96_config_type)                    :: config
    type(l96_model_type)                     :: model
    type(l96_reader_type)                    :: reader
    character(len=256)                      :: filename

    print *,"GET_TRU_VEC"
    reader = l96_reader_type("NETCDF")
    filename="truth"
    call reader%read(model, filename)
    tru_config = model%get_config()
    tru_nx = tru_config%get_nx()
    tru_len = tru_nx
    allocate (tru_vec(tru_len))
    allocate (tru_pos(tru_len))
    tru_pos(:) = model%get_location()
    tru_vec(:) = model%get_state()

    print *,"GET_TRU_VEC COMPLETE"

  end subroutine get_tru_vec


  ! THAT IS ALL FOLKS!
  !
end module module_verify_l96
