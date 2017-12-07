module Abstract_Model

  implicit none

  private

  public :: abstract_model_type

  type, abstract :: abstract_model_type
      private
  contains
    procedure(generic_adv_nsteps), deferred :: adv_nsteps
  end type abstract_model_type

  abstract interface
    subroutine generic_adv_nsteps(this, nsteps)
      import abstract_model_type
      class(abstract_model_type), intent(inout) :: this
      integer                   , intent(   in) :: nsteps
    end subroutine generic_adv_nsteps
  end interface

contains
  

end module Abstract_Model
