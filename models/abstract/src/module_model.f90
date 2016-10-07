module model

  use kind, only : r8kind

  implicit none

  private

  public :: model_type

  type, abstract ::  model_type
      private
!      integer, public :: size
!      real(r8kind) :: delta_t
!      real(r8kind) :: t
!      integer, public :: step
!      real(r8kind), allocatable, public :: state(:)
!      real(r8kind), allocatable, public :: trajectory(:)
!      real(r8kind), allocatable, public :: location(:)
  contains
!      procedure(adv_nsteps), deferred :: adv_nsteps
  end type model_type

!  abstract interface
!    subroutine adv_nsteps(this,nsteps)
!      import model_type
!      class(model_type) :: this
!      integer :: nsteps
!    end subroutine adv_nsteps
!  end interface

contains

  

end module model
