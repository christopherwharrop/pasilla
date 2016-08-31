module module_kind

  implicit none
  private 

  public :: i4kind, i8kind, r4kind, r8kind, ckind

  integer, parameter :: i4kind = SELECTED_INT_KIND(8)
  integer, parameter :: i8kind = SELECTED_INT_KIND(13)
  integer, parameter :: r4kind = SELECTED_REAL_KIND(6)
  integer, parameter :: r8kind = SELECTED_REAL_KIND(12)
  integer, parameter :: ckind  = SELECTED_CHAR_KIND('DEFAULT')

end module module_kind

