module QG_Util

  use kind

  implicit none

  private

  public :: copy, sub, addab, subab, add

  ! NOTE - SWITCHED VARIABLE "GJACOB" TO "ININAG", AS THEY MATCH

contains


  !-------------------------------------------------------------------------------
  ! Computes the sum of aa and bb, the result is put in aa.  
  !-------------------------------------------------------------------------------
  subroutine addab(aa, bb)

    real(r8kind), intent(inout) ::  aa(:)
    real(r8kind), intent(   in) ::  bb(:)

    integer :: k

    do k = 1, size(aa)
      aa(K) = aa(k) + bb(K)
    end do

  end subroutine


  !-------------------------------------------------------------------------------
  ! Aubtracts bb from aa, the result is put in aa.
  !-------------------------------------------------------------------------------
  subroutine subab(aa,bb)

    real(r8kind), intent(inout) :: aa(:)
    real(r8kind), intent(   in) :: bb(:)

    integer :: k

    do k = 1, size(aa)
      aa(k) = aa(k) - bb(k)
    end do

  end subroutine


  !-------------------------------------------------------------------------------
  ! Negate aa
  !-------------------------------------------------------------------------------
  subroutine neg(aa)

    real(r8kind), intent(inout) :: aa(:)

    integer :: k

    do k = 1, size(aa)
      aa(k) = -aa(k)
    end do

  end subroutine


  !-------------------------------------------------------------------------------
  ! Computes the sum of aa and bb, the result is put in CC.
  !-------------------------------------------------------------------------------
  subroutine add(aa, bb, cc)

    real(r8kind), intent( in) :: aa(:)
    real(r8kind), intent( in) :: bb(:)
    real(r8kind), intent(out) :: cc(:)

    integer :: k

    do k = 1, size(aa)
      cc(k) = aa(k) + bb(k)
    end do

  end subroutine

  !-------------------------------------------------------------------------------
  ! Computes the sum of aa and bb, the result is put in cc.
  !-------------------------------------------------------------------------------
  subroutine add3(aa, bb, cc)

    real(r8kind), intent( in) :: aa(:, :)
    real(r8kind), intent( in) :: bb(:, :)
    real(r8kind), intent(out) :: cc(:, :)

    integer k, l

    do l = 1, 3
      do k = 1, size(aa)
        cc(k, l) = aa(k, l) + bb(k, l)
      enddo
    enddo

  end subroutine



  !-------------------------------------------------------------------------------
  ! Subtracts bb from aa, the result is put in CC.
  !-------------------------------------------------------------------------------
  subroutine sub(aa, bb, CC)

    real(r8kind), intent( in) :: aa(:)
    real(r8kind), intent( in) :: bb(:)
    real(r8kind), intent(out) :: cc(:)

    integer :: k

    do k = 1, size(aa)
      cc(k) = aa(k) - bb(k)
    end do

  end subroutine


  !-------------------------------------------------------------------------------
  ! Subtracts bb from aa, the result is put in cc.
  !-------------------------------------------------------------------------------
  subroutine sub3(aa, bb, cc)

    real(r8kind), intent( in) :: aa(:,:)
    real(r8kind), intent( in) :: bb(:,:)
    real(r8kind), intent(out) :: cc(:,:)

    integer :: k, l

    do l = 1, 3
      do k = 1, size(aa, 1)
        cc(k, l) = aa(k, l) - bb(k, l)
      enddo
    enddo

  end subroutine

  !-------------------------------------------------------------------------------
  ! Computes the product of aa with real c, the result is put in bb.
  !-------------------------------------------------------------------------------
  subroutine mult(aa, C, bb)

    real(r8kind), intent( in) ::  aa(:)
    real(r8kind), intent( in) ::  c
    real(r8kind), intent(out) ::  bb(:)

    integer :: k

    do k = 1, size(aa)
      bb(k) = aa(k) * c
    end do

  end subroutine

  !-------------------------------------------------------------------------------
  ! Computes the product of aa with real C, the result is put in bb.
  !-------------------------------------------------------------------------------
  subroutine mult3(aa, c, bb)

    real(r8kind), intent( in) :: aa(:,:)
    real(r8kind), intent( in) :: c
    real(r8kind), intent(out) :: bb(:,:)

    integer k, l

    do l = 1, 3
      do k = 1, size(aa, 1)
        bb(k, l) = aa(k, l) * c
      enddo
    enddo

  end subroutine

  !-------------------------------------------------------------------------------
  ! copies aa(1 level) TO bb(1 level)
  !-------------------------------------------------------------------------------
  subroutine cop(aa, bb)

    real(r8kind), intent( in) :: aa(:)
    real(r8kind), intent(out) :: bb(:)

    integer :: k

    do k = 1, size(aa)
      bb(k) = aa(k)
    end do

  end subroutine

  !-------------------------------------------------------------------------------
  ! Copies aa into bb.
  !-------------------------------------------------------------------------------
  subroutine copy(aa, bb)

    real(r8kind), intent( in) :: aa(:,:)
    real(r8kind), intent(out) :: bb(:,:)

    integer :: k, l

    do l = 1, 3
      do k = 1, size(aa, 1)
        bb(k, l) = aa(k, l)
      end do
    end do

  end subroutine


end module QG_Util
