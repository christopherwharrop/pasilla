!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of c06fqq in reverse (adjoint) mode:
!   gradient     of useful results: a
!   with respect to varying inputs: a
SUBROUTINE C06FQQ_B(a, ab, m, n)
  IMPLICIT NONE
!VD$R VECTOR
! CVD$R NOLSTVAL
! c CVD$R STRIP
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!     .. Scalar Arguments ..
  INTEGER :: m, n
!     .. Array Arguments ..
  DOUBLE PRECISION :: a(0:m-1, 0:n-1)
  DOUBLE PRECISION :: ab(0:m-1, 0:n-1)
!     .. Local Scalars ..
  INTEGER :: l
!     .. Intrinsic Functions ..
  INTRINSIC MOD
  IF (MOD(n, 2) .EQ. 0) THEN
    DO l=m-1,0,-1
      ab(l, n/2) = 0.5d0*ab(l, n/2)
    END DO
  END IF
  DO l=m-1,0,-1
    ab(l, 0) = 0.5d0*ab(l, 0)
  END DO
END SUBROUTINE C06FQQ_B
