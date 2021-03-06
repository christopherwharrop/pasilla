!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of c06fqq in forward (tangent) mode:
!   variations   of useful results: a
!   with respect to varying inputs: a
SUBROUTINE C06FQQ_D(a, ad, m, n)
  IMPLICIT NONE
!VD$R VECTOR
! CVD$R NOLSTVAL
! c CVD$R STRIP
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!     .. Scalar Arguments ..
  INTEGER :: m, n
!     .. Array Arguments ..
  DOUBLE PRECISION :: a(0:m-1, 0:n-1)
  DOUBLE PRECISION :: ad(0:m-1, 0:n-1)
!     .. Local Scalars ..
  INTEGER :: l
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!     .. Executable Statements ..
  DO l=0,m-1
    ad(l, 0) = 0.5d0*ad(l, 0)
    a(l, 0) = 0.5d0*a(l, 0)
  END DO
  IF (MOD(n, 2) .EQ. 0) THEN
    DO l=0,m-1
      ad(l, n/2) = 0.5d0*ad(l, n/2)
      a(l, n/2) = 0.5d0*a(l, n/2)
    END DO
  END IF
!
  RETURN
END SUBROUTINE C06FQQ_D
