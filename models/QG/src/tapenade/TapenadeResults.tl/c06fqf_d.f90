!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of c06fqf in forward (tangent) mode:
!   variations   of useful results: x work
!   with respect to varying inputs: x work
SUBROUTINE C06FQF_D(m, n, x, xd, init, trig, work, workd, ifail)
  IMPLICIT NONE
!VD$R NOVECTOR
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!     .. Parameters ..
  CHARACTER(len=6) :: srname
  PARAMETER (srname='C06FQF')
!     .. Scalar Arguments ..
  INTEGER :: ifail, m, n
  CHARACTER(len=1) :: init
!     .. Array Arguments ..
  DOUBLE PRECISION :: trig(2*n), work(m*n), x(m*n)
  DOUBLE PRECISION :: workd(m*n), xd(m*n)
!     .. Local Scalars ..
  INTEGER :: ierror, nq, nrec
!     .. Local Arrays ..
  INTEGER :: q(30)
  CHARACTER(len=80) :: rec(1)
  EXTERNAL P01ABF
!     .. External Functions ..
  INTEGER :: P01ABF
!     .. External Subroutines ..
  EXTERNAL C06FPQ, C06FQX
  EXTERNAL C06FQX_D
!     .. Executable Statements ..
  CALL C06FPQ(m, n, init, trig, q, nq, ierror)
  IF (ierror .EQ. 0) CALL C06FQX_D(x, xd, work, workd, m, n, q, nq, trig&
&                           )
!CWH      ELSE IF (IERROR.EQ.1) THEN
!CWH         WRITE (REC(1),FMT=99999) M
!CWH      ELSE IF (IERROR.EQ.2) THEN
!CWH         WRITE (REC(1),FMT=99998) N
!CWH      ELSE IF (IERROR.EQ.3) THEN
!CWH         WRITE (REC(1),FMT=99997) INIT
!CWH      ELSE IF (IERROR.EQ.4) THEN
!CWH         WRITE (REC(1),FMT=99996) INIT
!CWH      ELSE IF (IERROR.EQ.5) THEN
!CWH         WRITE (REC(1),FMT=99995) INIT
!
  nrec = 1
  ifail = P01ABF(ifail, ierror, srname, nrec, rec)
!
  RETURN
!
99999 FORMAT(' ** M must be at least 1: M = ',i16)
99998 FORMAT(' ** N must be at least 1: N = ',i16)
99997 FORMAT(' ** ',a1,' is an invalid value of INIT')
99996 FORMAT(' ** INIT = ',a1,', but TRIG array never initialized')
99995 FORMAT(' ** INIT = ',a1,', but N and TRIG array incompatible')
END SUBROUTINE C06FQF_D
