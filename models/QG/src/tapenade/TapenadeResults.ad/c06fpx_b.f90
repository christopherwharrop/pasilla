!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of c06fpx in reverse (adjoint) mode:
!   gradient     of useful results: a b
!   with respect to varying inputs: a b
SUBROUTINE C06FPX_B(a, ab, b, bb, m, n, q, nq, trig)
  IMPLICIT NONE
!VD$R VECTOR
! CVD$R NOLSTVAL
! c CVD$R STRIP
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     Real to Hermitian Fast Fourier Transform Kernel Driver
!
!     Mixed-radix, self-sorting, decimation in time
!
!     .. Scalar Arguments ..
  INTEGER :: m, n, nq
!     .. Array Arguments ..
  DOUBLE PRECISION :: a(0:m*n-1), b(0:m*n-1), trig(0:2*n-1)
  DOUBLE PRECISION :: ab(0:m*n-1), bb(0:m*n-1)
  INTEGER :: q(nq)
!     .. Local Scalars ..
  DOUBLE PRECISION :: factor
  INTEGER :: i, p, qi, r
  LOGICAL :: ina
!     .. External Subroutines ..
  EXTERNAL C06FPR, C06FPS, C06FPT, C06FPU, C06FPV, C06FPW
  EXTERNAL C06FPR_B, C06FPS_B, C06FPT_B, C06FPU_B, C06FPV_B, C06FPW_B
!     .. Intrinsic Functions ..
  INTRINSIC SQRT, DBLE
  INTEGER :: arg1
  INTEGER :: branch
!     .. Executable Statements ..
  ina = .true.
  p = n
  r = 1
  IF (n .NE. 1) THEN
    DO i=nq,1,-1
      qi = q(i)
      CALL PUSHINTEGER4(p)
      p = p/qi
      IF (ina) THEN
        IF (qi .EQ. 2) THEN
          CALL PUSHINTEGER4(arg1)
          arg1 = m*p
          CALL PUSHCONTROL4B(0)
        ELSE IF (qi .EQ. 3) THEN
          CALL PUSHINTEGER4(arg1)
          arg1 = m*p
          CALL PUSHCONTROL4B(1)
        ELSE IF (qi .EQ. 4) THEN
          CALL PUSHINTEGER4(arg1)
          arg1 = m*p
          CALL PUSHCONTROL4B(2)
        ELSE IF (qi .EQ. 5) THEN
          CALL PUSHINTEGER4(arg1)
          arg1 = m*p
          CALL PUSHCONTROL4B(3)
        ELSE IF (qi .EQ. 6) THEN
          CALL PUSHINTEGER4(arg1)
          arg1 = m*p
          CALL PUSHCONTROL4B(4)
        ELSE
          CALL PUSHINTEGER4(arg1)
          arg1 = m*p
          CALL PUSHCONTROL4B(5)
        END IF
      ELSE IF (qi .EQ. 2) THEN
        CALL PUSHINTEGER4(arg1)
        arg1 = m*p
        CALL PUSHCONTROL4B(6)
      ELSE IF (qi .EQ. 3) THEN
        CALL PUSHINTEGER4(arg1)
        arg1 = m*p
        CALL PUSHCONTROL4B(7)
      ELSE IF (qi .EQ. 4) THEN
        CALL PUSHINTEGER4(arg1)
        arg1 = m*p
        CALL PUSHCONTROL4B(8)
      ELSE IF (qi .EQ. 5) THEN
        CALL PUSHINTEGER4(arg1)
        arg1 = m*p
        CALL PUSHCONTROL4B(9)
      ELSE IF (qi .EQ. 6) THEN
        CALL PUSHINTEGER4(arg1)
        arg1 = m*p
        CALL PUSHCONTROL4B(10)
      ELSE
        CALL PUSHINTEGER4(arg1)
        arg1 = m*p
        CALL PUSHCONTROL4B(11)
      END IF
      ina = .NOT.ina
      CALL PUSHINTEGER4(r)
      r = r*qi
    END DO
!
    factor = 1.0d0/SQRT(DBLE(n))
    IF (ina) THEN
      DO i=m*n-1,0,-1
        ab(i) = factor*ab(i)
      END DO
    ELSE
      DO i=m*n-1,0,-1
        bb(i) = bb(i) + factor*ab(i)
        ab(i) = 0.D0
      END DO
    END IF
    DO i=1,nq,1
      qi = q(i)
      CALL POPINTEGER4(r)
      CALL POPCONTROL4B(branch)
      IF (branch .LT. 6) THEN
        IF (branch .LT. 3) THEN
          IF (branch .EQ. 0) THEN
            CALL C06FPW_B(a, ab, b, bb, arg1, r, trig((p-1)*qi*r), trig(&
&                   n+(p-1)*qi*r))
            CALL POPINTEGER4(arg1)
          ELSE IF (branch .EQ. 1) THEN
            CALL C06FPV_B(a, ab, b, bb, arg1, r, trig((p-1)*qi*r), trig(&
&                   n+(p-1)*qi*r))
            CALL POPINTEGER4(arg1)
          ELSE
            CALL C06FPU_B(a, ab, b, bb, arg1, r, trig((p-1)*qi*r), trig(&
&                   n+(p-1)*qi*r))
            CALL POPINTEGER4(arg1)
          END IF
        ELSE IF (branch .EQ. 3) THEN
          CALL C06FPT_B(a, ab, b, bb, arg1, r, trig((p-1)*qi*r), trig(n+&
&                 (p-1)*qi*r))
          CALL POPINTEGER4(arg1)
        ELSE IF (branch .EQ. 4) THEN
          CALL C06FPS_B(a, ab, b, bb, arg1, r, trig((p-1)*qi*r), trig(n+&
&                 (p-1)*qi*r))
          CALL POPINTEGER4(arg1)
        ELSE
          CALL C06FPR_B(a, ab, b, bb, arg1, qi, r, trig((p-1)*qi*r), &
&                 trig(n+(p-1)*qi*r))
          CALL POPINTEGER4(arg1)
        END IF
      ELSE IF (branch .LT. 9) THEN
        IF (branch .EQ. 6) THEN
          CALL C06FPW_B(b, bb, a, ab, arg1, r, trig((p-1)*qi*r), trig(n+&
&                 (p-1)*qi*r))
          CALL POPINTEGER4(arg1)
        ELSE IF (branch .EQ. 7) THEN
          CALL C06FPV_B(b, bb, a, ab, arg1, r, trig((p-1)*qi*r), trig(n+&
&                 (p-1)*qi*r))
          CALL POPINTEGER4(arg1)
        ELSE
          CALL C06FPU_B(b, bb, a, ab, arg1, r, trig((p-1)*qi*r), trig(n+&
&                 (p-1)*qi*r))
          CALL POPINTEGER4(arg1)
        END IF
      ELSE IF (branch .EQ. 9) THEN
        CALL C06FPT_B(b, bb, a, ab, arg1, r, trig((p-1)*qi*r), trig(n+(p&
&               -1)*qi*r))
        CALL POPINTEGER4(arg1)
      ELSE IF (branch .EQ. 10) THEN
        CALL C06FPS_B(b, bb, a, ab, arg1, r, trig((p-1)*qi*r), trig(n+(p&
&               -1)*qi*r))
        CALL POPINTEGER4(arg1)
      ELSE
        CALL C06FPR_B(b, bb, a, ab, arg1, qi, r, trig((p-1)*qi*r), trig(&
&               n+(p-1)*qi*r))
        CALL POPINTEGER4(arg1)
      END IF
      CALL POPINTEGER4(p)
    END DO
  END IF
END SUBROUTINE C06FPX_B