!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
!
!  Differentiation of c06fqt in forward (tangent) mode:
!   variations   of useful results: b
!   with respect to varying inputs: a b
SUBROUTINE C06FQT_D(a, ad, b, bd, p, r, cosine, sine)
  IMPLICIT NONE
!VD$R VECTOR
! CVD$R NOLSTVAL
! c CVD$R STRIP
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!     MARK 14 REVISED. IER-700 (DEC 1989).
!
!     Radix five Hermitian to real fast Fourier transform kernel
!
!     Self-sorting, decimation in frequency
!
!     .. Parameters ..
  DOUBLE PRECISION :: r54, sin36, sin72, s36s72
  PARAMETER (r54=0.559016994374947424102293417182819d0, sin36=&
&   0.587785252292473129168705954639073d0, sin72=&
&   0.951056516295153572116439333379382d0, s36s72=&
&   0.618033988749894848204586834365638d0)
!     .. Scalar Arguments ..
  INTEGER :: p, r
!     .. Array Arguments ..
  DOUBLE PRECISION :: a(0:p-1, 0:r-1, 0:4), b(0:p-1, 0:4, 0:r-1), cosine&
& (0:r-1, 4), sine(0:r-1, 4)
  DOUBLE PRECISION :: ad(0:p-1, 0:r-1, 0:4), bd(0:p-1, 0:4, 0:r-1)
!     .. Local Scalars ..
  DOUBLE PRECISION :: c2k, c3k, c4k, ck, s2k, s3k, s4k, sk, t1, t10, &
& t10i, t10r, t11, t11i, t11r, t1i, t1r, t2, t2i, t2r, t3, t3i, t3r, t4&
& , t4i, t4r, t5, t5i, t5r, t6, t6i, t6r, t7, t7i, t7r, t8, t8i, t8r, t9&
& , t9i, t9r, x0p, x1p, x2p, x3p, x4p, y0p, y1p, y2p, y3p, y4p
  DOUBLE PRECISION :: t1d, t10d, t10id, t10rd, t11d, t11id, t11rd, t1id&
& , t1rd, t2d, t2id, t2rd, t3d, t3id, t3rd, t4d, t4id, t4rd, t5d, t5id, &
& t5rd, t6d, t6id, t6rd, t7d, t7id, t7rd, t8d, t8id, t8rd, t9d, t9id, &
& t9rd, x0pd, x1pd, x2pd, x3pd, x4pd, y0pd, y1pd, y2pd, y3pd, y4pd
  INTEGER :: i, k, kp, r2
!     .. Intrinsic Functions ..
  INTRINSIC MOD
!     .. Executable Statements ..
  DO i=0,p-1
!
!        Code for K=0 --
!
    t1d = ad(i, 0, 1)
    t1 = a(i, 0, 1)
    t2d = ad(i, 0, 2)
    t2 = a(i, 0, 2)
    t3d = sin72*ad(i, 0, 4)
    t3 = sin72*a(i, 0, 4)
    t4d = sin72*ad(i, 0, 3)
    t4 = sin72*a(i, 0, 3)
    t5d = t1d + t2d
    t5 = t1 + t2
    t6d = r54*(t1d-t2d)
    t6 = r54*(t1-t2)
    t7d = ad(i, 0, 0) - 0.25d0*t5d
    t7 = a(i, 0, 0) - 0.25d0*t5
    t8d = t7d + t6d
    t8 = t7 + t6
    t9d = t7d - t6d
    t9 = t7 - t6
    t10d = t3d + s36s72*t4d
    t10 = t3 + s36s72*t4
    t11d = s36s72*t3d - t4d
    t11 = s36s72*t3 - t4
    bd(i, 0, 0) = ad(i, 0, 0) + t5d
    b(i, 0, 0) = a(i, 0, 0) + t5
    bd(i, 1, 0) = t8d + t10d
    b(i, 1, 0) = t8 + t10
    bd(i, 2, 0) = t9d + t11d
    b(i, 2, 0) = t9 + t11
    bd(i, 3, 0) = t9d - t11d
    b(i, 3, 0) = t9 - t11
    bd(i, 4, 0) = t8d - t10d
    b(i, 4, 0) = t8 - t10
  END DO
!
!     Code for general K --
!
  IF (p .LE. (r-1)/2) THEN
    DO i=0,p-1
!DIR$ IVDEP
      DO k=1,(r-1)/2
        t1rd = ad(i, k, 1) + ad(i, r-k, 0)
        t1r = a(i, k, 1) + a(i, r-k, 0)
        t1id = ad(i, r-k, 3) - ad(i, k, 4)
        t1i = a(i, r-k, 3) - a(i, k, 4)
        t2rd = ad(i, k, 2) + ad(i, r-k, 1)
        t2r = a(i, k, 2) + a(i, r-k, 1)
        t2id = ad(i, r-k, 2) - ad(i, k, 3)
        t2i = a(i, r-k, 2) - a(i, k, 3)
        t3rd = sin72*(ad(i, k, 1)-ad(i, r-k, 0))
        t3r = sin72*(a(i, k, 1)-a(i, r-k, 0))
        t3id = sin72*(ad(i, r-k, 3)+ad(i, k, 4))
        t3i = sin72*(a(i, r-k, 3)+a(i, k, 4))
        t4rd = sin72*(ad(i, k, 2)-ad(i, r-k, 1))
        t4r = sin72*(a(i, k, 2)-a(i, r-k, 1))
        t4id = sin72*(ad(i, r-k, 2)+ad(i, k, 3))
        t4i = sin72*(a(i, r-k, 2)+a(i, k, 3))
        t5rd = t1rd + t2rd
        t5r = t1r + t2r
        t5id = t1id + t2id
        t5i = t1i + t2i
        t6rd = r54*(t1rd-t2rd)
        t6r = r54*(t1r-t2r)
        t6id = r54*(t1id-t2id)
        t6i = r54*(t1i-t2i)
        t7rd = ad(i, k, 0) - 0.25d0*t5rd
        t7r = a(i, k, 0) - 0.25d0*t5r
        t7id = ad(i, r-k, 4) - 0.25d0*t5id
        t7i = a(i, r-k, 4) - 0.25d0*t5i
        t8rd = t7rd + t6rd
        t8r = t7r + t6r
        t8id = t7id + t6id
        t8i = t7i + t6i
        t9rd = t7rd - t6rd
        t9r = t7r - t6r
        t9id = t7id - t6id
        t9i = t7i - t6i
        t10rd = t3rd + s36s72*t4rd
        t10r = t3r + s36s72*t4r
        t10id = t3id + s36s72*t4id
        t10i = t3i + s36s72*t4i
        t11rd = s36s72*t3rd - t4rd
        t11r = s36s72*t3r - t4r
        t11id = s36s72*t3id - t4id
        t11i = s36s72*t3i - t4i
        x0pd = ad(i, k, 0) + t5rd
        x0p = a(i, k, 0) + t5r
        y0pd = ad(i, r-k, 4) + t5id
        y0p = a(i, r-k, 4) + t5i
        x1pd = t8rd + t10id
        x1p = t8r + t10i
        y1pd = t8id - t10rd
        y1p = t8i - t10r
        x2pd = t9rd + t11id
        x2p = t9r + t11i
        y2pd = t9id - t11rd
        y2p = t9i - t11r
        x3pd = t9rd - t11id
        x3p = t9r - t11i
        y3pd = t9id + t11rd
        y3p = t9i + t11r
        x4pd = t8rd - t10id
        x4p = t8r - t10i
        y4pd = t8id + t10rd
        y4p = t8i + t10r
        bd(i, 0, k) = x0pd
        b(i, 0, k) = x0p
        bd(i, 0, r-k) = y0pd
        b(i, 0, r-k) = y0p
        bd(i, 1, k) = cosine(k, 1)*x1pd - sine(k, 1)*y1pd
        b(i, 1, k) = cosine(k, 1)*x1p - sine(k, 1)*y1p
        bd(i, 1, r-k) = cosine(k, 1)*y1pd + sine(k, 1)*x1pd
        b(i, 1, r-k) = cosine(k, 1)*y1p + sine(k, 1)*x1p
        bd(i, 2, k) = cosine(k, 2)*x2pd - sine(k, 2)*y2pd
        b(i, 2, k) = cosine(k, 2)*x2p - sine(k, 2)*y2p
        bd(i, 2, r-k) = cosine(k, 2)*y2pd + sine(k, 2)*x2pd
        b(i, 2, r-k) = cosine(k, 2)*y2p + sine(k, 2)*x2p
        bd(i, 3, k) = cosine(k, 3)*x3pd - sine(k, 3)*y3pd
        b(i, 3, k) = cosine(k, 3)*x3p - sine(k, 3)*y3p
        bd(i, 3, r-k) = cosine(k, 3)*y3pd + sine(k, 3)*x3pd
        b(i, 3, r-k) = cosine(k, 3)*y3p + sine(k, 3)*x3p
        bd(i, 4, k) = cosine(k, 4)*x4pd - sine(k, 4)*y4pd
        b(i, 4, k) = cosine(k, 4)*x4p - sine(k, 4)*y4p
        bd(i, 4, r-k) = cosine(k, 4)*y4pd + sine(k, 4)*x4pd
        b(i, 4, r-k) = cosine(k, 4)*y4p + sine(k, 4)*x4p
      END DO
    END DO
  ELSE
    DO k=1,(r-1)/2
      kp = r - k
      ck = cosine(k, 1)
      sk = sine(k, 1)
      c2k = cosine(k, 2)
      s2k = sine(k, 2)
      c3k = cosine(k, 3)
      s3k = sine(k, 3)
      c4k = cosine(k, 4)
      s4k = sine(k, 4)
      DO i=0,p-1
        t1rd = ad(i, k, 1) + ad(i, kp, 0)
        t1r = a(i, k, 1) + a(i, kp, 0)
        t1id = ad(i, kp, 3) - ad(i, k, 4)
        t1i = a(i, kp, 3) - a(i, k, 4)
        t2rd = ad(i, k, 2) + ad(i, kp, 1)
        t2r = a(i, k, 2) + a(i, kp, 1)
        t2id = ad(i, kp, 2) - ad(i, k, 3)
        t2i = a(i, kp, 2) - a(i, k, 3)
        t3rd = sin72*(ad(i, k, 1)-ad(i, kp, 0))
        t3r = sin72*(a(i, k, 1)-a(i, kp, 0))
        t3id = sin72*(ad(i, kp, 3)+ad(i, k, 4))
        t3i = sin72*(a(i, kp, 3)+a(i, k, 4))
        t4rd = sin72*(ad(i, k, 2)-ad(i, kp, 1))
        t4r = sin72*(a(i, k, 2)-a(i, kp, 1))
        t4id = sin72*(ad(i, kp, 2)+ad(i, k, 3))
        t4i = sin72*(a(i, kp, 2)+a(i, k, 3))
        t5rd = t1rd + t2rd
        t5r = t1r + t2r
        t5id = t1id + t2id
        t5i = t1i + t2i
        t6rd = r54*(t1rd-t2rd)
        t6r = r54*(t1r-t2r)
        t6id = r54*(t1id-t2id)
        t6i = r54*(t1i-t2i)
        t7rd = ad(i, k, 0) - 0.25d0*t5rd
        t7r = a(i, k, 0) - 0.25d0*t5r
        t7id = ad(i, kp, 4) - 0.25d0*t5id
        t7i = a(i, kp, 4) - 0.25d0*t5i
        t8rd = t7rd + t6rd
        t8r = t7r + t6r
        t8id = t7id + t6id
        t8i = t7i + t6i
        t9rd = t7rd - t6rd
        t9r = t7r - t6r
        t9id = t7id - t6id
        t9i = t7i - t6i
        t10rd = t3rd + s36s72*t4rd
        t10r = t3r + s36s72*t4r
        t10id = t3id + s36s72*t4id
        t10i = t3i + s36s72*t4i
        t11rd = s36s72*t3rd - t4rd
        t11r = s36s72*t3r - t4r
        t11id = s36s72*t3id - t4id
        t11i = s36s72*t3i - t4i
        x0pd = ad(i, k, 0) + t5rd
        x0p = a(i, k, 0) + t5r
        y0pd = ad(i, kp, 4) + t5id
        y0p = a(i, kp, 4) + t5i
        x1pd = t8rd + t10id
        x1p = t8r + t10i
        y1pd = t8id - t10rd
        y1p = t8i - t10r
        x2pd = t9rd + t11id
        x2p = t9r + t11i
        y2pd = t9id - t11rd
        y2p = t9i - t11r
        x3pd = t9rd - t11id
        x3p = t9r - t11i
        y3pd = t9id + t11rd
        y3p = t9i + t11r
        x4pd = t8rd - t10id
        x4p = t8r - t10i
        y4pd = t8id + t10rd
        y4p = t8i + t10r
        bd(i, 0, k) = x0pd
        b(i, 0, k) = x0p
        bd(i, 0, kp) = y0pd
        b(i, 0, kp) = y0p
        bd(i, 1, k) = ck*x1pd - sk*y1pd
        b(i, 1, k) = ck*x1p - sk*y1p
        bd(i, 1, kp) = ck*y1pd + sk*x1pd
        b(i, 1, kp) = ck*y1p + sk*x1p
        bd(i, 2, k) = c2k*x2pd - s2k*y2pd
        b(i, 2, k) = c2k*x2p - s2k*y2p
        bd(i, 2, kp) = c2k*y2pd + s2k*x2pd
        b(i, 2, kp) = c2k*y2p + s2k*x2p
        bd(i, 3, k) = c3k*x3pd - s3k*y3pd
        b(i, 3, k) = c3k*x3p - s3k*y3p
        bd(i, 3, kp) = c3k*y3pd + s3k*x3pd
        b(i, 3, kp) = c3k*y3p + s3k*x3p
        bd(i, 4, k) = c4k*x4pd - s4k*y4pd
        b(i, 4, k) = c4k*x4p - s4k*y4p
        bd(i, 4, kp) = c4k*y4pd + s4k*x4pd
        b(i, 4, kp) = c4k*y4p + s4k*x4p
      END DO
    END DO
  END IF
!
!     Code for K=R/2 when R is even not needed
!
!     IF (MOD(R,2).EQ.0) THEN
!        R2 = R/2
!        DO 120 I = 0, P - 1
!           T1 = A(I,R2,0) + A(I,R2,1)
!           T2 = 0.25D0*T1 - A(I,R2,2)
!           T3 = R54*(A(I,R2,0)-A(I,R2,1))
!           T4 = SIN36*A(I,R2,4) + SIN72*A(I,R2,3)
!           T5 = SIN72*A(I,R2,4) - SIN36*A(I,R2,3)
!           T6 = T2 + T3
!           T7 = T2 - T3
!           B(I,0,R2) = T1 + A(I,R2,2)
!           B(I,1,R2) = T4 + T6
!           B(I,2,R2) = T5 - T7
!           B(I,3,R2) = T5 + T7
!           B(I,4,R2) = T4 - T6
! 120    CONTINUE
!     END IF
!
  RETURN
END SUBROUTINE C06FQT_D