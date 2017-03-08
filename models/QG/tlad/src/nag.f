      SUBROUTINE C06FPF(M,N,X,INIT,TRIG,WORK,IFAIL)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FPF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N), WORK(M*N), X(M*N)
C     .. Local Scalars ..
      INTEGER           IERROR, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FPX
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0) THEN
         CALL C06FPX(X,WORK,M,N,Q,NQ,TRIG)
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,', but TRIG array never initialized')
99995 FORMAT (' ** INIT = ',A1,', but N and TRIG array incompatible')
      END
      SUBROUTINE C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           IERROR, M, N, NQ
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N)
      INTEGER           Q(30)
C     .. Local Scalars ..
      INTEGER           NCHECK
C     .. External Subroutines ..
      EXTERNAL          C06FPY, C06FPZ
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     .. Save statement ..
      SAVE              NCHECK
C     .. Data statements ..
      DATA              NCHECK/-1/
C     .. Executable Statements ..
      IERROR = 0
C
      IF (M.LT.1) THEN
         IERROR = 1
         RETURN
      ELSE IF (N.LT.1) THEN
         IERROR = 2
         RETURN
      END IF
      IF (INIT.NE.'I' .AND. INIT.NE.'i' .AND. INIT.NE.'S' .AND. INIT.NE.
     *    's' .AND. INIT.NE.'R' .AND. INIT.NE.'r') THEN
         IERROR = 3
         RETURN
      END IF
      IF (INIT.EQ.'S' .OR. INIT.EQ.'s') THEN
         IF (NCHECK.EQ.-1) THEN
            IERROR = 4
            RETURN
         ELSE IF (NINT(TRIG(N)).NE.N .OR. NINT(TRIG(2*N)).NE.N) THEN
            IERROR = 5
            RETURN
         END IF
      END IF
      IF (INIT.EQ.'R' .OR. INIT.EQ.'r') THEN
         IF (NINT(TRIG(N)).NE.N .OR. NINT(TRIG(2*N)).NE.N) THEN
            IERROR = 5
            RETURN
         END IF
      END IF
      CALL C06FPZ(N,NQ,Q)
      IF (INIT.EQ.'I' .OR. INIT.EQ.'i') THEN
         CALL C06FPY(N,NQ,Q,TRIG(1),TRIG(N+1))
      END IF
      NCHECK = N
      RETURN
      END
      SUBROUTINE C06FPR(A,B,P,Q,R,COSINE,SINE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Real to Hermitian fast Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:Q-1,0:R-1), B(0:P-1,0:R-1,0:Q-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, AR, TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER           I, INDX, J, K, KP, L, Q2, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      Q2 = (Q-1)/2
      IF (P.GE.R/2) THEN
C
C        Code for K=0 --
C
         DO 40 J = 1, Q2
            DO 20 I = 0, P - 1
               TEMP1 = A(I,J,0)
               A(I,J,0) = TEMP1 + A(I,Q-J,0)
               A(I,Q-J,0) = TEMP1 - A(I,Q-J,0)
   20       CONTINUE
   40    CONTINUE
         DO 120 L = 1, Q2
            DO 60 I = 0, P - 1
               B(I,0,L) = A(I,0,0)
               B(I,0,Q-L) = 0.0D0
   60       CONTINUE
            DO 100 J = 1, Q2
               INDX = MOD(J*L,Q)
               DO 80 I = 0, P - 1
                  B(I,0,L) = B(I,0,L) + A(I,J,0)*COSINE(0,INDX)
                  B(I,0,Q-L) = B(I,0,Q-L) + A(I,Q-J,0)*SINE(0,INDX)
   80          CONTINUE
  100       CONTINUE
  120    CONTINUE
         DO 140 I = 0, P - 1
            B(I,0,0) = A(I,0,0)
  140    CONTINUE
         DO 180 J = 1, Q2
            DO 160 I = 0, P - 1
               B(I,0,0) = B(I,0,0) + A(I,J,0)
  160       CONTINUE
  180    CONTINUE
C
C        Code for general K --
C
         DO 460 K = 1, (R-1)/2
            KP = R - K
            DO 220 J = 1, Q - 1
               DO 200 I = 0, P - 1
                  AR = A(I,J,K)
                  AI = A(I,J,KP)
                  A(I,J,K) = COSINE(K,J)*AR - SINE(K,J)*AI
                  A(I,J,KP) = COSINE(K,J)*AI + SINE(K,J)*AR
  200          CONTINUE
  220       CONTINUE
            DO 260 J = 1, Q2
               DO 240 I = 0, P - 1
                  TEMP1 = A(I,J,K)
                  TEMP2 = A(I,J,KP)
                  A(I,J,K) = TEMP1 + A(I,Q-J,K)
                  A(I,J,KP) = TEMP2 + A(I,Q-J,KP)
                  A(I,Q-J,K) = TEMP1 - A(I,Q-J,K)
                  A(I,Q-J,KP) = TEMP2 - A(I,Q-J,KP)
  240          CONTINUE
  260       CONTINUE
            DO 340 L = 1, Q2
               DO 280 I = 0, P - 1
                  B(I,K,L) = A(I,0,K)
                  B(I,KP,Q-L-1) = A(I,0,KP)
                  B(I,KP,L-1) = 0.0D0
                  B(I,K,Q-L) = 0.0D0
  280          CONTINUE
               DO 320 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 300 I = 0, P - 1
                     B(I,K,L) = B(I,K,L) + A(I,J,K)*COSINE(0,INDX)
                     B(I,KP,Q-L-1) = B(I,KP,Q-L-1) + A(I,J,KP)*COSINE(0,
     *                               INDX)
                     B(I,KP,L-1) = B(I,KP,L-1) - A(I,Q-J,K)*SINE(0,INDX)
                     B(I,K,Q-L) = B(I,K,Q-L) + A(I,Q-J,KP)*SINE(0,INDX)
  300             CONTINUE
  320          CONTINUE
  340       CONTINUE
            DO 360 I = 0, P - 1
               B(I,K,0) = A(I,0,K)
               B(I,KP,Q-1) = A(I,0,KP)
  360       CONTINUE
            DO 400 J = 1, Q2
               DO 380 I = 0, P - 1
                  B(I,K,0) = B(I,K,0) + A(I,J,K)
                  B(I,KP,Q-1) = B(I,KP,Q-1) + A(I,J,KP)
  380          CONTINUE
  400       CONTINUE
            DO 440 L = 1, Q2
               DO 420 I = 0, P - 1
                  TEMP1 = B(I,K,L)
                  TEMP2 = B(I,KP,Q-L-1)
                  TEMP3 = B(I,KP,L-1)
                  B(I,K,L) = TEMP1 - B(I,K,Q-L)
                  B(I,KP,L-1) = TEMP1 + B(I,K,Q-L)
                  B(I,KP,Q-L-1) = TEMP2 - TEMP3
                  B(I,K,Q-L) = -TEMP2 - TEMP3
  420          CONTINUE
  440       CONTINUE
  460    CONTINUE
C
C        Code for K=R/2 if R is even --
C
         IF (MOD(R,2).EQ.0) THEN
            R2 = R/2
            DO 500 J = 1, Q2 - 1, 2
               DO 480 I = 0, P - 1
                  TEMP1 = A(I,J,R2)
                  TEMP2 = A(I,Q-J,R2)
                  A(I,J,R2) = -TEMP1 + TEMP2
                  A(I,Q-J,R2) = -TEMP1 - TEMP2
                  TEMP3 = A(I,J+1,R2)
                  TEMP4 = A(I,Q-J-1,R2)
                  A(I,J+1,R2) = TEMP3 - TEMP4
                  A(I,Q-J-1,R2) = TEMP3 + TEMP4
  480          CONTINUE
  500       CONTINUE
            IF (MOD(Q2,2).EQ.1) THEN
               DO 520 I = 0, P - 1
                  TEMP1 = A(I,Q2,R2)
                  TEMP2 = A(I,Q2+1,R2)
                  A(I,Q2,R2) = -TEMP1 + TEMP2
                  A(I,Q2+1,R2) = -TEMP1 - TEMP2
  520          CONTINUE
            END IF
            DO 600 L = 0, Q2 - 1
               DO 540 I = 0, P - 1
                  B(I,R2,L) = A(I,0,R2)
                  B(I,R2,Q-L-1) = 0.0D0
  540          CONTINUE
               DO 580 J = 1, Q2
                  INDX = MOD(J*(L+Q2+1),Q)
                  DO 560 I = 0, P - 1
                     B(I,R2,L) = B(I,R2,L) + A(I,J,R2)*COSINE(0,INDX)
                     B(I,R2,Q-L-1) = B(I,R2,Q-L-1) + A(I,Q-J,R2)*SINE(0,
     *                               INDX)
  560             CONTINUE
  580          CONTINUE
  600       CONTINUE
            DO 620 I = 0, P - 1
               B(I,R2,Q2) = A(I,0,R2)
  620       CONTINUE
            DO 660 J = 1, Q2
               DO 640 I = 0, P - 1
                  B(I,R2,Q2) = B(I,R2,Q2) + A(I,J,R2)
  640          CONTINUE
  660       CONTINUE
         END IF
C
      ELSE
C
         DO 1100 I = 0, P - 1
C
C           Code for K=0 --
C
            DO 680 J = 1, Q2
               TEMP1 = A(I,J,0)
               A(I,J,0) = TEMP1 + A(I,Q-J,0)
               A(I,Q-J,0) = TEMP1 - A(I,Q-J,0)
  680       CONTINUE
            DO 720 L = 1, Q2
               B(I,0,L) = A(I,0,0)
               B(I,0,Q-L) = 0.0D0
               DO 700 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  B(I,0,L) = B(I,0,L) + A(I,J,0)*COSINE(0,INDX)
                  B(I,0,Q-L) = B(I,0,Q-L) + A(I,Q-J,0)*SINE(0,INDX)
  700          CONTINUE
  720       CONTINUE
            B(I,0,0) = A(I,0,0)
            DO 740 J = 1, Q2
               B(I,0,0) = B(I,0,0) + A(I,J,0)
  740       CONTINUE
C
C           Code for general K --
C
            DO 780 J = 1, Q - 1
               DO 760 K = 1, (R-1)/2
                  AR = A(I,J,K)
                  AI = A(I,J,R-K)
                  A(I,J,K) = COSINE(K,J)*AR - SINE(K,J)*AI
                  A(I,J,R-K) = COSINE(K,J)*AI + SINE(K,J)*AR
  760          CONTINUE
  780       CONTINUE
            DO 820 J = 1, Q2
               DO 800 K = 1, (R-1)/2
                  TEMP1 = A(I,J,K)
                  TEMP2 = A(I,J,R-K)
                  A(I,J,K) = TEMP1 + A(I,Q-J,K)
                  A(I,J,R-K) = TEMP2 + A(I,Q-J,R-K)
                  A(I,Q-J,K) = TEMP1 - A(I,Q-J,K)
                  A(I,Q-J,R-K) = TEMP2 - A(I,Q-J,R-K)
  800          CONTINUE
  820       CONTINUE
            DO 900 L = 1, Q2
               DO 840 K = 1, (R-1)/2
                  B(I,K,L) = A(I,0,K)
                  B(I,R-K,Q-L-1) = A(I,0,R-K)
                  B(I,R-K,L-1) = 0.0D0
                  B(I,K,Q-L) = 0.0D0
  840          CONTINUE
               DO 880 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 860 K = 1, (R-1)/2
                     B(I,K,L) = B(I,K,L) + A(I,J,K)*COSINE(0,INDX)
                     B(I,R-K,Q-L-1) = B(I,R-K,Q-L-1) + A(I,J,R-K)
     *                                *COSINE(0,INDX)
                     B(I,R-K,L-1) = B(I,R-K,L-1) - A(I,Q-J,K)*SINE(0,
     *                              INDX)
                     B(I,K,Q-L) = B(I,K,Q-L) + A(I,Q-J,R-K)*SINE(0,INDX)
  860             CONTINUE
  880          CONTINUE
  900       CONTINUE
            DO 920 K = 1, (R-1)/2
               B(I,K,0) = A(I,0,K)
               B(I,R-K,Q-1) = A(I,0,R-K)
  920       CONTINUE
            DO 960 J = 1, Q2
               DO 940 K = 1, (R-1)/2
                  B(I,K,0) = B(I,K,0) + A(I,J,K)
                  B(I,R-K,Q-1) = B(I,R-K,Q-1) + A(I,J,R-K)
  940          CONTINUE
  960       CONTINUE
            DO 1000 L = 1, Q2
               DO 980 K = 1, (R-1)/2
                  TEMP1 = B(I,K,L)
                  TEMP2 = B(I,R-K,Q-L-1)
                  TEMP3 = B(I,R-K,L-1)
                  B(I,K,L) = TEMP1 - B(I,K,Q-L)
                  B(I,R-K,L-1) = TEMP1 + B(I,K,Q-L)
                  B(I,R-K,Q-L-1) = TEMP2 - TEMP3
                  B(I,K,Q-L) = -TEMP2 - TEMP3
  980          CONTINUE
 1000       CONTINUE
C
C           Code for K=R/2 if R is even --
C
            IF (MOD(R,2).EQ.0) THEN
               R2 = R/2
               DO 1020 J = 1, Q2 - 1, 2
                  TEMP1 = A(I,J,R2)
                  TEMP2 = A(I,Q-J,R2)
                  A(I,J,R2) = -TEMP1 + TEMP2
                  A(I,Q-J,R2) = -TEMP1 - TEMP2
                  TEMP3 = A(I,J+1,R2)
                  TEMP4 = A(I,Q-J-1,R2)
                  A(I,J+1,R2) = TEMP3 - TEMP4
                  A(I,Q-J-1,R2) = TEMP3 + TEMP4
 1020          CONTINUE
               IF (MOD(Q2,2).EQ.1) THEN
                  TEMP1 = A(I,Q2,R2)
                  TEMP2 = A(I,Q2+1,R2)
                  A(I,Q2,R2) = -TEMP1 + TEMP2
                  A(I,Q2+1,R2) = -TEMP1 - TEMP2
               END IF
               DO 1060 L = 0, Q2 - 1
                  B(I,R2,L) = A(I,0,R2)
                  B(I,R2,Q-L-1) = 0.0D0
                  DO 1040 J = 1, Q2
                     INDX = MOD(J*(L+Q2+1),Q)
                     B(I,R2,L) = B(I,R2,L) + A(I,J,R2)*COSINE(0,INDX)
                     B(I,R2,Q-L-1) = B(I,R2,Q-L-1) + A(I,Q-J,R2)*SINE(0,
     *                               INDX)
 1040             CONTINUE
 1060          CONTINUE
               B(I,R2,Q2) = A(I,0,R2)
               DO 1080 J = 1, Q2
                  B(I,R2,Q2) = B(I,R2,Q2) + A(I,J,R2)
 1080          CONTINUE
            END IF
 1100    CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE C06FPS(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-694 (DEC 1989).
C
C     Radix six real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:5,0:R-1), B(0:P-1,0:R-1,0:5),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1, T1I, T1R, T2, T2I, T2R, T3, T3I, T3R, T4,
     *                  T5, T6, U0, U0I, U0R, U1I, U1R, U2I, U2R, UI,
     *                  UR, V0, V0I, V0R, V1I, V1R, V2I, V2R, VI, VR,
     *                  X1P, X2P, X3P, X4P, X5P, Y1P, Y2P, Y3P, Y4P, Y5P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,2,0) + A(I,4,0)
         UR = A(I,0,0) - 0.5D0*T1
         UI = -SIN60*(A(I,2,0)-A(I,4,0))
         U0 = A(I,0,0) + T1
         T1 = A(I,5,0) + A(I,1,0)
         VR = A(I,3,0) - 0.5D0*T1
         VI = -SIN60*(A(I,5,0)-A(I,1,0))
         V0 = A(I,3,0) + T1
         B(I,0,0) = U0 + V0
         B(I,0,1) = UR - VR
         B(I,0,2) = UR + VR
         B(I,0,3) = U0 - V0
         B(I,0,4) = -UI - VI
         B(I,0,5) = UI - VI
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               X4P = COSINE(K,4)*A(I,4,K) - SINE(K,4)*A(I,4,R-K)
               Y4P = COSINE(K,4)*A(I,4,R-K) + SINE(K,4)*A(I,4,K)
               X5P = COSINE(K,5)*A(I,5,K) - SINE(K,5)*A(I,5,R-K)
               Y5P = COSINE(K,5)*A(I,5,R-K) + SINE(K,5)*A(I,5,K)
               T1R = X2P + X4P
               T1I = Y2P + Y4P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,R-K) - 0.5D0*T1I
               T3R = SIN60*(X2P-X4P)
               T3I = SIN60*(Y2P-Y4P)
               U0R = A(I,0,K) + T1R
               U0I = A(I,0,R-K) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X5P + X1P
               T1I = Y5P + Y1P
               T2R = X3P - 0.5D0*T1R
               T2I = Y3P - 0.5D0*T1I
               T3R = SIN60*(X5P-X1P)
               T3I = SIN60*(Y5P-Y1P)
               V0R = X3P + T1R
               V0I = Y3P + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               B(I,K,0) = U0R + V0R
               B(I,R-K,0) = U2R - V2R
               B(I,K,1) = U1R - V1R
               B(I,R-K,1) = U1R + V1R
               B(I,K,2) = U2R + V2R
               B(I,R-K,2) = U0R - V0R
               B(I,K,3) = -U0I + V0I
               B(I,R-K,3) = U2I + V2I
               B(I,K,4) = -U1I - V1I
               B(I,R-K,4) = U1I - V1I
               B(I,K,5) = -U2I + V2I
               B(I,R-K,5) = U0I + V0I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            C5K = COSINE(K,5)
            S5K = SINE(K,5)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               X4P = C4K*A(I,4,K) - S4K*A(I,4,KP)
               Y4P = C4K*A(I,4,KP) + S4K*A(I,4,K)
               X5P = C5K*A(I,5,K) - S5K*A(I,5,KP)
               Y5P = C5K*A(I,5,KP) + S5K*A(I,5,K)
               T1R = X2P + X4P
               T1I = Y2P + Y4P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,KP) - 0.5D0*T1I
               T3R = SIN60*(X2P-X4P)
               T3I = SIN60*(Y2P-Y4P)
               U0R = A(I,0,K) + T1R
               U0I = A(I,0,KP) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X5P + X1P
               T1I = Y5P + Y1P
               T2R = X3P - 0.5D0*T1R
               T2I = Y3P - 0.5D0*T1I
               T3R = SIN60*(X5P-X1P)
               T3I = SIN60*(Y5P-Y1P)
               V0R = X3P + T1R
               V0I = Y3P + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               B(I,K,0) = U0R + V0R
               B(I,KP,0) = U2R - V2R
               B(I,K,1) = U1R - V1R
               B(I,KP,1) = U1R + V1R
               B(I,K,2) = U2R + V2R
               B(I,KP,2) = U0R - V0R
               B(I,K,3) = -U0I + V0I
               B(I,KP,3) = U2I + V2I
               B(I,K,4) = -U1I - V1I
               B(I,KP,4) = U1I - V1I
               B(I,K,5) = -U2I + V2I
               B(I,KP,5) = U0I + V0I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = A(I,2,R2) - A(I,4,R2)
            T2 = A(I,0,R2) + 0.5D0*T1
            T3 = SIN60*(A(I,2,R2)+A(I,4,R2))
            T4 = A(I,1,R2) + A(I,5,R2)
            T5 = -A(I,3,R2) - 0.5D0*T4
            T6 = SIN60*(A(I,1,R2)-A(I,5,R2))
            B(I,R2,0) = T2 + T6
            B(I,R2,1) = A(I,0,R2) - T1
            B(I,R2,2) = T2 - T6
            B(I,R2,3) = T5 + T3
            B(I,R2,4) = A(I,3,R2) - T4
            B(I,R2,5) = T5 - T3
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FPT(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-695 (DEC 1989).
C
C     Radix five real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  R54, SIN72, S36S72
      PARAMETER         (R54=0.559016994374947424102293417182819D0,
     *                  SIN72=0.951056516295153572116439333379382D0,
     *                  S36S72=0.618033988749894848204586834365638D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:4,0:R-1), B(0:P-1,0:R-1,0:4),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T1, T10I,
     *                  T10R, T11I, T11R, T1I, T1R, T2, T2I, T2R, T3,
     *                  T3I, T3R, T4, T4I, T4R, T5, T5I, T5R, T6, T6I,
     *                  T6R, T7, T7I, T7R, T8I, T8R, T9I, T9R, X1P, X2P,
     *                  X3P, X4P, Y1P, Y2P, Y3P, Y4P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,1,0) + A(I,4,0)
         T2 = A(I,2,0) + A(I,3,0)
         T3 = SIN72*(A(I,1,0)-A(I,4,0))
         T4 = SIN72*(A(I,2,0)-A(I,3,0))
         T5 = T1 + T2
         T6 = R54*(T1-T2)
         T7 = A(I,0,0) - 0.25D0*T5
         B(I,0,0) = A(I,0,0) + T5
         B(I,0,1) = T7 + T6
         B(I,0,2) = T7 - T6
         B(I,0,3) = -S36S72*T3 + T4
         B(I,0,4) = -T3 - S36S72*T4
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               X4P = COSINE(K,4)*A(I,4,K) - SINE(K,4)*A(I,4,R-K)
               Y4P = COSINE(K,4)*A(I,4,R-K) + SINE(K,4)*A(I,4,K)
               T1R = X1P + X4P
               T1I = Y1P + Y4P
               T2R = X2P + X3P
               T2I = Y2P + Y3P
               T3R = SIN72*(X1P-X4P)
               T3I = SIN72*(Y1P-Y4P)
               T4R = SIN72*(X2P-X3P)
               T4I = SIN72*(Y2P-Y3P)
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,0,K) - 0.25D0*T5R
               T7I = A(I,0,R-K) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               B(I,K,0) = A(I,0,K) + T5R
               B(I,R-K,0) = T8R - T10I
               B(I,K,1) = T8R + T10I
               B(I,R-K,1) = T9R - T11I
               B(I,K,2) = T9R + T11I
               B(I,R-K,2) = T9I - T11R
               B(I,K,3) = -T9I - T11R
               B(I,R-K,3) = T8I - T10R
               B(I,K,4) = -T8I - T10R
               B(I,R-K,4) = A(I,0,R-K) + T5I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               X4P = C4K*A(I,4,K) - S4K*A(I,4,KP)
               Y4P = C4K*A(I,4,KP) + S4K*A(I,4,K)
               T1R = X1P + X4P
               T1I = Y1P + Y4P
               T2R = X2P + X3P
               T2I = Y2P + Y3P
               T3R = SIN72*(X1P-X4P)
               T3I = SIN72*(Y1P-Y4P)
               T4R = SIN72*(X2P-X3P)
               T4I = SIN72*(Y2P-Y3P)
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,0,K) - 0.25D0*T5R
               T7I = A(I,0,KP) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               B(I,K,0) = A(I,0,K) + T5R
               B(I,KP,0) = T8R - T10I
               B(I,K,1) = T8R + T10I
               B(I,KP,1) = T9R - T11I
               B(I,K,2) = T9R + T11I
               B(I,KP,2) = T9I - T11R
               B(I,K,3) = -T9I - T11R
               B(I,KP,3) = T8I - T10R
               B(I,K,4) = -T8I - T10R
               B(I,KP,4) = A(I,0,KP) + T5I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = SIN72*(A(I,1,R2)+A(I,4,R2))
C           T2 = SIN72*(A(I,2,R2)+A(I,3,R2))
C           T3 = A(I,1,R2) - A(I,4,R2)
C           T4 = A(I,2,R2) - A(I,3,R2)
C           T5 = T4 - T3
C           T6 = R54*(T4+T3)
C           T7 = A(I,0,R2) - 0.25D0*T5
C           B(I,R2,0) = T7 + T6
C           B(I,R2,1) = T7 - T6
C           B(I,R2,2) = A(I,0,R2) + T5
C           B(I,R2,3) = -T1 + S36S72*T2
C           B(I,R2,4) = -S36S72*T1 - T2
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FPU(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-696 (DEC 1989).
C
C     Radix four real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  ROOT2I
      PARAMETER         (ROOT2I=0.707106781186547524400844362104849D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:3,0:R-1), B(0:P-1,0:R-1,0:3),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, CK, S2K, S3K, SK, T1, T1I, T1R, T2,
     *                  T2I, T2R, T3I, T3R, T4I, T4R, X1P, X2P, X3P,
     *                  Y1P, Y2P, Y3P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,0) + A(I,2,0)
         T2 = A(I,1,0) + A(I,3,0)
         B(I,0,0) = T1 + T2
         B(I,0,1) = A(I,0,0) - A(I,2,0)
         B(I,0,2) = T1 - T2
         B(I,0,3) = -A(I,1,0) + A(I,3,0)
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               T1R = A(I,0,K) + X2P
               T1I = A(I,0,R-K) + Y2P
               T2R = X1P + X3P
               T2I = Y1P + Y3P
               T3R = A(I,0,K) - X2P
               T3I = A(I,0,R-K) - Y2P
               T4R = X1P - X3P
               T4I = Y1P - Y3P
               B(I,K,0) = T1R + T2R
               B(I,R-K,0) = T3R - T4I
               B(I,K,1) = T3R + T4I
               B(I,R-K,1) = T1R - T2R
               B(I,K,2) = T2I - T1I
               B(I,R-K,2) = T3I - T4R
               B(I,K,3) = -T3I - T4R
               B(I,R-K,3) = T1I + T2I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               T1R = A(I,0,K) + X2P
               T1I = A(I,0,KP) + Y2P
               T2R = X1P + X3P
               T2I = Y1P + Y3P
               T3R = A(I,0,K) - X2P
               T3I = A(I,0,KP) - Y2P
               T4R = X1P - X3P
               T4I = Y1P - Y3P
               B(I,K,0) = T1R + T2R
               B(I,KP,0) = T3R - T4I
               B(I,K,1) = T3R + T4I
               B(I,KP,1) = T1R - T2R
               B(I,K,2) = T2I - T1I
               B(I,KP,2) = T3I - T4R
               B(I,K,3) = -T3I - T4R
               B(I,KP,3) = T1I + T2I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = ROOT2I*(A(I,1,R2)-A(I,3,R2))
            T2 = ROOT2I*(A(I,1,R2)+A(I,3,R2))
            B(I,R2,0) = A(I,0,R2) + T1
            B(I,R2,1) = A(I,0,R2) - T1
            B(I,R2,2) = A(I,2,R2) - T2
            B(I,R2,3) = -A(I,2,R2) - T2
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FPV(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-697 (DEC 1989).
C
C     Radix three Real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:2,0:R-1), B(0:P-1,0:R-1,0:2),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, CK, S2K, SK, T1, T1I, T1R, T2I, T2R, T3I,
     *                  T3R, X1P, X2P, Y1P, Y2P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,1,0) + A(I,2,0)
         B(I,0,0) = A(I,0,0) + T1
         B(I,0,1) = A(I,0,0) - 0.5D0*T1
         B(I,0,2) = -SIN60*(A(I,1,0)-A(I,2,0))
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               T1R = X1P + X2P
               T1I = Y1P + Y2P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,R-K) - 0.5D0*T1I
               T3R = SIN60*(X1P-X2P)
               T3I = SIN60*(Y1P-Y2P)
               B(I,K,0) = A(I,0,K) + T1R
               B(I,R-K,0) = T2R - T3I
               B(I,K,1) = T2R + T3I
               B(I,R-K,1) = T2I - T3R
               B(I,K,2) = -(T2I+T3R)
               B(I,R-K,2) = A(I,0,R-K) + T1I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               T1R = X1P + X2P
               T1I = Y1P + Y2P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,KP) - 0.5D0*T1I
               T3R = SIN60*(X1P-X2P)
               T3I = SIN60*(Y1P-Y2P)
               B(I,K,0) = A(I,0,K) + T1R
               B(I,KP,0) = T2R - T3I
               B(I,K,1) = T2R + T3I
               B(I,KP,1) = T2I - T3R
               B(I,K,2) = -(T2I+T3R)
               B(I,KP,2) = A(I,0,KP) + T1I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,1,R2) - A(I,2,R2)
C           B(I,R2,0) = A(I,0,R2) + 0.5D0*T1
C           B(I,R2,1) = A(I,0,R2) - T1
C           B(I,R2,2) = -SIN60*(A(I,1,R2)+A(I,2,R2))
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FPW(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-698 (DEC 1989).
C
C     Radix two Real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:1,0:R-1), B(0:P-1,0:R-1,0:1),
     *                  COSINE(0:R-1), SINE(0:R-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  CK, SK, X1, Y1
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         B(I,0,0) = A(I,0,0) + A(I,1,0)
         B(I,0,1) = A(I,0,0) - A(I,1,0)
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1 = COSINE(K)*A(I,1,K) - SINE(K)*A(I,1,R-K)
               Y1 = COSINE(K)*A(I,1,R-K) + SINE(K)*A(I,1,K)
               B(I,K,0) = A(I,0,K) + X1
               B(I,R-K,0) = A(I,0,K) - X1
               B(I,K,1) = -A(I,0,R-K) + Y1
               B(I,R-K,1) = A(I,0,R-K) + Y1
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K)
            SK = SINE(K)
            DO 80 I = 0, P - 1
               X1 = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1 = CK*A(I,1,KP) + SK*A(I,1,K)
               B(I,K,0) = A(I,0,K) + X1
               B(I,KP,0) = A(I,0,K) - X1
               B(I,K,1) = -A(I,0,KP) + Y1
               B(I,KP,1) = A(I,0,KP) + Y1
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           B(I,R2,0) = A(I,0,R2)
C           B(I,R2,1) = -A(I,1,R2)
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FPX(A,B,M,N,Q,NQ,TRIG)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Real to Hermitian Fast Fourier Transform Kernel Driver
C
C     Mixed-radix, self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:M*N-1), B(0:M*N-1), TRIG(0:2*N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR
      INTEGER           I, P, QI, R
      LOGICAL           INA
C     .. External Subroutines ..
      EXTERNAL          C06FPR, C06FPS, C06FPT, C06FPU, C06FPV, C06FPW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = N
      R = 1
      IF (N.EQ.1) RETURN
      DO 20 I = NQ, 1, -1
         QI = Q(I)
         P = P/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FPW(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FPV(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FPU(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FPT(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FPS(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FPR(A,B,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FPW(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FPV(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FPU(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FPT(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FPS(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FPR(B,A,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         END IF
         INA = .NOT. INA
         R = R*QI
   20 CONTINUE
C
      FACTOR = 1.0D0/SQRT(DBLE(N))
      IF (INA) THEN
         DO 40 I = 0, M*N - 1
            A(I) = A(I)*FACTOR
   40    CONTINUE
      ELSE
         DO 60 I = 0, M*N - 1
            A(I) = B(I)*FACTOR
   60    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FPY(N,NQ,Q,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Trig function initialisation subroutine
C
C     .. Scalar Arguments ..
      INTEGER           N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  COSINE(0:N-1), SINE(0:N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  TWOPI, Z
      INTEGER           I, J, K, L, L1, QI, R
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, SIN, DBLE
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      Z = TWOPI/DBLE(N)
C
      R = N
      L = 0
C
      DO 80 I = 1, NQ
         QI = Q(I)
         R = R/QI
         L1 = L
         DO 40 J = 1, QI - 1
            DO 20 K = 0, R - 1
               COSINE(L) = Z*DBLE(J*K)
               L = L + 1
   20       CONTINUE
   40    CONTINUE
         IF (QI.GE.7) THEN
            L = L1
            DO 60 J = 1, QI - 1
               COSINE(L) = Z*J*R
               L = L + R
   60       CONTINUE
         END IF
         Z = Z*QI
   80 CONTINUE
C
      DO 100 I = 0, N - 2
         SINE(I) = -SIN(COSINE(I))
         COSINE(I) = COS(COSINE(I))
  100 CONTINUE
C
C     Check on consistency of N and TRIG array --
C
      COSINE(N-1) = DBLE(N)
      SINE(N-1) = DBLE(N)
C
      RETURN
      END
      SUBROUTINE C06FPZ(N,NQ,Q)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           N, NQ
C     .. Array Arguments ..
      INTEGER           Q(30)
C     .. Local Scalars ..
      INTEGER           I, K, L, NN
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      NN = N
      K = 0
C
C     Trap the special case N = 1
C
      IF (N.EQ.1) THEN
         NQ = 1
         Q(1) = 1
         RETURN
      END IF
C
C     Factors of 6 --
C
   20 IF (MOD(NN,6).NE.0) GO TO 40
      K = K + 1
      Q(K) = 6
      NN = NN/6
      IF (NN.EQ.1) GO TO 160
      GO TO 20
C
C     Factors of 4 --
C
   40 IF (MOD(NN,4).NE.0) GO TO 60
      K = K + 1
      Q(K) = 4
      NN = NN/4
      IF (NN.EQ.1) GO TO 160
      GO TO 40
C
C     Factors of 2 --
C
   60 IF (MOD(NN,2).NE.0) GO TO 80
      K = K + 1
      Q(K) = 2
      NN = NN/2
      IF (NN.EQ.1) GO TO 160
      GO TO 60
C
C     Factors of 3 --
C
   80 IF (MOD(NN,3).NE.0) GO TO 100
      K = K + 1
      Q(K) = 3
      NN = NN/3
      IF (NN.EQ.1) GO TO 160
      GO TO 80
C
C     Remaining odd factors --
C
  100 L = 5
      I = 2
C
C     I is alternatively 2 or 4 --
C
  120 IF (MOD(NN,L).NE.0) GO TO 140
      K = K + 1
      Q(K) = L
      NN = NN/L
      IF (NN.EQ.1) GO TO 160
      GO TO 120
  140 L = L + I
      I = 6 - I
      GO TO 120
  160 NQ = K
C
      RETURN
      END
      SUBROUTINE C06FQF(M,N,X,INIT,TRIG,WORK,IFAIL)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FQF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N), WORK(M*N), X(M*N)
C     .. Local Scalars ..
      INTEGER           IERROR, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FQX
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0) THEN
         CALL C06FQX(X,WORK,M,N,Q,NQ,TRIG)
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,', but TRIG array never initialized')
99995 FORMAT (' ** INIT = ',A1,', but N and TRIG array incompatible')
      END
      SUBROUTINE C06FQQ(A,M,N)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:M-1,0:N-1)
C     .. Local Scalars ..
      INTEGER           L
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 L = 0, M - 1
         A(L,0) = 0.5D0*A(L,0)
   20 CONTINUE
      IF (MOD(N,2).EQ.0) THEN
         DO 40 L = 0, M - 1
            A(L,N/2) = 0.5D0*A(L,N/2)
   40    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FQR(A,B,P,Q,R,COSINE,SINE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Hermitian to real fast Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:Q-1), B(0:P-1,0:Q-1,0:R-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  BI, BR, TEMP, TEMP1, TEMPI, TEMPR
      INTEGER           I, INDX, J, K, KP, L, Q2, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      Q2 = (Q-1)/2
      IF (P.GE.R/2) THEN
C
C        Code for K=0 --
C
         DO 80 L = 1, Q2
            DO 20 I = 0, P - 1
               B(I,L,0) = A(I,0,0)
               B(I,Q-L,0) = 0.0D0
   20       CONTINUE
            DO 60 J = 1, Q2
               INDX = MOD(J*L,Q)
               DO 40 I = 0, P - 1
                  B(I,L,0) = B(I,L,0) + A(I,0,J)*COSINE(0,INDX)
                  B(I,Q-L,0) = B(I,Q-L,0) - A(I,0,Q-J)*SINE(0,INDX)
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
         DO 100 I = 0, P - 1
            B(I,0,0) = A(I,0,0)
  100    CONTINUE
         DO 140 J = 1, Q2
            DO 120 I = 0, P - 1
               B(I,0,0) = B(I,0,0) + A(I,0,J)
  120       CONTINUE
  140    CONTINUE
         DO 180 J = 1, Q2
            DO 160 I = 0, P - 1
               TEMP = B(I,J,0)
               B(I,J,0) = B(I,J,0) + B(I,Q-J,0)
               B(I,Q-J,0) = TEMP - B(I,Q-J,0)
  160       CONTINUE
  180    CONTINUE
C
C        Code for general K --
C
         DO 460 K = 1, (R-1)/2
            KP = R - K
            DO 220 J = 1, Q2
               DO 200 I = 0, P - 1
                  TEMPR = A(I,K,J)
                  TEMPI = A(I,KP,Q-J-1)
                  A(I,K,J) = TEMPR + A(I,KP,J-1)
                  A(I,KP,Q-J-1) = TEMPI - A(I,K,Q-J)
                  A(I,KP,J-1) = TEMPR - A(I,KP,J-1)
                  A(I,K,Q-J) = -TEMPI - A(I,K,Q-J)
  200          CONTINUE
  220       CONTINUE
            DO 300 L = 1, Q2
               DO 240 I = 0, P - 1
                  B(I,L,K) = A(I,K,0)
                  B(I,L,KP) = A(I,KP,Q-1)
                  B(I,Q-L,K) = 0.0D0
                  B(I,Q-L,KP) = 0.0D0
  240          CONTINUE
               DO 280 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 260 I = 0, P - 1
                     B(I,L,K) = B(I,L,K) + A(I,K,J)*COSINE(0,INDX)
                     B(I,L,KP) = B(I,L,KP) + A(I,KP,Q-J-1)*COSINE(0,
     *                           INDX)
                     B(I,Q-L,K) = B(I,Q-L,K) + A(I,KP,J-1)*SINE(0,INDX)
                     B(I,Q-L,KP) = B(I,Q-L,KP) - A(I,K,Q-J)*SINE(0,INDX)
  260             CONTINUE
  280          CONTINUE
  300       CONTINUE
            DO 320 I = 0, P - 1
               B(I,0,K) = A(I,K,0)
               B(I,0,KP) = A(I,KP,Q-1)
  320       CONTINUE
            DO 360 J = 1, Q2
               DO 340 I = 0, P - 1
                  B(I,0,K) = B(I,0,K) + A(I,K,J)
                  B(I,0,KP) = B(I,0,KP) + A(I,KP,Q-J-1)
  340          CONTINUE
  360       CONTINUE
            DO 400 J = 1, Q2
               DO 380 I = 0, P - 1
                  TEMPR = B(I,J,K)
                  TEMPI = B(I,J,KP)
                  B(I,J,K) = TEMPR - B(I,Q-J,KP)
                  B(I,J,KP) = TEMPI + B(I,Q-J,K)
                  TEMP1 = B(I,Q-J,K)
                  B(I,Q-J,K) = TEMPR + B(I,Q-J,KP)
                  B(I,Q-J,KP) = TEMPI - TEMP1
  380          CONTINUE
  400       CONTINUE
            DO 440 J = 1, Q - 1
               DO 420 I = 0, P - 1
                  BR = B(I,J,K)
                  BI = B(I,J,KP)
                  B(I,J,K) = COSINE(K,J)*BR - SINE(K,J)*BI
                  B(I,J,KP) = COSINE(K,J)*BI + SINE(K,J)*BR
  420          CONTINUE
  440       CONTINUE
  460    CONTINUE
C
C        Code for K=R/2 when R is even --
C
         IF (MOD(R,2).EQ.0) THEN
            R2 = R/2
            DO 500 L = 1, Q2
               DO 480 I = 0, P - 1
                  B(I,L,R2) = A(I,R2,Q2)
                  B(I,Q-L,R2) = 0.0D0
  480          CONTINUE
  500       CONTINUE
            DO 560 L = 1, Q2
               DO 540 J = 0, Q2 - 1
                  INDX = MOD(L*(J+Q2+1),Q)
                  DO 520 I = 0, P - 1
                     B(I,L,R2) = B(I,L,R2) + A(I,R2,J)*COSINE(0,INDX)
                     B(I,Q-L,R2) = B(I,Q-L,R2) + A(I,R2,Q-J-1)*SINE(0,
     *                             INDX)
  520             CONTINUE
  540          CONTINUE
  560       CONTINUE
            DO 580 I = 0, P - 1
               B(I,0,R2) = A(I,R2,Q2)
  580       CONTINUE
            DO 620 J = 0, Q2 - 1
               DO 600 I = 0, P - 1
                  B(I,0,R2) = B(I,0,R2) + A(I,R2,J)
  600          CONTINUE
  620       CONTINUE
            DO 660 J = 1, Q2 - 1, 2
               DO 640 I = 0, P - 1
                  TEMPR = B(I,J,R2)
                  TEMPI = B(I,Q-J,R2)
                  B(I,J,R2) = TEMPI - TEMPR
                  B(I,Q-J,R2) = TEMPR + TEMPI
                  TEMPR = B(I,J+1,R2)
                  TEMPI = B(I,Q-J-1,R2)
                  B(I,J+1,R2) = TEMPR - TEMPI
                  B(I,Q-J-1,R2) = -TEMPR - TEMPI
  640          CONTINUE
  660       CONTINUE
            IF (MOD(Q2,2).EQ.1) THEN
               DO 680 I = 0, P - 1
                  TEMPR = B(I,Q2,R2)
                  TEMPI = B(I,Q2+1,R2)
                  B(I,Q2,R2) = TEMPI - TEMPR
                  B(I,Q2+1,R2) = TEMPR + TEMPI
  680          CONTINUE
            END IF
         END IF
C
      ELSE
C
         DO 1140 I = 0, P - 1
C
C           Code for K=0 --
C
            DO 720 L = 1, Q2
               B(I,L,0) = A(I,0,0)
               B(I,Q-L,0) = 0.0D0
               DO 700 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  B(I,L,0) = B(I,L,0) + A(I,0,J)*COSINE(0,INDX)
                  B(I,Q-L,0) = B(I,Q-L,0) - A(I,0,Q-J)*SINE(0,INDX)
  700          CONTINUE
  720       CONTINUE
            B(I,0,0) = A(I,0,0)
            DO 740 J = 1, Q2
               B(I,0,0) = B(I,0,0) + A(I,0,J)
  740       CONTINUE
            DO 760 J = 1, Q2
               TEMP = B(I,J,0)
               B(I,J,0) = B(I,J,0) + B(I,Q-J,0)
               B(I,Q-J,0) = TEMP - B(I,Q-J,0)
  760       CONTINUE
C
C           Code for general K --
C
            DO 800 J = 1, Q2
               DO 780 K = 1, (R-1)/2
                  TEMPR = A(I,K,J)
                  TEMPI = A(I,R-K,Q-J-1)
                  A(I,K,J) = TEMPR + A(I,R-K,J-1)
                  A(I,R-K,Q-J-1) = TEMPI - A(I,K,Q-J)
                  A(I,R-K,J-1) = TEMPR - A(I,R-K,J-1)
                  A(I,K,Q-J) = -TEMPI - A(I,K,Q-J)
  780          CONTINUE
  800       CONTINUE
            DO 880 L = 1, Q2
               DO 820 K = 1, (R-1)/2
                  B(I,L,K) = A(I,K,0)
                  B(I,L,R-K) = A(I,R-K,Q-1)
                  B(I,Q-L,K) = 0.0D0
                  B(I,Q-L,R-K) = 0.0D0
  820          CONTINUE
               DO 860 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 840 K = 1, (R-1)/2
                     B(I,L,K) = B(I,L,K) + A(I,K,J)*COSINE(0,INDX)
                     B(I,L,R-K) = B(I,L,R-K) + A(I,R-K,Q-J-1)*COSINE(0,
     *                            INDX)
                     B(I,Q-L,K) = B(I,Q-L,K) + A(I,R-K,J-1)*SINE(0,INDX)
                     B(I,Q-L,R-K) = B(I,Q-L,R-K) - A(I,K,Q-J)*SINE(0,
     *                              INDX)
  840             CONTINUE
  860          CONTINUE
  880       CONTINUE
            DO 900 K = 1, (R-1)/2
               B(I,0,K) = A(I,K,0)
               B(I,0,R-K) = A(I,R-K,Q-1)
  900       CONTINUE
            DO 940 J = 1, Q2
               DO 920 K = 1, (R-1)/2
                  B(I,0,K) = B(I,0,K) + A(I,K,J)
                  B(I,0,R-K) = B(I,0,R-K) + A(I,R-K,Q-J-1)
  920          CONTINUE
  940       CONTINUE
            DO 980 J = 1, Q2
               DO 960 K = 1, (R-1)/2
                  TEMPR = B(I,J,K)
                  TEMPI = B(I,J,R-K)
                  B(I,J,K) = TEMPR - B(I,Q-J,R-K)
                  B(I,J,R-K) = TEMPI + B(I,Q-J,K)
                  TEMP1 = B(I,Q-J,K)
                  B(I,Q-J,K) = TEMPR + B(I,Q-J,R-K)
                  B(I,Q-J,R-K) = TEMPI - TEMP1
  960          CONTINUE
  980       CONTINUE
            DO 1020 J = 1, Q - 1
               DO 1000 K = 1, (R-1)/2
                  BR = B(I,J,K)
                  BI = B(I,J,R-K)
                  B(I,J,K) = COSINE(K,J)*BR - SINE(K,J)*BI
                  B(I,J,R-K) = COSINE(K,J)*BI + SINE(K,J)*BR
 1000          CONTINUE
 1020       CONTINUE
C
C           Code for K=R/2 when R is even --
C
            IF (MOD(R,2).EQ.0) THEN
               R2 = R/2
               DO 1040 L = 1, Q2
                  B(I,L,R2) = A(I,R2,Q2)
                  B(I,Q-L,R2) = 0.0D0
 1040          CONTINUE
               DO 1080 L = 1, Q2
                  DO 1060 J = 0, Q2 - 1
                     INDX = MOD(L*(J+Q2+1),Q)
                     B(I,L,R2) = B(I,L,R2) + A(I,R2,J)*COSINE(0,INDX)
                     B(I,Q-L,R2) = B(I,Q-L,R2) + A(I,R2,Q-J-1)*SINE(0,
     *                             INDX)
 1060             CONTINUE
 1080          CONTINUE
               B(I,0,R2) = A(I,R2,Q2)
               DO 1100 J = 0, Q2 - 1
                  B(I,0,R2) = B(I,0,R2) + A(I,R2,J)
 1100          CONTINUE
               DO 1120 J = 1, Q2 - 1, 2
                  TEMPR = B(I,J,R2)
                  TEMPI = B(I,Q-J,R2)
                  B(I,J,R2) = TEMPI - TEMPR
                  B(I,Q-J,R2) = TEMPR + TEMPI
                  TEMPR = B(I,J+1,R2)
                  TEMPI = B(I,Q-J-1,R2)
                  B(I,J+1,R2) = TEMPR - TEMPI
                  B(I,Q-J-1,R2) = -TEMPR - TEMPI
 1120          CONTINUE
               IF (MOD(Q2,2).EQ.1) THEN
                  TEMPR = B(I,Q2,R2)
                  TEMPI = B(I,Q2+1,R2)
                  B(I,Q2,R2) = TEMPI - TEMPR
                  B(I,Q2+1,R2) = TEMPR + TEMPI
               END IF
            END IF
 1140    CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE C06FQS(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-699 (DEC 1989).
C
C     Radix six Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:5), B(0:P-1,0:5,0:R-1),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1, T1I, T1R, T2, T2I, T2R, T3, T3I, T3R, T4,
     *                  T5, T6, U0, U0I, U0R, U1, U1I, U1R, U2, U2I,
     *                  U2R, V0, V0I, V0R, V1, V1I, V1R, V2, V2I, V2R,
     *                  X0P, X1P, X2P, X3P, X4P, X5P, Y0P, Y1P, Y2P,
     *                  Y3P, Y4P, Y5P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,2)
         T2 = A(I,0,0) - 0.5D0*T1
         T3 = SIN60*A(I,0,4)
         U0 = A(I,0,0) + T1
         U1 = T2 + T3
         U2 = T2 - T3
         T1 = A(I,0,1)
         T2 = A(I,0,3) - 0.5D0*T1
         T3 = -SIN60*A(I,0,5)
         V0 = A(I,0,3) + T1
         V1 = T2 + T3
         V2 = T2 - T3
         B(I,0,0) = U0 + V0
         B(I,1,0) = U1 - V1
         B(I,2,0) = U2 + V2
         B(I,3,0) = U0 - V0
         B(I,4,0) = U1 + V1
         B(I,5,0) = U2 - V2
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,2) + A(I,R-K,1)
               T1I = A(I,R-K,3) - A(I,K,4)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,R-K,5) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,2)-A(I,R-K,1))
               T3I = SIN60*(A(I,R-K,3)+A(I,K,4))
               U0R = A(I,K,0) + T1R
               U0I = A(I,R-K,5) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = A(I,R-K,0) + A(I,K,1)
               T1I = -A(I,K,5) + A(I,R-K,4)
               T2R = A(I,R-K,2) - 0.5D0*T1R
               T2I = -A(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(A(I,R-K,0)-A(I,K,1))
               T3I = SIN60*(-A(I,K,5)-A(I,R-K,4))
               V0R = A(I,R-K,2) + T1R
               V0I = -A(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               X0P = U0R + V0R
               Y0P = U0I + V0I
               X1P = U1R - V1R
               Y1P = U1I - V1I
               X2P = U2R + V2R
               Y2P = U2I + V2I
               X3P = U0R - V0R
               Y3P = U0I - V0I
               X4P = U1R + V1R
               Y4P = U1I + V1I
               X5P = U2R - V2R
               Y5P = U2I - V2I
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
               B(I,4,K) = COSINE(K,4)*X4P - SINE(K,4)*Y4P
               B(I,4,R-K) = COSINE(K,4)*Y4P + SINE(K,4)*X4P
               B(I,5,K) = COSINE(K,5)*X5P - SINE(K,5)*Y5P
               B(I,5,R-K) = COSINE(K,5)*Y5P + SINE(K,5)*X5P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            C5K = COSINE(K,5)
            S5K = SINE(K,5)
            DO 80 I = 0, P - 1
               T1R = A(I,K,2) + A(I,KP,1)
               T1I = A(I,KP,3) - A(I,K,4)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,KP,5) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,2)-A(I,KP,1))
               T3I = SIN60*(A(I,KP,3)+A(I,K,4))
               U0R = A(I,K,0) + T1R
               U0I = A(I,KP,5) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = A(I,KP,0) + A(I,K,1)
               T1I = -A(I,K,5) + A(I,KP,4)
               T2R = A(I,KP,2) - 0.5D0*T1R
               T2I = -A(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(A(I,KP,0)-A(I,K,1))
               T3I = SIN60*(-A(I,K,5)-A(I,KP,4))
               V0R = A(I,KP,2) + T1R
               V0I = -A(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               X0P = U0R + V0R
               Y0P = U0I + V0I
               X1P = U1R - V1R
               Y1P = U1I - V1I
               X2P = U2R + V2R
               Y2P = U2I + V2I
               X3P = U0R - V0R
               Y3P = U0I - V0I
               X4P = U1R + V1R
               Y4P = U1I + V1I
               X5P = U2R - V2R
               Y5P = U2I - V2I
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
               B(I,4,K) = C4K*X4P - S4K*Y4P
               B(I,4,KP) = C4K*Y4P + S4K*X4P
               B(I,5,K) = C5K*X5P - S5K*Y5P
               B(I,5,KP) = C5K*Y5P + S5K*X5P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = A(I,R2,0) + A(I,R2,2)
            T2 = A(I,R2,5) + A(I,R2,3)
            T3 = A(I,R2,1) - 0.5D0*T1
            T4 = A(I,R2,4) + 0.5D0*T2
            T5 = SIN60*(A(I,R2,0)-A(I,R2,2))
            T6 = SIN60*(A(I,R2,5)-A(I,R2,3))
            B(I,0,R2) = A(I,R2,1) + T1
            B(I,1,R2) = T4 + T5
            B(I,2,R2) = T6 - T3
            B(I,3,R2) = T2 - A(I,R2,4)
            B(I,4,R2) = T3 + T6
            B(I,5,R2) = T4 - T5
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FQT(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-700 (DEC 1989).
C
C     Radix five Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  R54, SIN36, SIN72, S36S72
      PARAMETER         (R54=0.559016994374947424102293417182819D0,
     *                  SIN36=0.587785252292473129168705954639073D0,
     *                  SIN72=0.951056516295153572116439333379382D0,
     *                  S36S72=0.618033988749894848204586834365638D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:4), B(0:P-1,0:4,0:R-1),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T1, T10,
     *                  T10I, T10R, T11, T11I, T11R, T1I, T1R, T2, T2I,
     *                  T2R, T3, T3I, T3R, T4, T4I, T4R, T5, T5I, T5R,
     *                  T6, T6I, T6R, T7, T7I, T7R, T8, T8I, T8R, T9,
     *                  T9I, T9R, X0P, X1P, X2P, X3P, X4P, Y0P, Y1P,
     *                  Y2P, Y3P, Y4P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,1)
         T2 = A(I,0,2)
         T3 = SIN72*A(I,0,4)
         T4 = SIN72*A(I,0,3)
         T5 = T1 + T2
         T6 = R54*(T1-T2)
         T7 = A(I,0,0) - 0.25D0*T5
         T8 = T7 + T6
         T9 = T7 - T6
         T10 = T3 + S36S72*T4
         T11 = S36S72*T3 - T4
         B(I,0,0) = A(I,0,0) + T5
         B(I,1,0) = T8 + T10
         B(I,2,0) = T9 + T11
         B(I,3,0) = T9 - T11
         B(I,4,0) = T8 - T10
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,1) + A(I,R-K,0)
               T1I = A(I,R-K,3) - A(I,K,4)
               T2R = A(I,K,2) + A(I,R-K,1)
               T2I = A(I,R-K,2) - A(I,K,3)
               T3R = SIN72*(A(I,K,1)-A(I,R-K,0))
               T3I = SIN72*(A(I,R-K,3)+A(I,K,4))
               T4R = SIN72*(A(I,K,2)-A(I,R-K,1))
               T4I = SIN72*(A(I,R-K,2)+A(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,K,0) - 0.25D0*T5R
               T7I = A(I,R-K,4) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               X0P = A(I,K,0) + T5R
               Y0P = A(I,R-K,4) + T5I
               X1P = T8R + T10I
               Y1P = T8I - T10R
               X2P = T9R + T11I
               Y2P = T9I - T11R
               X3P = T9R - T11I
               Y3P = T9I + T11R
               X4P = T8R - T10I
               Y4P = T8I + T10R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
               B(I,4,K) = COSINE(K,4)*X4P - SINE(K,4)*Y4P
               B(I,4,R-K) = COSINE(K,4)*Y4P + SINE(K,4)*X4P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            DO 80 I = 0, P - 1
               T1R = A(I,K,1) + A(I,KP,0)
               T1I = A(I,KP,3) - A(I,K,4)
               T2R = A(I,K,2) + A(I,KP,1)
               T2I = A(I,KP,2) - A(I,K,3)
               T3R = SIN72*(A(I,K,1)-A(I,KP,0))
               T3I = SIN72*(A(I,KP,3)+A(I,K,4))
               T4R = SIN72*(A(I,K,2)-A(I,KP,1))
               T4I = SIN72*(A(I,KP,2)+A(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,K,0) - 0.25D0*T5R
               T7I = A(I,KP,4) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               X0P = A(I,K,0) + T5R
               Y0P = A(I,KP,4) + T5I
               X1P = T8R + T10I
               Y1P = T8I - T10R
               X2P = T9R + T11I
               Y2P = T9I - T11R
               X3P = T9R - T11I
               Y3P = T9I + T11R
               X4P = T8R - T10I
               Y4P = T8I + T10R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
               B(I,4,K) = C4K*X4P - S4K*Y4P
               B(I,4,KP) = C4K*Y4P + S4K*X4P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,R2,0) + A(I,R2,1)
C           T2 = 0.25D0*T1 - A(I,R2,2)
C           T3 = R54*(A(I,R2,0)-A(I,R2,1))
C           T4 = SIN36*A(I,R2,4) + SIN72*A(I,R2,3)
C           T5 = SIN72*A(I,R2,4) - SIN36*A(I,R2,3)
C           T6 = T2 + T3
C           T7 = T2 - T3
C           B(I,0,R2) = T1 + A(I,R2,2)
C           B(I,1,R2) = T4 + T6
C           B(I,2,R2) = T5 - T7
C           B(I,3,R2) = T5 + T7
C           B(I,4,R2) = T4 - T6
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FQU(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-701 (DEC 1989).
C
C     Radix four Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  ROOT2I
      PARAMETER         (ROOT2I=0.707106781186547524400844362104849D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:3), B(0:P-1,0:3,0:R-1),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, CK, S2K, S3K, SK, T1, T1I, T1R, T2,
     *                  T2I, T2R, T3, T3I, T3R, T4, T4I, T4R, X0P, X1P,
     *                  X2P, X3P, Y0P, Y1P, Y2P, Y3P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,0) + A(I,0,2)
         T2 = A(I,0,1)
         T3 = A(I,0,0) - A(I,0,2)
         T4 = A(I,0,3)
         B(I,0,0) = T1 + T2
         B(I,1,0) = T3 + T4
         B(I,2,0) = T1 - T2
         B(I,3,0) = T3 - T4
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,0) + A(I,R-K,1)
               T1I = A(I,R-K,3) - A(I,K,2)
               T2R = A(I,K,1) + A(I,R-K,0)
               T2I = A(I,R-K,2) - A(I,K,3)
               T3R = A(I,K,0) - A(I,R-K,1)
               T3I = A(I,R-K,3) + A(I,K,2)
               T4R = A(I,K,1) - A(I,R-K,0)
               T4I = A(I,R-K,2) + A(I,K,3)
               X0P = T1R + T2R
               Y0P = T1I + T2I
               X1P = T3R + T4I
               Y1P = T3I - T4R
               X2P = T1R - T2R
               Y2P = T1I - T2I
               X3P = T3R - T4I
               Y3P = T3I + T4R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            DO 80 I = 0, P - 1
               T1R = A(I,K,0) + A(I,KP,1)
               T1I = A(I,KP,3) - A(I,K,2)
               T2R = A(I,K,1) + A(I,KP,0)
               T2I = A(I,KP,2) - A(I,K,3)
               T3R = A(I,K,0) - A(I,KP,1)
               T3I = A(I,KP,3) + A(I,K,2)
               T4R = A(I,K,1) - A(I,KP,0)
               T4I = A(I,KP,2) + A(I,K,3)
               X0P = T1R + T2R
               Y0P = T1I + T2I
               X1P = T3R + T4I
               Y1P = T3I - T4R
               X2P = T1R - T2R
               Y2P = T1I - T2I
               X3P = T3R - T4I
               Y3P = T3I + T4R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            B(I,0,R2) = A(I,R2,0) + A(I,R2,1)
            B(I,2,R2) = A(I,R2,3) - A(I,R2,2)
            T3 = A(I,R2,0) - A(I,R2,1)
            T4 = A(I,R2,3) + A(I,R2,2)
            B(I,1,R2) = ROOT2I*(T3+T4)
            B(I,3,R2) = -ROOT2I*(T3-T4)
  120    CONTINUE
      END IF
C
      RETURN
      END
      SUBROUTINE C06FQV(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-702 (DEC 1989).
C
C     Radix three Hermitian to real Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:2), B(0:P-1,0:2,0:R-1),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, CK, S2K, SK, T1, T1I, T1R, T2, T2I, T2R,
     *                  T3, T3I, T3R, X0P, X1P, X2P, Y0P, Y1P, Y2P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,1)
         T2 = A(I,0,0) - 0.5D0*T1
         T3 = SIN60*A(I,0,2)
         B(I,0,0) = A(I,0,0) + T1
         B(I,1,0) = T2 + T3
         B(I,2,0) = T2 - T3
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,1) + A(I,R-K,0)
               T1I = A(I,R-K,1) - A(I,K,2)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,R-K,2) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,1)-A(I,R-K,0))
               T3I = SIN60*(A(I,R-K,1)+A(I,K,2))
               X0P = A(I,K,0) + T1R
               Y0P = A(I,R-K,2) + T1I
               X1P = T2R + T3I
               Y1P = T2I - T3R
               X2P = T2R - T3I
               Y2P = T2I + T3R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = SINE(K,1)*X1P + COSINE(K,1)*Y1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = SINE(K,2)*X2P + COSINE(K,2)*Y2P
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            DO 80 I = 0, P - 1
               T1R = A(I,K,1) + A(I,KP,0)
               T1I = A(I,KP,1) - A(I,K,2)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,KP,2) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,1)-A(I,KP,0))
               T3I = SIN60*(A(I,KP,1)+A(I,K,2))
               X0P = A(I,K,0) + T1R
               Y0P = A(I,KP,2) + T1I
               X1P = T2R + T3I
               Y1P = T2I - T3R
               X2P = T2R - T3I
               Y2P = T2I + T3R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = SK*X1P + CK*Y1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = S2K*X2P + C2K*Y2P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,R2,0) + A(I,R2,1)
C           T2 = 0.5D0*A(I,R2,0) - A(I,R2,1)
C           T3 = SIN60*A(I,R2,2)
C           B(I,0,R2) = T1
C           B(I,1,R2) = T2 + T3
C           B(I,2,R2) = -T2 + T3
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FQW(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-703 (DEC 1989).
C
C     Radix two Hermitian to real Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:1), B(0:P-1,0:1,0:R-1),
     *                  COSINE(0:R-1), SINE(0:R-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  CK, SK, X1HAT, Y1HAT
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         B(I,0,0) = A(I,0,0) + A(I,0,1)
         B(I,1,0) = A(I,0,0) - A(I,0,1)
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1HAT = A(I,K,0) - A(I,R-K,0)
               Y1HAT = A(I,R-K,1) + A(I,K,1)
               B(I,0,K) = A(I,K,0) + A(I,R-K,0)
               B(I,0,R-K) = A(I,R-K,1) - A(I,K,1)
               B(I,1,K) = COSINE(K)*X1HAT - SINE(K)*Y1HAT
               B(I,1,R-K) = COSINE(K)*Y1HAT + SINE(K)*X1HAT
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K)
            SK = SINE(K)
            DO 80 I = 0, P - 1
               X1HAT = A(I,K,0) - A(I,KP,0)
               Y1HAT = A(I,KP,1) + A(I,K,1)
               B(I,0,K) = A(I,K,0) + A(I,KP,0)
               B(I,0,KP) = A(I,KP,1) - A(I,K,1)
               B(I,1,K) = CK*X1HAT - SK*Y1HAT
               B(I,1,KP) = CK*Y1HAT + SK*X1HAT
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           B(I,0,R2) = A(I,R2,0)
C           B(I,1,R2) = A(I,R2,1)
C 120    CONTINUE
C     END IF
C
      RETURN
      END
      SUBROUTINE C06FQX(A,B,M,N,Q,NQ,TRIG)
CVD$R VECTOR
c CVD$R NOLSTVAL
c c CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Hermitian to Real Fast Fourier Transform Kernel Driver
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:M*N-1), B(0:M*N-1), TRIG(0:2*N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR
      INTEGER           I, P, QI, R
      LOGICAL           INA
C     .. External Subroutines ..
      EXTERNAL          C06FQQ, C06FQR, C06FQS, C06FQT, C06FQU, C06FQV,
     *                  C06FQW
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT, DBLE
C     .. Executable Statements ..
      INA = .TRUE.
      P = 1
      R = N
      IF (N.EQ.1) RETURN
      CALL C06FQQ(A,M,N)
      DO 20 I = 1, NQ
         QI = Q(I)
         R = R/QI
         IF (INA) THEN
            IF (QI.EQ.2) THEN
               CALL C06FQW(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FQV(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FQU(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FQT(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FQS(A,B,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FQR(A,B,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         ELSE
            IF (QI.EQ.2) THEN
               CALL C06FQW(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.3) THEN
               CALL C06FQV(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.4) THEN
               CALL C06FQU(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.5) THEN
               CALL C06FQT(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE IF (QI.EQ.6) THEN
               CALL C06FQS(B,A,M*P,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)*QI*R)
     *                     )
            ELSE
               CALL C06FQR(B,A,M*P,QI,R,TRIG((P-1)*QI*R),TRIG(N+(P-1)
     *                     *QI*R))
            END IF
         END IF
         INA = .NOT. INA
         P = P*QI
   20 CONTINUE
C
      FACTOR = 2.0D0/SQRT(DBLE(N))
      IF (INA) THEN
         DO 40 I = 0, M*N - 1
            A(I) = A(I)*FACTOR
   40    CONTINUE
      ELSE
         DO 60 I = 0, M*N - 1
            A(I) = B(I)*FACTOR
   60    CONTINUE
      END IF
C
      RETURN
      END

      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
      DOUBLE PRECISION FUNCTION X01AAF(X)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE REAL CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Executable Statements ..
      X01AAF = 3.14159265358979323846264338328D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02ABF(X)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C     I.E. RETURNS THE SAME VALUE AS X02AKF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION          X
C     .. External Functions ..
      DOUBLE PRECISION          X02AKF
      EXTERNAL                  X02AKF
C     .. Executable Statements ..
      X02ABF = X02AKF()
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /1.11022302462516D-16 /
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AKF = X02CON
      RETURN
      END
      DOUBLE PRECISION FUNCTION X02AMF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE 'SAFE RANGE' PARAMETER
C     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
C     FOR ANY X WHICH SATISFIES X.GE.Z AND X.LE.1/Z
C     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
C     ERROR
C
C        -X
C        1.0/X
C        SQRT(X)
C        LOG(X)
C        EXP(LOG(X))
C        Y**(LOG(X)/LOG(Y)) FOR ANY Y
C
      DOUBLE PRECISION X02CON
      DATA X02CON /2.22507385850721D-308 /
C     .. Executable Statements ..
      X02AMF = X02CON
      RETURN
      END
      INTEGER FUNCTION X02BHF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, B.
C
C     .. Executable Statements ..
      X02BHF =     2
      RETURN
      END
      INTEGER FUNCTION X02BJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS THE MODEL PARAMETER, p.
C
C     .. Executable Statements ..
      X02BJF =    53
      RETURN
      END
      LOGICAL FUNCTION X02DAF(X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS .FALSE. IF THE SYSTEM SETS UNDERFLOWING QUANTITIES
C     TO ZERO, WITHOUT ANY ERROR INDICATION OR UNDESIRABLE WARNING
C     OR SYSTEM OVERHEAD.
C     RETURNS .TRUE. OTHERWISE, IN WHICH CASE CERTAIN LIBRARY
C     ROUTINES WILL TAKE SPECIAL PRECAUTIONS TO AVOID UNDERFLOW
C     (USUALLY AT SOME COST IN EFFICIENCY).
C
C     X IS A DUMMY ARGUMENT
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        X
C     .. Executable Statements ..
      X02DAF = .FALSE.
      RETURN
      END

      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
      INTEGER FUNCTION G05DYF(M,N)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 7 REVISED IER-135 (DEC 1978)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS AN INTEGER RESULT, UNIFORMLY DISTRIBUTED
C     BETWEEN M AND N INCLUSIVE.
C     THE CONTORTED PROGRAMMING IS TO GIVE CORRECT RESULTS AND TO
C     AVOID DIAGNOSTICS IN CASES WHERE:
C     1) ROUNDING CAUSES OVERFLOW OF THE RANGE (M,N).
C     2) INTEGER/REAL CONVERSION IS NOT EXACT.
C     3) INT TRUNCATES TOWARDS MINUS INFINITY.
C     .. Scalar Arguments ..
      INTEGER                 M, N
C     .. Local Scalars ..
      DOUBLE PRECISION        ONE, X, Y, Z
      INTEGER                 I, J, K
C     .. External Functions ..
      DOUBLE PRECISION        G05CAF
      EXTERNAL                G05CAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MAX, MIN, DBLE, INT
C     .. Data statements ..
      DATA                    ONE/1.0D0/
C     .. Executable Statements ..
      Y = DBLE(M)
      Z = DBLE(N)
      X = MIN(Y,Z) + (ABS(Z-Y)+ONE)*G05CAF(X)
      I = INT(X)
      IF (DBLE(I).GT.X) I = I - 1
C     ON A MACHINE WHICH OVERFLOWS WHEN COMPARING THE LARGEST (MOST
C     POSITIVE) AND SMALLEST (MOST NEGATIVE) INTEGERS, THE
C     FOLLOWING STATEMENTS SHOULD BE INSERTED HERE:
C     IF (M) 50, 140, 55
C     50 IF (N) 140, 60, 60
C     55 IF (N) 80, 80, 140
C     60 J = M
C     K = N
C     GO TO 100
C     80 J = N
C     K = M
C     100 IF (I.GT.0) GO TO 120
C     G05DYF = MAX0(I,J)
C     RETURN
C     120 G05DYF = MIN0(I,K)
C     RETURN
C     140 CONTINUE
      IF (M.GT.N) GO TO 20
      J = M
      K = N
      GO TO 40
   20 J = N
      K = M
   40 G05DYF = MIN(MAX(I,J),K)
      RETURN
      END
      SUBROUTINE G05EJF(IA,N,IZ,M,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G05EJF TAKES A RANDOM SAMPLE OF SIZE M FROM IA (OF
C     SIZE N) AND PUTS IT INTO IZ.
C
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      INTEGER           IA(N), IZ(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, P, Q
      INTEGER           I, IERR, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              ONE/1.0D0/
C     .. Executable Statements ..
      IERR = 1
      IF (N.LT.1) GO TO 60
      IERR = 2
      IF ((M.LT.1) .OR. (M.GT.N)) GO TO 60
      P = DBLE(M)
      Q = DBLE(N)
      J = 1
      DO 40 I = 1, N
         IF (Q*G05CAF(0.0D0).GT.P) GO TO 20
         IZ(J) = IA(I)
         J = J + 1
         P = P - ONE
   20    Q = Q - ONE
   40 CONTINUE
      IFAIL = 0
      RETURN
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE G05CAY(REINIT)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     called by G05CAF, G05FAF, G05FBF or G05FDF when needed, to fill
C     the internal array RV in COMMON block CG05CA with new
C     pseudo-random numbers.
C
C     G05CAY uses a multiplicative congruential algorithm
C
C     N := N * 13**13 modulo 2**59
C
C     where N is a notional variable internal to G05CAY. The value of N
C     is converted to a real number in the range 0.0 to 1.0 by scaling
C     by 2**(-59), with care taken that the result lies strictly
C     between 0.0 and 1.0.
C
C     N is initially set to 123456789*(2**32+1) but can be changed
C     by a call to G05CBF or G05CCF.
C
C     G05CAY generates number 63 at a time, in order to achieve
C     efficiency on vector-processing machines. The first call of
C     G05CAY generates 63 consecutive values of N, N(i), i = 1,...,63.
C     Subsequent calls generate the next set of 63 values of N by
C
C     N(i) := N(i) * (13**13)**63 modulo 2**59, for i = 1,...,63.
C
C     The value 63 is defined as the symbol LV in a parameter statement
C     in each routine which needs it. The particular value 63 was
C     chosen because of special properties of the multiplier
C     (13**13)**63 modulo 2**59, which permit efficient multi-length
C     arithmetic when ILIM = 4 (see below). Only a few values of LV
C     have such properties.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     These notes are intended to guide implementors through the text
C     changes necessary to implement the basic random number generator
C     routines G05CAY, G05CAZ, G05CBF, G05CCF, G05CFZ, G05CGZ. Please
C     follow these guidelines, and consult NAG Central Office if in any
C     doubt or difficulty. Please send a listing of your final text for
C     these routines to Central Office.
C
C     1.  Read "DETAILS-NOTE-1" below.
C         Decide the relevant value of ILIM, say nn, taking account of
C         the suggestion for 'long' integers.
C
C     2.  Activate all lines beginning CAnn.
C
C     3.  Read "DETAILS-NOTE-2" below.
C         Check whether your compiler has the functions ISHFT and IAND
C         (or equivalent functions) and compiles inline code for them.
C
C     4.  If ISHFT and IAND or equivalent functions are available as
C         inline functions, activate all lines beginning CYnn. If
C         necessary, change the function names. Otherwise activate all
C         lines beginning CXnn.
C
C     ******************************************************************
C
C     ************************ DETAILS-NOTE-1 **************************
C
C     The algorithm requires that the values of N and of the multi-
C     plier 13**13 be stored as 59-bit unsigned integers and that
C     the least significant 59 bits of their product be computed. On
C     most machines this can be done much more efficiently in
C     machine code than in Fortran. The Fortran code given here is
C     intended to give guidance on a machine code implementation,
C     and to provide a less efficient implementation as a fall-back.
C
C     The 59-bit integer N is stored as a multiple-length integer in
C     the array B. In fact for convenience the 60-bit integer 2*N is
C     stored. The multiplier 13**13 is stored in the array M.
C     The multiplier (13**13)**63 modulo 2**59 is stored in the array
C     MLV in exactly the same way as the basic multiplier is stored in
C     the array M.
C
C     The number of elements in N and M (ILIM) and the number of bits
C     used in each element of N and M (IBITS) depend on the number
C     of bits (including sign) in an integer variable as follows -
C
C        ILIM     IBITS     number of bits in integer variable
C          4        15                 .ge. 32
C          3        20                 .ge. 41
C          2        30                 .ge. 60
C
C     For greatest efficiency ILIM should be chosen as small as
C     possible.
C
C     N.B. the most significant bits of N are stored in B(I,ILIM),
C     the next most significant bits in B(I,ILIM-1), . . . , and
C     the least significant bits in B(I,1). The multiplier is stored
C     in M(ILIM), M(ILIM-1), . . . , M(1) in the same way.
C
C     Note -
C
C     1) in the above table the value of IBITS is less than half the
C     number of bits in an integer variable. This ensures that the
C     necessary integer products can be formed and summed correctly
C     without integer overflow. However many machines have instruc-
C     tions for forming double-length integer products. A machine
C     code implementation can take advantage of this and allow IBITS
C     to be as large (or almost as large) as the number of bits in
C     an integer variable and ILIM to be correspondingly smaller.
C     This should be much more efficient.
C
C     2) the figures in the rightmost column in the above table are
C     correct for the specific value of the multiplier. They are
C     certainly not correct for arbitrary 60-bit arithmetic.
C
C     3) it may well be advantageous to use 'long' integers, if
C     available, within G05CAY, even if they are not used
C     elsewhere in the library.
C
C     Variant code for the array declarations and data statements
C     is supplied in comments beginning CAnn where the digits nn are
C     the value of ILIM.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LV
      PARAMETER         (LV=63)
      INTEGER           ILIM
      PARAMETER         (ILIM=4)
CA03  PARAMETER         (ILIM=3)
CA02  PARAMETER         (ILIM=2)
      DOUBLE PRECISION  ONE, R2
      PARAMETER         (ONE=1.0D0,R2=0.5D0)
      DOUBLE PRECISION  RP1, RP2
      PARAMETER         (RP1=R2**60,RP2=R2**30)
CA03  DOUBLE PRECISION  RP1, RP2
CA03  PARAMETER         (RP1=R2**60,RP2=R2**40)
CA02  DOUBLE PRECISION  RP1
CA02  PARAMETER         (RP1=R2**60)
C     .. Scalar Arguments ..
      LOGICAL           REINIT
C     .. Scalars in Common ..
      INTEGER           DEFOPT, OPTION, POSSOP, KV
C     .. Arrays in Common ..
      DOUBLE PRECISION  RV(LV)
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONEM
      INTEGER           I, T1, T2, T3, U, V
CX03  INTEGER           I, T1, T2, U, V
CX02  INTEGER           I, T1, U, V
CY04  INTEGER           I, T1, T2, T3, T4
CY03  INTEGER           I, T1, T2, T3
CY02  INTEGER           I, T1, T2
      LOGICAL           INIT
C     .. Local Arrays ..
      INTEGER           M(ILIM), MLV(ILIM)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /CG05CA/RV, KV
C     .. Save statement ..
      SAVE              /AG05CA/, /CG05CA/, ONEM, INIT
C     .. Data statements ..
      DATA              INIT / .TRUE. /
      DATA              M /
     *                  17917, 13895, 19930,     8 /
CA03 *                  247293, 485810,    275 /
CA02 *                  455329277,    282074 /
      DATA              MLV /
     *                  85,  3703,  6070,  6822 /
CA03 *                  753749, 972915, 218309 /
CA02 *                  121339989, 223549366 /
C     .. Executable Statements ..
C
C     ************************ DETAILS-NOTE-2 **************************
C
C     It is advantageous to use non-standard Fortran intrinsic
C     functions for shifting and masking if these are available and if
C     they are compiled as in-line code without the overhead of a
C     subroutine call. Alternative code is given which uses the integer
C     functions:
C
C     ISHFT(I,J) to shift I J bits to the left (a negative value of
C                 J indicating a right shift)
C     IAND(I,J)  to form the logical and of I and J
C
C     It may be necesssary to replace these by calls to different
C     intrinsic functions provided by the fortran compiler.
C
C     Variant code for this computation is supplied in comments
C     beginning CXnn (using only arithmetic operations) or in
C     comments beginning CYnn (using shifting and masking functions)
C     where the digits nn are the value of ILIM.
C
C     ******************************************************************
C
      IF (INIT.OR.REINIT) THEN
         INIT = .FALSE.
         ONEM = ONE - X02AJF()
C
C        Generate first buffer of LV integers by multiplying
C        recursively by M modulo 2**59.
C        This loop cannot be vectorized.
C
         DO 20 I = 1, LV
            V = B(I-1,1)*M(1)
            U = V/32768
            B(I,1) = V - 32768*U
            V = U + B(I-1,2)*M(1) + B(I-1,1)*M(2)
            U = V/32768
            B(I,2) = V - 32768*U
            V = U + B(I-1,3)*M(1) + B(I-1,2)*M(2) + B(I-1,1)*M(3)
            U = V/32768
            B(I,3) = V - 32768*U
            V = U + B(I-1,4)*M(1) + B(I-1,3)*M(2) + B(I-1,2)*M(3)
     *            + B(I-1,1)*M(4)
            U = V/32768
            B(I,4) = V - 32768*U
CX03        V = B(I-1,1)*M(1)
CX03        U = V/1048576
CX03        B(I,1) = V - 1048576*U
CX03        V = U + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CX03        U = V/1048576
CX03        B(I,2) = V - 1048576*U
CX03        V = U + B(I-1,3)*M(1) + B(I-1,2)*M(2) + B(I-1,1)*M(3)
CX03        U = V/1048576
CX03        B(I,3) = V - 1048576*U
CX02        V = B(I-1,1)*M(1)
CX02        U = V/1073741824
CX02        B(I,1) = V - 1073741824*U
CX02        V = U + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CX02        U = V/1073741824
CX02        B(I,2) = V - 1073741824*U
CY04        T1 = B(I-1,1)*M(1)
CY04        T2 = ISHFT(T1,-15) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CY04        T3 = ISHFT(T2,-15) + B(I-1,3)*M(1) + B(I-1,2)*M(2)
CY04 *                         + B(I-1,1)*M(3)
CY04        T4 = ISHFT(T3,-15) + B(I-1,4)*M(1) + B(I-1,3)*M(2)
CY04 *                         + B(I-1,2)*M(3) + B(I-1,1)*M(4)
CY04        B(I,4) = IAND(T4,32767)
CY04        B(I,3) = IAND(T3,32767)
CY04        B(I,2) = IAND(T2,32767)
CY04        B(I,1) = IAND(T1,32767)
CY03        T1 = B(I-1,1)*M(1)
CY03        T2 = ISHFT(T1,-20) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CY03        T3 = ISHFT(T2,-20) + B(I-1,3)*M(1) + B(I-1,2)*M(2)
CY03 *                         + B(I-1,1)*M(3)
CY03        B(I,3) = IAND(T3,1048575)
CY03        B(I,2) = IAND(T2,1048575)
CY03        B(I,1) = IAND(T1,1048575)
CY02        T1 = B(I-1,1)*M(1)
CY02        T2 = ISHFT(T1,-30) + B(I-1,2)*M(1) + B(I-1,1)*M(2)
CY02        B(I,2) = IAND(T2,1073741823)
CY02        B(I,1) = IAND(T1,1073741823)
   20    CONTINUE
      ELSE
C
C        Generate next buffer of LV integers by multiplying in
C        parallel by M**LV modulo 2**59.
C
         DO 40 I = 1, LV
            V = B(I,1)*MLV(1)
            U = V/32768
            T1 = V - 32768*U
            V = U + B(I,2)*MLV(1) + B(I,1)*MLV(2)
            U = V/32768
            T2 = V - 32768*U
            V = U + B(I,3)*MLV(1) + B(I,2)*MLV(2) + B(I,1)*MLV(3)
            U = V/32768
            T3 = V - 32768*U
            V = U + B(I,4)*MLV(1) + B(I,3)*MLV(2) + B(I,2)*MLV(3)
     *            + B(I,1)*MLV(4)
            U = V/32768
            B(I,4) = V - 32768*U
            B(I,3) = T3
            B(I,2) = T2
            B(I,1) = T1
CX03        V = B(I,1)*MLV(1)
CX03        U = V/1048576
CX03        T1 = V - 1048576*U
CX03        V = U + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CX03        U = V/1048576
CX03        T2 = V - 1048576*U
CX03        V = U + B(I,3)*MLV(1) + B(I,2)*MLV(2) + B(I,1)*MLV(3)
CX03        U = V/1048576
CX03        B(I,3) = V - 1048576*U
CX03        B(I,2) = T2
CX03        B(I,1) = T1
CX02        V = B(I,1)*MLV(1)
CX02        U = V/1073741824
CX02        T1 = V - 1073741824*U
CX02        V = U + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CX02        U = V/1073741824
CX02        B(I,2) = V - 1073741824*U
CX02        B(I,1) = T1
CY04        T1 = B(I,1)*MLV(1)
CY04        T2 = ISHFT(T1,-15) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CY04        T3 = ISHFT(T2,-15) + B(I,3)*MLV(1) + B(I,2)*MLV(2)
CY04 *                         + B(I,1)*MLV(3)
CY04        T4 = ISHFT(T3,-15) + B(I,4)*MLV(1) + B(I,3)*MLV(2)
CY04 *                         + B(I,2)*MLV(3) + B(I,1)*MLV(4)
CY04        B(I,4) = IAND(T4,32767)
CY04        B(I,3) = IAND(T3,32767)
CY04        B(I,2) = IAND(T2,32767)
CY04        B(I,1) = IAND(T1,32767)
CY03        T1 = B(I,1)*MLV(1)
CY03        T2 = ISHFT(T1,-20) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CY03        T3 = ISHFT(T2,-20) + B(I,3)*MLV(1) + B(I,2)*MLV(2)
CY03 *                         + B(I,1)*MLV(3)
CY03        B(I,3) = IAND(T3,1048575)
CY03        B(I,2) = IAND(T2,1048575)
CY03        B(I,1) = IAND(T1,1048575)
CY02        T1 = B(I,1)*MLV(1)
CY02        T2 = ISHFT(T1,-30) + B(I,2)*MLV(1) + B(I,1)*MLV(2)
CY02        B(I,2) = IAND(T2,1073741823)
CY02        B(I,1) = IAND(T1,1073741823)
   40    CONTINUE
      END IF
C
C     Convert integers in B to real numbers in (0.0,1.0) stored in RV.
C
      DO 60 I = 1, LV
         RV(I) = MIN(ONEM,(B(I,4)*32768+B(I,3))*RP2
     *                     +(B(I,2)*32768+B(I,1))*RP1)
CX03     RV(I) = MIN(ONEM,(B(I,3)*1048576+B(I,2))*RP2 + B(I,1)*RP1)
CX02     RV(I) = MIN(ONEM,(B(I,2)*1073741824+B(I,1))*RP1)
CY04     RV(I) = MIN(ONEM,(ISHFT(B(I,4),15)+B(I,3))*RP2
CY04 *                     +(ISHFT(B(I,2),15)+B(I,1))*RP1)
CY03     RV(I) = MIN(ONEM,(ISHFT(B(I,3),20)+B(I,2))*RP2 + B(I,1)*RP1)
CY02     RV(I) = MIN(ONEM,(ISHFT(B(I,2),30)+B(I,1))*RP1)
   60 CONTINUE
      KV = 0
C
      RETURN
      END
      SUBROUTINE G05CAZ(INIT)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     called by G05CAF, G05CBF, G05CCF, G05CFZ, G05CGZ, G05DGF, G05FAF,
C     G05FBF AND G05FDF to ensure that the contents of common blocks
C     /AG05CA/, /BG05CA/, /CG05CA/ and /DG05CA/ are initialized.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     This version of G05CAZ must be used in conjunction with the
C     new auxiliary routine G05CAY which has been introduced at Mark 14.
C
C     These notes are intended to guide implementors through the text
C     changes necessary to implement the basic random number generator
C     routines G05CAY, G05CAZ, G05CBF, G05CCF, G05CFZ, G05CGZ. Please
C     follow these guidelines, and consult NAG Central Office if in any
C     doubt or difficulty. Please send a listing of your final text for
C     these routines to Central Office.
C
C     1.  Prepare code for G05CAY following guidelines supplied there.
C
C     2.  Read "DETAILS-NOTE-1" below.
C
C     3.  Activate all lines beginning CAnn, where nn is the value of
C         ILIM used in G05CAY.
C
C     ******************************************************************
C
C     ************************ DETAILS-NOTE-1 **************************
C
C     G05CAZ must be implemented consistently with G05CAY.
C
C     If G05CAY has been implemented simply by selecting suitable
C     variant code according to the value of ILIM, then a consistent
C     implementation of G05CAY may be obtained by using the variant
C     code supplied in comments beginning CAnn where the digits nn
C     are the value of ILIM.
C
C     If G05CAY has been implemented in machine code, it will still
C     be possible on many machines to implement G05CAZ in Fortran
C     and this will be satisfactory since it is not important for
C     G05CAZ to be particularly efficient. Essentially the code for
C     G05CAZ depends only on how the internal variable N is stored in
C     the array B in the common block /AG05CA/ and the code given
C     below should be applicable provided that N is stored in
C     accordance with a particular value of ILIM as defined in the
C     text of G05CAY.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LV
      PARAMETER         (LV=63)
      INTEGER           ILIM
      PARAMETER         (ILIM=4)
CA03  PARAMETER         (ILIM=3)
CA02  PARAMETER         (ILIM=2)
C     .. Scalar Arguments ..
      LOGICAL           INIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  GAMMA, NORMAL, VNORML
      INTEGER           DEFOPT, OPTION, POSSOP
C     .. Arrays in Common ..
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      LOGICAL           INIT2
C     .. External Subroutines ..
      EXTERNAL          G05CAY
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /BG05CA/NORMAL, GAMMA
      COMMON            /DG05CA/VNORML
C     .. Save statement ..
      SAVE              INIT2, /AG05CA/, /BG05CA/, /DG05CA/
C     .. Data statements ..
      DATA              INIT2/.TRUE./
C     .. Executable Statements ..
C
C     If INIT2 is not already .FALSE. , initialize /AG05CA/, /BG05CA/
C     and /DG05CA/ and set INIT2 to .FALSE.
C
      IF (INIT2) THEN
C
         B(0,1) =  6698
         B(0,2) =  7535
         B(0,3) = 26792
         B(0,4) = 30140
CA03     B(0,1) = 498218
CA03     B(0,2) = 172267
CA03     B(0,3) = 964506
CA02     B(0,1) = 246913578
CA02     B(0,2) = 987654312
         OPTION = 0
         DEFOPT = 0
         POSSOP = 0
C
         NORMAL = 1.0D0
         GAMMA = -1.0D0
         VNORML = 256.0D0
C
         INIT2 = .FALSE.
C
C        Initialize the buffer
C
         CALL G05CAY(.TRUE.)
      END IF
C
C     Set INIT to .FALSE. in any case
C
      INIT = .FALSE.
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION G05CAF(X)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     Returns a pseudo-random number uniformly distributed between
C     A and B.
C
C     Pseudo-random numbers are generated by the auxiliary routine
C     G05CAY, 63 at a time, and stored in the array RV in common block
C     CG05CA. G05CAF copies one number from the array RV into X,
C     calling G05CAY to replenish RV when necessary.
C
C     This revised version of G05CAF has been introduced for
C     compatibility with the new routines G05FAF, G05FBF and G05FDF,
C     introduced at Mark 14.
C
C     Jeremy Du Croz, NAG Ltd, June 1989.
C
C     .. Parameters ..
      INTEGER                          LV
      PARAMETER                        (LV=63)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Scalars in Common ..
      INTEGER                          KV
C     .. Arrays in Common ..
      DOUBLE PRECISION                 RV(LV)
C     .. Local Scalars ..
      LOGICAL                          INIT
C     .. External Subroutines ..
      EXTERNAL                         G05CAY, G05CAZ
C     .. Common blocks ..
      COMMON                           /CG05CA/RV, KV
C     .. Save statement ..
      SAVE                             INIT, /CG05CA/
C     .. Data statements ..
      DATA                             INIT/.TRUE./
C     .. Executable Statements ..
C
C     Ensure that KV in common block /CG05CA/ has been initialized
C
      IF (INIT) CALL G05CAZ(INIT)
C
C     Replenish the buffer if necessary
C
      IF (KV.GE.LV) CALL G05CAY(.FALSE.)
C
      KV = KV + 1
      G05CAF = RV(KV)
      RETURN
      END

