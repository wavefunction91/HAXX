      Subroutine HGEMMRR(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     $  B, LDB, BETA, C, LDC)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      CHARACTER TRANSA, TRANSB
      LOGICAL   NOTA, NOTB, CONJA, CONJB, AZERO, AONE,
     $          BZERO, BONE
      INTEGER*4 M,N,K, LDA,LDB,LDC, I,J,L
      REAL*8    ALPHA, BETA
      REAL*8    A(4,LDA,*), B(4,LDB,*), C(4,LDC,*)
C
      REAL*8    HTMPS, HTMPI, HTMPJ, HTMPK
C
      REAL*8    ONE, ZERO
      PARAMETER (ONE=1.0D+0, ZERO=0.0D+0)
C

      NOTA  = LSAME(TRANSA,'N')
      NOTB  = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')

      AZERO = ALPHA.EQ.ZERO
      AONE  = ALPHA.EQ.ONE
      BZERO = BETA.EQ.ZERO
      BONE  = BETA.EQ.ONE

C
C     Quick Return
C
      If ( (M.EQ.0) .OR. (N.EQ.0) ) RETURN
      If ( (AZERO .OR. (K.EQ.0)) .AND. BONE ) RETURN

C
C     ALPHA.EQ.ZERO
C
      IF ( AZERO ) THEN
        IF ( BZERO ) THEN
          DO 10 J = 1,N
          DO 10 I = 1,M
            C(1,I,J) = ZERO
            C(2,I,J) = ZERO
            C(3,I,J) = ZERO
            C(4,I,J) = ZERO
  10      CONTINUE
        ELSE
          DO 20 J = 1,N
          DO 20 I = 1,M
            C(1,I,J) = BETA * C(1,I,J) 
            C(2,I,J) = BETA * C(2,I,J) 
            C(3,I,J) = BETA * C(3,I,J) 
            C(4,I,J) = BETA * C(4,I,J) 
  20      CONTINUE
        ENDIF
        RETURN ! Nothing more to do
      ENDIF
C
      IF ( NOTB ) THEN
        IF ( NOTA ) THEN
C
C       C_{IJ} = ALPHA * A_{IK} * B_{KJ} + BETA * C_{IJ}
C
          DO 30 J = 1,N
            IF ( BZERO ) THEN
              DO 40 I = 1,M
                C(1,I,J) = ZERO
                C(2,I,J) = ZERO
                C(3,I,J) = ZERO
                C(4,I,J) = ZERO
  40          CONTINUE
            ELSEIF ( .NOT. BONE ) THEN
              DO 50 I = 1,M
                C(1,I,J) = BETA * C(1,I,J) 
                C(2,I,J) = BETA * C(2,I,J) 
                C(3,I,J) = BETA * C(3,I,J) 
                C(4,I,J) = BETA * C(4,I,J) 
  50          CONTINUE
            ENDIF
C
            DO 60 L = 1,K
              HTMPS = ALPHA * B(1,L,J)
              HTMPI = ALPHA * B(2,L,J)
              HTMPJ = ALPHA * B(3,L,J)
              HTMPK = ALPHA * B(4,L,J)
              DO 70 I = 1,M
C               Hamilton Product
                C(1,I,J) = C(1,I,J) + 
     $            A(1,I,L)*HTMPS - A(2,I,L)*HTMPI - 
     $            A(3,I,L)*HTMPJ - A(4,I,L)*HTMPK

                C(2,I,J) = C(2,I,J) + 
     $            A(1,I,L)*HTMPI + A(2,I,L)*HTMPS +
     $            A(3,I,L)*HTMPK - A(4,I,L)*HTMPJ

                C(3,I,J) = C(3,I,J) + 
     $            A(1,I,L)*HTMPJ - A(2,I,L)*HTMPK +
     $            A(3,I,L)*HTMPS + A(4,I,L)*HTMPI

                C(4,I,J) = C(4,I,J) + 
     $            A(1,I,L)*HTMPK + A(2,I,L)*HTMPJ -
     $            A(3,I,L)*HTMPI + A(4,I,L)*HTMPS
  70          CONTINUE
  60        CONTINUE
  30      CONTINUE
        ELSE IF ( CONJA ) THEN
C
C         C_{IJ} = ALPHA * CONJ(A_{KI}) * B_{KJ} + BETA * C_{IJ}
C
          DO 80 J = 1,N
          DO 80 I = 1,M
            HTMPS = ZERO 
            HTMPI = ZERO 
            HTMPJ = ZERO 
            HTMPK = ZERO 
            DO 90 L = 1,K
C             Hamilton Product (IMPLIED CONJ(A))
              HTMPS = HTMPS + A(1,L,I)*B(1,L,J) + A(2,L,I)*B(2,L,J) +
     $                        A(3,L,I)*B(3,L,J) + A(4,L,I)*B(4,L,J)

              HTMPI = HTMPI + A(1,L,I)*B(2,L,J) - A(2,L,I)*B(1,L,J) -
     $                        A(3,L,I)*B(4,L,J) + A(4,L,I)*B(3,L,J)

              HTMPJ = HTMPJ + A(1,L,I)*B(3,L,J) + A(2,L,I)*B(4,L,J) -
     $                        A(3,L,I)*B(1,L,J) - A(4,L,I)*B(2,L,J)

              HTMPK = HTMPK + A(1,L,I)*B(4,L,J) - A(2,L,I)*B(3,L,J) +
     $                        A(3,L,I)*B(2,L,J) - A(4,L,I)*B(1,L,J)
  90        CONTINUE
C
            IF( BZERO ) THEN
              C(1,I,J) = ALPHA * HTMPS
              C(2,I,J) = ALPHA * HTMPI
              C(3,I,J) = ALPHA * HTMPJ
              C(4,I,J) = ALPHA * HTMPK
            ELSE
              C(1,I,J) = ALPHA * HTMPS + BETA * C(1,I,J)
              C(2,I,J) = ALPHA * HTMPI + BETA * C(2,I,J)
              C(3,I,J) = ALPHA * HTMPJ + BETA * C(3,I,J)
              C(4,I,J) = ALPHA * HTMPK + BETA * C(4,I,J)
            ENDIF
  80      CONTINUE
        ELSE
C         IMPLIES TRANSA == 'T'
C
C         C_{IJ} = ALPHA * A_{KI} * B_{KJ} + BETA * C_{IJ}
C
          DO 100 J = 1,N
          DO 100 I = 1,M
            HTMPS = ZERO 
            HTMPI = ZERO 
            HTMPJ = ZERO 
            HTMPK = ZERO 
            DO 110 L = 1,K
C             Hamilton Product
              HTMPS = HTMPS + A(1,L,I)*B(1,L,J) - A(2,L,I)*B(2,L,J) -
     $                        A(3,L,I)*B(3,L,J) - A(4,L,I)*B(4,L,J)

              HTMPI = HTMPI + A(1,L,I)*B(2,L,J) + A(2,L,I)*B(1,L,J) +
     $                        A(3,L,I)*B(4,L,J) - A(4,L,I)*B(3,L,J)

              HTMPJ = HTMPJ + A(1,L,I)*B(3,L,J) - A(2,L,I)*B(4,L,J) +
     $                        A(3,L,I)*B(1,L,J) + A(4,L,I)*B(2,L,J)

              HTMPK = HTMPK + A(1,L,I)*B(4,L,J) + A(2,L,I)*B(3,L,J) -
     $                        A(3,L,I)*B(2,L,J) + A(4,L,I)*B(1,L,J)
  110       CONTINUE

            IF( BZERO ) THEN
              C(1,I,J) = ALPHA * HTMPS
              C(2,I,J) = ALPHA * HTMPI
              C(3,I,J) = ALPHA * HTMPJ
              C(4,I,J) = ALPHA * HTMPK
            ELSE
              C(1,I,J) = ALPHA * HTMPS + BETA * C(1,I,J)
              C(2,I,J) = ALPHA * HTMPI + BETA * C(2,I,J)
              C(3,I,J) = ALPHA * HTMPJ + BETA * C(3,I,J)
              C(4,I,J) = ALPHA * HTMPK + BETA * C(4,I,J)
            ENDIF
  100     CONTINUE
        ENDIF
      ELSE IF ( NOTA ) THEN
        IF ( CONJB ) THEN
C
C         C_{IJ} = ALPHA * A_{IK} * CONJ(B_{JK}) + BETA * C_{IJ}
C
          DO 120 J = 1,N
            IF ( BZERO ) THEN
              DO 130 I = 1,M
                C(1,I,J) = ZERO
                C(2,I,J) = ZERO
                C(3,I,J) = ZERO
                C(4,I,J) = ZERO
  130         CONTINUE
            ELSEIF ( .NOT. BONE ) THEN
              DO 140 I = 1,M
                C(1,I,J) = BETA * C(1,I,J) 
                C(2,I,J) = BETA * C(2,I,J) 
                C(3,I,J) = BETA * C(3,I,J) 
                C(4,I,J) = BETA * C(4,I,J) 
  140         CONTINUE
            ENDIF
C
            DO 150 L = 1,K
              HTMPS =  ALPHA * B(1,J,L)
              HTMPI = -ALPHA * B(2,J,L)
              HTMPJ = -ALPHA * B(3,J,L)
              HTMPK = -ALPHA * B(4,J,L)
              DO 160 I = 1,M
C               Hamilton Product
                C(1,I,J) = C(1,I,J) + 
     $            A(1,I,L)*HTMPS - A(2,I,L)*HTMPI -
     $            A(3,I,L)*HTMPJ - A(4,I,L)*HTMPK

                C(2,I,J) = C(2,I,J) + 
     $            A(1,I,L)*HTMPI + A(2,I,L)*HTMPS +
     $            A(3,I,L)*HTMPK - A(4,I,L)*HTMPJ

                C(3,I,J) = C(3,I,J) + 
     $            A(1,I,L)*HTMPJ - A(2,I,L)*HTMPK +
     $            A(3,I,L)*HTMPS + A(4,I,L)*HTMPI

                C(4,I,J) = C(4,I,J) + 
     $            A(1,I,L)*HTMPK + A(2,I,L)*HTMPJ -
     $            A(3,I,L)*HTMPI + A(4,I,L)*HTMPS
  160         CONTINUE
  150       CONTINUE
  120     CONTINUE
        ELSE
C         Implies TRANSB == 'T'
C
C         C_{IJ} = ALPHA * A_{IK} * B_{JK} + BETA * C_{IJ}
C
          DO 170 J = 1,N
            IF ( BZERO ) THEN
              DO 180 I = 1,M
                C(1,I,J) = ZERO
                C(2,I,J) = ZERO
                C(3,I,J) = ZERO
                C(4,I,J) = ZERO
  180         CONTINUE
            ELSEIF ( .NOT. BONE ) THEN
              DO 190 I = 1,M
                C(1,I,J) = BETA * C(1,I,J) 
                C(2,I,J) = BETA * C(2,I,J) 
                C(3,I,J) = BETA * C(3,I,J) 
                C(4,I,J) = BETA * C(4,I,J) 
  190         CONTINUE
            ENDIF
C
            DO 200 L = 1,K
              HTMPS = ALPHA * B(1,J,L)
              HTMPI = ALPHA * B(2,J,L)
              HTMPJ = ALPHA * B(3,J,L)
              HTMPK = ALPHA * B(4,J,L)
              DO 210 I = 1,M
C               Hamilton Product
                C(1,I,J) = C(1,I,J) + 
     $            A(1,I,L)*HTMPS - A(2,I,L)*HTMPI -
     $            A(3,I,L)*HTMPJ - A(4,I,L)*HTMPK

                C(2,I,J) = C(2,I,J) + 
     $            A(1,I,L)*HTMPI + A(2,I,L)*HTMPS +
     $            A(3,I,L)*HTMPK - A(4,I,L)*HTMPJ

                C(3,I,J) = C(3,I,J) + 
     $            A(1,I,L)*HTMPJ - A(2,I,L)*HTMPK +
     $            A(3,I,L)*HTMPS + A(4,I,L)*HTMPI

                C(4,I,J) = C(4,I,J) + 
     $            A(1,I,L)*HTMPK + A(2,I,L)*HTMPJ -
     $            A(3,I,L)*HTMPI + A(4,I,L)*HTMPS
  210         CONTINUE
  200       CONTINUE
  170     CONTINUE
        ENDIF
      ELSE IF ( CONJA ) THEN
        IF ( CONJB ) THEN
C
C         C_{IJ} = ALPHA * CONJ(A_{KI}) * CONJ(B_{JK}) + BETA * C_{IJ}
C
          DO 220 J = 1,N
          DO 220 I = 1,M
            HTMPS = ZERO 
            HTMPI = ZERO 
            HTMPJ = ZERO 
            HTMPK = ZERO 
            DO 230 L = 1,K
C             Hamilton Product (IMPLIED CONJ(A) + CONJ(B))
              HTMPS = HTMPS + A(1,L,I)*B(1,J,L) - A(2,L,I)*B(2,J,L) -
     $                        A(3,L,I)*B(3,J,L) - A(4,L,I)*B(4,J,L)

              HTMPI = HTMPI - A(1,L,I)*B(2,J,L) - A(2,L,I)*B(1,J,L) +
     $                        A(3,L,I)*B(4,J,L) - A(4,L,I)*B(3,J,L)

              HTMPJ = HTMPJ - A(1,L,I)*B(3,J,L) - A(2,L,I)*B(4,J,L) -
     $                        A(3,L,I)*B(1,J,L) + A(4,L,I)*B(2,J,L)

              HTMPK = HTMPK - A(1,L,I)*B(4,J,L) + A(2,L,I)*B(3,J,L) -
     $                        A(3,L,I)*B(2,J,L) - A(4,L,I)*B(1,J,L)
  230       CONTINUE
C
            IF( BZERO ) THEN
              C(1,I,J) = ALPHA * HTMPS
              C(2,I,J) = ALPHA * HTMPI
              C(3,I,J) = ALPHA * HTMPJ
              C(4,I,J) = ALPHA * HTMPK
            ELSE
              C(1,I,J) = ALPHA * HTMPS + BETA * C(1,I,J)
              C(2,I,J) = ALPHA * HTMPI + BETA * C(2,I,J)
              C(3,I,J) = ALPHA * HTMPJ + BETA * C(3,I,J)
              C(4,I,J) = ALPHA * HTMPK + BETA * C(4,I,J)
            ENDIF
  220     CONTINUE
        ELSE
C         Implies TRANSB == 'T'
C
C         C_{IJ} = ALPHA * CONJ(A_{KI}) * B_{JK} + BETA * C_{IJ}
C
          DO 240 J = 1,N
          DO 240 I = 1,M
            HTMPS = ZERO 
            HTMPI = ZERO 
            HTMPJ = ZERO 
            HTMPK = ZERO 
            DO 250 L = 1,K
C             Hamilton Product (IMPLIED CONJ(A))
              HTMPS = HTMPS + A(1,L,I)*B(1,J,L) + A(2,L,I)*B(2,J,L) +
     $                        A(3,L,I)*B(3,J,L) + A(4,L,I)*B(4,J,L)
                                               
              HTMPI = HTMPI + A(1,L,I)*B(2,J,L) - A(2,L,I)*B(1,J,L) -
     $                        A(3,L,I)*B(4,J,L) + A(4,L,I)*B(3,J,L)
                                               
              HTMPJ = HTMPJ + A(1,L,I)*B(3,J,L) + A(2,L,I)*B(4,J,L) -
     $                        A(3,L,I)*B(1,J,L) - A(4,L,I)*B(2,J,L)
                                               
              HTMPK = HTMPK + A(1,L,I)*B(4,J,L) - A(2,L,I)*B(3,J,L) +
     $                        A(3,L,I)*B(2,J,L) - A(4,L,I)*B(1,J,L)
  250       CONTINUE
C
            IF( BZERO ) THEN
              C(1,I,J) = ALPHA * HTMPS
              C(2,I,J) = ALPHA * HTMPI
              C(3,I,J) = ALPHA * HTMPJ
              C(4,I,J) = ALPHA * HTMPK
            ELSE
              C(1,I,J) = ALPHA * HTMPS + BETA * C(1,I,J)
              C(2,I,J) = ALPHA * HTMPI + BETA * C(2,I,J)
              C(3,I,J) = ALPHA * HTMPJ + BETA * C(3,I,J)
              C(4,I,J) = ALPHA * HTMPK + BETA * C(4,I,J)
            ENDIF
  240     CONTINUE
        ENDIF
      ELSE
C       Implies TRANSA == 'T'
        IF ( CONJB ) THEN
C
C         C_{IJ} = ALPHA * A_{KI} * CONJ(B_{JK}) + BETA * C_{IJ}
C
          DO 260 J = 1,N
          DO 260 I = 1,M
            HTMPS = ZERO 
            HTMPI = ZERO 
            HTMPJ = ZERO 
            HTMPK = ZERO 
            DO 270 L = 1,K
C             Hamilton Product
              HTMPS = HTMPS + A(1,L,I)*B(1,J,L) + A(2,L,I)*B(2,J,L) +
     $                        A(3,L,I)*B(3,J,L) + A(4,L,I)*B(4,J,L)

              HTMPI = HTMPI - A(1,L,I)*B(2,J,L) + A(2,L,I)*B(1,J,L) -
     $                        A(3,L,I)*B(4,J,L) + A(4,L,I)*B(3,J,L)

              HTMPJ = HTMPJ - A(1,L,I)*B(3,J,L) + A(2,L,I)*B(4,J,L) +
     $                        A(3,L,I)*B(1,J,L) - A(4,L,I)*B(2,J,L)

              HTMPK = HTMPK - A(1,L,I)*B(4,J,L) - A(2,L,I)*B(3,J,L) +
     $                        A(3,L,I)*B(2,J,L) + A(4,L,I)*B(1,J,L)
  270       CONTINUE

            IF( BZERO ) THEN
              C(1,I,J) = ALPHA * HTMPS
              C(2,I,J) = ALPHA * HTMPI
              C(3,I,J) = ALPHA * HTMPJ
              C(4,I,J) = ALPHA * HTMPK
            ELSE
              C(1,I,J) = ALPHA * HTMPS + BETA * C(1,I,J)
              C(2,I,J) = ALPHA * HTMPI + BETA * C(2,I,J)
              C(3,I,J) = ALPHA * HTMPJ + BETA * C(3,I,J)
              C(4,I,J) = ALPHA * HTMPK + BETA * C(4,I,J)
            ENDIF
  260     CONTINUE
        ELSE
C
C         C_{IJ} = ALPHA * A_{KI} * B_{JK} + BETA * C_{IJ}
C
          DO 280 J = 1,N
          DO 280 I = 1,M
            HTMPS = ZERO 
            HTMPI = ZERO 
            HTMPJ = ZERO 
            HTMPK = ZERO 
            DO 290 L = 1,K
C             Hamilton Product
              HTMPS = HTMPS + A(1,L,I)*B(1,J,L) - A(2,L,I)*B(2,J,L) -
     $                        A(3,L,I)*B(3,J,L) - A(4,L,I)*B(4,J,L)

              HTMPI = HTMPI + A(1,L,I)*B(2,J,L) + A(2,L,I)*B(1,J,L) +
     $                        A(3,L,I)*B(4,J,L) - A(4,L,I)*B(3,J,L)

              HTMPJ = HTMPJ + A(1,L,I)*B(3,J,L) - A(2,L,I)*B(4,J,L) +
     $                        A(3,L,I)*B(1,J,L) + A(4,L,I)*B(2,J,L)

              HTMPK = HTMPK + A(1,L,I)*B(4,J,L) + A(2,L,I)*B(3,J,L) -
     $                        A(3,L,I)*B(2,J,L) + A(4,L,I)*B(1,J,L)
  290       CONTINUE

            IF( BZERO ) THEN
              C(1,I,J) = ALPHA * HTMPS
              C(2,I,J) = ALPHA * HTMPI
              C(3,I,J) = ALPHA * HTMPJ
              C(4,I,J) = ALPHA * HTMPK
            ELSE
              C(1,I,J) = ALPHA * HTMPS + BETA * C(1,I,J)
              C(2,I,J) = ALPHA * HTMPI + BETA * C(2,I,J)
              C(3,I,J) = ALPHA * HTMPJ + BETA * C(3,I,J)
              C(4,I,J) = ALPHA * HTMPK + BETA * C(4,I,J)
            ENDIF
  280     CONTINUE
        ENDIF
      ENDIF
      END
