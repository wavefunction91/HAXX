C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HGEMMHZ(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
     $  B, LDB, BETA, C, LDC)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      CHARACTER TRANSA, TRANSB
      LOGICAL   NOTA, NOTB, CONJA, CONJB, AZERO, AONE,
     $          BZERO, BONE
      INTEGER*4 M,N,K, LDA,LDB,LDC, I,J,L
      Complex*16    A(2,LDA,*), B(2,LDB,*), C(2,LDC,*), ALPHA(2), BETA
C
      COMPLEX*16    HTMP1S, HTMP1J, HTMP2S, HTMP2J
      COMPLEX*16    HTMP3S, HTMP3J, HTMP4S, HTMP4J
C
      COMPLEX*16    ONE, ZERO
      PARAMETER (ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0))
C

      NOTA  = LSAME(TRANSA,'N')
      NOTB  = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')

      AZERO = (ALPHA(1).EQ.ZERO).AND.((ALPHA(2).EQ.ZERO))
      AONE  = (ALPHA(1).EQ.ONE).AND.((ALPHA(2).EQ.ZERO))
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
  10      CONTINUE
        ELSE
          DO 20 J = 1,N
          DO 20 I = 1,M
            C(1,I,J) = BETA * C(1,I,J) 
            C(2,I,J) = BETA * C(2,I,J) 
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
  40          CONTINUE
            ELSEIF ( .NOT. BONE ) THEN
              DO 50 I = 1,M
                C(1,I,J) = BETA * C(1,I,J) 
                C(2,I,J) = BETA * C(2,I,J) 
  50          CONTINUE
            ENDIF
C
            DO 60 L = 1,K
              HTMP1S = ALPHA(1) * B(1,L,J) ! AE
              HTMP1J = ALPHA(1) * B(2,L,J) ! AF
              HTMP2S = ALPHA(1) * CONJG(B(1,L,J)) ! A CONJ(E)
              HTMP2J = ALPHA(1) * CONJG(B(2,L,J)) ! A CONJ(F)
              HTMP3S = ALPHA(2) * B(1,L,J) ! BE
              HTMP3J = ALPHA(2) * B(2,L,J) ! BF
              HTMP4S = ALPHA(2) * CONJG(B(1,L,J)) ! B CONJ(E)
              HTMP4J = ALPHA(2) * CONJG(B(2,L,J)) ! B CONJ(F)
              DO 70 I = 1,M
C               Hamilton Product
                C(1,I,J) = C(1,I,J) + 
     $            A(1,I,L)*HTMP1S - A(2,I,L)*HTMP2J - 
     $            CONJG(A(1,I,L))*HTMP4J - CONJG(A(2,I,L))*HTMP3S 
                C(2,I,J) = C(2,I,J) + 
     $            A(1,I,L)*HTMP1J + A(2,I,L)*HTMP2S +
     $            CONJG(A(1,I,L))*HTMP4S - CONJG(A(2,I,L))*HTMP3J
  70          CONTINUE
  60        CONTINUE
  30      CONTINUE
        ELSE IF ( CONJA ) THEN
C
C         C_{IJ} = ALPHA * CONJ(A_{KI}) * B_{KJ} + BETA * C_{IJ}
C
          DO 80 J = 1,N
          DO 80 I = 1,M
            HTMP1S = ZERO 
            HTMP1J = ZERO 
            DO 90 L = 1,K
C             Hamilton Product (IMPLIED CONJ(A))
              HTMP1S = HTMP1S + 
     $          CONJG(A(1,L,I))*B(1,L,J) + A(2,L,I)*CONJG(B(2,L,J))
              HTMP1J = HTMP1J + 
     $          CONJG(A(1,L,I))*B(2,L,J) - A(2,L,I)*CONJG(B(1,L,J))
  90        CONTINUE
C
            IF( BZERO ) THEN
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            ELSE
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J) + 
     $          BETA * C(1,I,J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S) + 
     $          BETA * C(2,I,J)
            ENDIF
  80      CONTINUE
        ELSE
C         IMPLIES TRANSA == 'T'
C
C         C_{IJ} = ALPHA * A_{KI} * B_{KJ} + BETA * C_{IJ}
C
          DO 100 J = 1,N
          DO 100 I = 1,M
            HTMP1S = ZERO 
            HTMP1J = ZERO 
            DO 110 L = 1,K
C             Hamilton Product
              HTMP1S = HTMP1S + 
     $          A(1,L,I)*B(1,L,J) - A(2,L,I)*CONJG(B(2,L,J))
              HTMP1J = HTMP1J + 
     $          A(1,L,I)*B(2,L,J) + A(2,L,I)*CONJG(B(1,L,J))
  110       CONTINUE

            IF( BZERO ) THEN
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            ELSE
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J) + 
     $          BETA * C(1,I,J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S) + 
     $          BETA * C(2,I,J)
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
  130         CONTINUE
            ELSEIF ( .NOT. BONE ) THEN
              DO 140 I = 1,M
                C(1,I,J) = BETA * C(1,I,J) 
                C(2,I,J) = BETA * C(2,I,J) 
  140         CONTINUE
            ENDIF
C
            DO 150 L = 1,K
              HTMP1S = ALPHA(1) * B(1,J,L) ! AE
              HTMP1J = -ALPHA(1) * B(2,J,L) ! AF
              HTMP2S = ALPHA(1) * CONJG(B(1,J,L)) ! A CONJ(E)
              HTMP2J = -ALPHA(1) * CONJG(B(2,J,L)) ! A CONJ(F)
              HTMP3S = ALPHA(2) * B(1,J,L) ! BE
              HTMP3J = -ALPHA(2) * B(2,J,L) ! BF
              HTMP4S = ALPHA(2) * CONJG(B(1,J,L)) ! B CONJ(E)
              HTMP4J = -ALPHA(2) * CONJG(B(2,J,L)) ! B CONJ(F)
              DO 160 I = 1,M
C               Hamilton Product (IMPLIES CONJB)
                C(1,I,J) = C(1,I,J) + 
     $            A(1,I,L)*HTMP2S - A(2,I,L)*HTMP2J - 
     $            CONJG(A(1,I,L))*HTMP4J - CONJG(A(2,I,L))*HTMP4S 
                C(2,I,J) = C(2,I,J) + 
     $            A(1,I,L)*HTMP1J + A(2,I,L)*HTMP1S +
     $            CONJG(A(1,I,L))*HTMP3S - CONJG(A(2,I,L))*HTMP3J
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
  180         CONTINUE
            ELSEIF ( .NOT. BONE ) THEN
              DO 190 I = 1,M
                C(1,I,J) = BETA * C(1,I,J) 
                C(2,I,J) = BETA * C(2,I,J) 
  190         CONTINUE
            ENDIF
C
            DO 200 L = 1,K
              HTMP1S = ALPHA(1) * B(1,J,L) ! AE
              HTMP1J = ALPHA(1) * B(2,J,L) ! AF
              HTMP2S = ALPHA(1) * CONJG(B(1,J,L)) ! A CONJ(E)
              HTMP2J = ALPHA(1) * CONJG(B(2,J,L)) ! A CONJ(F)
              HTMP3S = ALPHA(2) * B(1,J,L) ! BE
              HTMP3J = ALPHA(2) * B(2,J,L) ! BF
              HTMP4S = ALPHA(2) * CONJG(B(1,J,L)) ! B CONJ(E)
              HTMP4J = ALPHA(2) * CONJG(B(2,J,L)) ! B CONJ(F)
              DO 210 I = 1,M
C               Hamilton Product
                C(1,I,J) = C(1,I,J) + 
     $            A(1,I,L)*HTMP1S - A(2,I,L)*HTMP2J - 
     $            CONJG(A(1,I,L))*HTMP4J - CONJG(A(2,I,L))*HTMP3S 
                C(2,I,J) = C(2,I,J) + 
     $            A(1,I,L)*HTMP1J + A(2,I,L)*HTMP2S +
     $            CONJG(A(1,I,L))*HTMP4S - CONJG(A(2,I,L))*HTMP3J
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
            HTMP1S = ZERO 
            HTMP1J = ZERO 
            DO 230 L = 1,K
C             Hamilton Product (IMPLIED CONJ(A) + CONJ(B))
              HTMP1S = HTMP1S + 
     $          CONJG(A(1,L,I))*CONJG(B(1,J,L)) - 
     $          A(2,L,I)*CONJG(B(2,J,L))
              HTMP1J = HTMP1J - 
     $          CONJG(A(1,L,I))*B(2,J,L) - A(2,L,I)*B(1,J,L)
  230       CONTINUE
C
            IF( BZERO ) THEN
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            ELSE
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J) + 
     $          BETA * C(1,I,J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S) + 
     $          BETA * C(2,I,J)
            ENDIF
  220     CONTINUE
        ELSE
C         Implies TRANSB == 'T'
C
C         C_{IJ} = ALPHA * CONJ(A_{KI}) * B_{JK} + BETA * C_{IJ}
C
          DO 240 J = 1,N
          DO 240 I = 1,M
            HTMP1S = ZERO 
            HTMP1J = ZERO 
            DO 250 L = 1,K
C             Hamilton Product (IMPLIED CONJ(A))
              HTMP1S = HTMP1S + 
     $          CONJG(A(1,L,I))*B(1,J,L) + A(2,L,I)*CONJG(B(2,J,L))
              HTMP1J = HTMP1J + 
     $          CONJG(A(1,L,I))*B(2,J,L) - A(2,L,I)*CONJG(B(1,J,L))
  250       CONTINUE
C
            IF( BZERO ) THEN
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            ELSE
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J) + 
     $          BETA * C(1,I,J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S) + 
     $          BETA * C(2,I,J)
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
            HTMP1S = ZERO 
            HTMP1J = ZERO 
            DO 270 L = 1,K
C             Hamilton Product (IMPLIED CONJ(B))
              HTMP1S = HTMP1S + 
     $          A(1,L,I)*CONJG(B(1,J,L)) + A(2,L,I)*CONJG(B(2,J,L))
              HTMP1J = HTMP1J - 
     $          A(1,L,I)*B(2,J,L) + A(2,L,I)*B(1,J,L)
  270       CONTINUE

            IF( BZERO ) THEN
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            ELSE
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J) + 
     $          BETA * C(1,I,J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S) + 
     $          BETA * C(2,I,J)
            ENDIF
  260     CONTINUE
        ELSE
C
C         C_{IJ} = ALPHA * A_{KI} * B_{JK} + BETA * C_{IJ}
C
          DO 280 J = 1,N
          DO 280 I = 1,M
            HTMP1S = ZERO 
            HTMP1J = ZERO 
            DO 290 L = 1,K
C             Hamilton Product
              HTMP1S = HTMP1S + 
     $          A(1,L,I)*B(1,J,L) - A(2,L,I)*CONJG(B(2,J,L))
              HTMP1J = HTMP1J + 
     $          A(1,L,I)*B(2,J,L) + A(2,L,I)*CONJG(B(1,J,L))
  290       CONTINUE

            IF( BZERO ) THEN
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            ELSE
              C(1,I,J) = ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J) + 
     $          BETA * C(1,I,J)
              C(2,I,J) = ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S) + 
     $          BETA * C(2,I,J)
            ENDIF
  280     CONTINUE
        ENDIF
      ENDIF
      END
