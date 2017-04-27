      Subroutine HGEMVHD(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA,
     $  Y, INCY)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      CHARACTER TRANS
      LOGICAL   NOCONJ, AZERO,AONE,BZERO,BONE
      INTEGER*4 M,N, LDA,INCX,INCY, LENX,LENY, I,J,IY,JX,IX,JY, KX,KY
      REAL*8    BETA
      COMPLEX*16    A(2,LDA,*), X(2,*), Y(2,*), ALPHA(2)
C
      COMPLEX*16    HTMP1S, HTMP1J, HTMP2S, HTMP2J
      COMPLEX*16    HTMP3S, HTMP3J, HTMP4S, HTMP4J
C
      COMPLEX*16    ONE, ZERO
      PARAMETER (ONE=(1.0D+0,0.0D+0), ZERO=(0.0D+0,0.0D+0))
C
      NOCONJ = LSAME('T',TRANS)
      IF ( LSAME('T',TRANS) ) THEN
        LENX = N
        LENY = M
      ELSE
        LENX = M
        LENY = N
      ENDIF

      AZERO = (ALPHA(1).EQ.ZERO).AND.((ALPHA(2).EQ.ZERO))
      AONE  = (ALPHA(1).EQ.ONE).AND.((ALPHA(2).EQ.ZERO))
      BZERO = BETA.EQ.ZERO
      BONE  = BETA.EQ.ONE
C
      IF ( (M.EQ.0) .OR. (N.EQ.0) .OR. 
     $     ( AZERO .AND. BONE ) ) RETURN
C
C
      IF ( INCX.GT.0 ) THEN
        KX = 1
      ELSE
        KX = 1 - (LENX-1)*INCX
      ENDIF

      IF ( INCY.GT.0 ) THEN
        KY = 1
      ELSE
        KY = 1 - (LENY-1)*INCY
      ENDIF
C

      IF ( .NOT.BONE ) THEN
        IF ( INCY.EQ.1 ) THEN
          IF ( BZERO ) THEN
            DO 10 I = 1,LENY
              Y(1,I) = ZERO
              Y(2,I) = ZERO
  10        CONTINUE
          ELSE
            DO 20 I = 1,LENY
              Y(1,I) = BETA * Y(1,I) 
              Y(2,I) = BETA * Y(2,I) 
  20        CONTINUE
          ENDIF
        ELSE
          IY = KY
          IF ( BZERO ) THEN
            DO 30 I = 1,LENY
              Y(1,IY) = ZERO
              Y(2,IY) = ZERO
              IY = IY + INCY
  30        CONTINUE
          ELSE
            DO 40 I = 1,LENY
              Y(1,IY) = BETA * Y(1,IY) 
              Y(2,IY) = BETA * Y(2,IY) 
              IY = IY + INCY
  40        CONTINUE
          ENDIF
        ENDIF
      ENDIF
C
      IF ( AZERO )  RETURN ! Nothing left to do
C
      IF ( LSAME(TRANS,'N') ) THEN
        JX = KX
        IF ( INCY.EQ.1 ) THEN
          DO 60 J = 1,N
            HTMP1S = ALPHA(1) * X(1,JX)
            HTMP1J = ALPHA(1) * X(2,JX)
            HTMP2S = ALPHA(1) * CONJG(X(1,JX))
            HTMP2J = ALPHA(1) * CONJG(X(2,JX))
            HTMP3S = ALPHA(2) * X(1,JX)
            HTMP3J = ALPHA(2) * X(2,JX)
            HTMP4S = ALPHA(2) * CONJG(X(1,JX))
            HTMP4J = ALPHA(2) * CONJG(X(2,JX))
            Do 50 I = 1,M
C             Hamilton Product
              Y(1,I) = Y(1,I) + 
     $          A(1,I,J)*HTMP1S - A(2,I,J)*HTMP2J -
     $          CONJG(A(1,I,J))*HTMP4J - CONJG(A(2,I,J))*HTMP3S
              Y(2,I) = Y(2,I) + 
     $          A(1,I,J)*HTMP1J + A(2,I,J)*HTMP2S +
     $          CONJG(A(1,I,J))*HTMP4S - CONJG(A(2,I,J))*HTMP3J
  50        CONTINUE 
            JX = JX + INCX
  60      CONTINUE
        ELSE
          DO 80 J = 1,N
            HTMP1S = ALPHA(1) * X(1,JX)
            HTMP1J = ALPHA(1) * X(2,JX)
            HTMP2S = ALPHA(1) * CONJG(X(1,JX))
            HTMP2J = ALPHA(1) * CONJG(X(2,JX))
            HTMP3S = ALPHA(2) * X(1,JX)
            HTMP3J = ALPHA(2) * X(2,JX)
            HTMP4S = ALPHA(2) * CONJG(X(1,JX))
            HTMP4J = ALPHA(2) * CONJG(X(2,JX))
            IY = KY
            Do 70 I = 1,M
C             Hamilton Product
              Y(1,IY) = Y(1,IY) + 
     $          A(1,I,J)*HTMP1S - A(2,I,J)*HTMP2J -
     $          CONJG(A(1,I,J))*HTMP4J - CONJG(A(2,I,J))*HTMP3S
              Y(2,IY) = Y(2,IY) + 
     $          A(1,I,J)*HTMP1J + A(2,I,J)*CONJG(HTMP2S) +
     $          CONJG(A(1,I,J))*HTMP4S - CONJG(A(2,I,J))*HTMP3J
              IY = IY + INCY
  70        CONTINUE 
            JX = JX + INCX
  80      CONTINUE
        ENDIF
      ELSE
        JY = KY
        IF ( INCX.EQ.1 ) THEN
          DO 110 J = 1,N
            HTMP1S = ZERO
            HTMP1J = ZERO
            IF ( NOCONJ ) THEN
              DO 90 I = 1,M
C             Hamilton Product
                HTMP1S = HTMP1S + 
     $            A(1,I,J)*X(1,I) - A(2,I,J)*CONJG(X(2,I))
                HTMP1J = HTMP1J + 
     $            A(1,I,J)*X(2,I) + A(2,I,J)*CONJG(X(1,I))
  90          CONTINUE
            ELSE
              DO 100 I = 1,M
C               Hamilton Product (IMPLIED CONJ(A))
                HTMP1S = HTMP1S + 
     $            CONJG(A(1,I,J))*X(1,I) + A(2,I,J)*CONJG(X(2,I))
                HTMP1J = HTMP1J + 
     $            CONJG(A(1,I,J))*X(2,I) - A(2,I,J)*CONJG(X(1,I))
  100         CONTINUE
            ENDIF
            Y(1,JY) = Y(1,JY) + ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
            Y(2,JY) = Y(2,JY) + ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            JY = JY + INCY
  110     CONTINUE
        ELSE
          DO 140 J = 1,N
            HTMP1S = ZERO
            HTMP1J = ZERO
            IX = KX
            IF ( NOCONJ ) THEN
              DO 120 I = 1,M
C             Hamilton Product
                HTMP1S = HTMP1S + 
     $            A(1,I,J)*X(1,IX) - A(2,I,J)*CONJG(X(2,IX))
                HTMP1J = HTMP1J + 
     $            A(1,I,J)*X(2,IX) + A(2,I,J)*CONJG(X(1,IX))
                IX = IX + INCX
  120         CONTINUE
            ELSE
              DO 130 I = 1,M
C               Hamilton Product (IMPLIED CONJ(A))
                HTMP1S = HTMP1S + 
     $            CONJG(A(1,I,J))*X(1,IX) + A(2,I,J)*CONJG(X(2,IX))
                HTMP1J = HTMP1J + 
     $            CONJG(A(1,I,J))*X(2,IX) - A(2,I,J)*CONJG(X(1,IX))
                IX = IX + INCX
  130         CONTINUE
            ENDIF
            Y(1,JY) = Y(1,JY) + ALPHA(1)*HTMP1S - ALPHA(2)*CONJG(HTMP1J)
            Y(2,JY) = Y(2,JY) + ALPHA(1)*HTMP1J + ALPHA(2)*CONJG(HTMP1S)
            JY = JY + INCY
  140     CONTINUE
        ENDIF
      ENDIF
C
      END
