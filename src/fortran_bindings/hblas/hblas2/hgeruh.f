      Subroutine HGERUH(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
C      
      LOGICAL    AZERO,AONE
      INTEGER*4  M,N, LDA,INCX,INCY, I,J,IY,JX,IX,JY, KX
      COMPLEX*16 X(2,*), Y(2,*), A(2,LDA,*), ALPHA(2)
      COMPLEX*16    HTMP1S, HTMP1J, HTMP2S, HTMP2J
      COMPLEX*16    HTMP3S, HTMP3J, HTMP4S, HTMP4J
C
      REAL*8        ONED, ZEROD
      COMPLEX*16    ONEZ, ZEROZ
      PARAMETER (ONED=1.0D+0, ZEROD=0.0D+0)
      PARAMETER (ONEZ=(1.0D+0,0.0D+0), ZEROZ=(0.0D+0,0.0D+0))
C

      AZERO = (ALPHA(1).EQ.ZEROZ).AND.(ALPHA(2).EQ.ZEROZ)
      AONE  = (ALPHA(1).EQ.ONEZ).AND.(ALPHA(2).EQ.ZEROZ)

      IF ( M.EQ.0 .OR. N.EQ.0 .OR. AZERO ) RETURN
C
      IF ( INCY.GT.0 ) THEN
        JY = 1
      ELSE
        JY = 1 - (N-1)*INCY
      ENDIF
C
      IF ( INCX.EQ.1 ) THEN
        DO 20 J = 1,N
          IF ( Y(1,JY).NE.ZEROZ .OR. Y(2,JY).NE.ZEROZ ) THEN
            HTMP1S = ALPHA(1) * Y(1,JY)
            HTMP1J = ALPHA(1) * Y(2,JY)
            HTMP2S = ALPHA(1) * CONJG(Y(1,JY))
            HTMP2J = ALPHA(1) * CONJG(Y(2,JY))
            HTMP3S = ALPHA(2) * Y(1,JY)
            HTMP3J = ALPHA(2) * Y(2,JY)
            HTMP4S = ALPHA(2) * CONJG(Y(1,JY))
            HTMP4J = ALPHA(2) * CONJG(Y(2,JY))
            DO 10 I = 1,M
              A(1,I,J) = A(1,I,J) + 
     $          X(1,I)*HTMP1S - X(2,I)*HTMP2J -
     $          CONJG(X(1,I))*HTMP4J - CONJG(X(2,I))*HTMP3S
              A(2,I,J) = A(2,I,J) + 
     $          X(1,I)*HTMP1J + X(2,I)*HTMP2S +
     $          CONJG(X(1,I))*HTMP4S - CONJG(X(2,I))*HTMP3J
  10        CONTINUE
          ENDIF
          JY = JY + INCY
  20    CONTINUE
      ELSE
        IF ( INCX.GT.0 ) THEN
          KX = 1
        ELSE
          KX = 1 - (M-1)*INCX
        ENDIF
        DO 40 J = 1,N
          IF ( Y(1,JY).NE.ZEROZ .OR. Y(2,JY).NE.ZEROZ ) THEN
            HTMP1S = ALPHA(1) * Y(1,JY)
            HTMP1J = ALPHA(1) * Y(2,JY)
            HTMP2S = ALPHA(1) * CONJG(Y(1,JY))
            HTMP2J = ALPHA(1) * CONJG(Y(2,JY))
            HTMP3S = ALPHA(2) * Y(1,JY)
            HTMP3J = ALPHA(2) * Y(2,JY)
            HTMP4S = ALPHA(2) * CONJG(Y(1,JY))
            HTMP4J = ALPHA(2) * CONJG(Y(2,JY))
            IX = KX
            DO 30 I = 1,M
              A(1,I,J) = A(1,I,J) + 
     $          X(1,IX)*HTMP1S - X(2,IX)*HTMP2J -
     $          CONJG(X(1,IX))*HTMP4J - CONJG(X(2,IX))*HTMP3S
              A(2,I,J) = A(2,I,J) + 
     $          X(1,IX)*HTMP1J + X(2,IX)*HTMP2S +
     $          CONJG(X(1,IX))*HTMP4S - CONJG(X(2,IX))*HTMP3J
              IX = IX + INCX
  30        CONTINUE
          ENDIF
          JY = JY + INCY
  40    CONTINUE
      ENDIF 
C     
      END
