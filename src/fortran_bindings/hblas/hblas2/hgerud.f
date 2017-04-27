      Subroutine HGERUD(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
C      
      LOGICAL    AZERO,AONE
      INTEGER*4  M,N, LDA,INCX,INCY, I,J,IY,JX,IX,JY, KX
      REAL*8     ALPHA
      COMPLEX*16 X(2,*), Y(2,*), A(2,LDA,*)
      COMPLEX*16    HTMPS, HTMPJ
C
      REAL*8        ONED, ZEROD
      COMPLEX*16    ONEZ, ZEROZ
      PARAMETER (ONED=1.0D+0, ZEROD=0.0D+0)
      PARAMETER (ONEZ=(1.0D+0,0.0D+0), ZEROZ=(0.0D+0,0.0D+0))
C

      AZERO = ALPHA.EQ.ZEROD
      AONE  = ALPHA.EQ.ONED

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
            HTMPS = ALPHA * Y(1,JY)
            HTMPJ = ALPHA * Y(2,JY)
            DO 10 I = 1,M
              A(1,I,J) = A(1,I,J) + X(1,I)*HTMPS - X(2,I)*CONJG(HTMPJ)
              A(2,I,J) = A(2,I,J) + X(1,I)*HTMPJ + X(2,I)*CONJG(HTMPS)
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
            HTMPI = ALPHA * Y(1,JY)
            HTMPJ = ALPHA * Y(2,JY)
            IX = KX
            DO 30 I = 1,M
              A(1,I,J) = A(1,I,J) + X(1,IX)*HTMPS - X(2,IX)*CONJG(HTMPJ)
              A(2,I,J) = A(2,I,J) + X(1,IX)*HTMPJ + X(2,IX)*CONJG(HTMPS)
              IX = IX + INCX
  30        CONTINUE
          ENDIF
          JY = JY + INCY
  40    CONTINUE
      ENDIF 
C     
      END
