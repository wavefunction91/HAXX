      Subroutine HAXPYDH(SIDE, N, ALPHA, X, INCX, Y, INCY)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      CHARACTER SIDE
      INTEGER*4 N, INCX, INCY, I, IX, IY
      REAL*8    ALPHA, X(4,*), Y(4,*)
C
      IF ( N.LE.0 ) RETURN
      IF ( ALPHA.EQ.0.0D0 ) RETURN
C
      IF ( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
        DO 10 I = 1,N
          Y(1,I) = Y(1,I) + ALPHA*X(1,I)
          Y(2,I) = Y(2,I) + ALPHA*X(2,I)
          Y(3,I) = Y(3,I) + ALPHA*X(3,I)
          Y(4,I) = Y(4,I) + ALPHA*X(4,I)
  10    CONTINUE
      ELSE

        IX = 1
        IY = 1
        IF ( INCX.LT.0 ) IX = (-N + 1)*INCX + 1
        IF ( INCY.LT.0 ) IY = (-N + 1)*INCY + 1

        DO 20 I = 1,N
          Y(1,IY) = Y(1,IY) + ALPHA*X(1,IX)
          Y(2,IY) = Y(2,IY) + ALPHA*X(2,IX)
          Y(3,IY) = Y(3,IY) + ALPHA*X(3,IX)
          Y(4,IY) = Y(4,IY) + ALPHA*X(4,IX)
          IX = IX + INCX
          IY = IY + INCY
  20    CONTINUE
      ENDIF
C
      END
