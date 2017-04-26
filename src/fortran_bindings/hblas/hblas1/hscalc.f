      Subroutine HSCALC(SIDE, N, ALPHA, X, INCX)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      LOGICAL   MRIGHT
      CHARACTER SIDE
      INTEGER*4 N, INCX, I, IX
      REAL*8    ALPHA(2), X(4,*)
C
      REAL*8    HTMPS, HTMPI, HTMPJ
C
      IF ( N.LE.0 ) RETURN
      IF ( ALPHA(1).EQ.0.0D0 .AND. ALPHA(2).EQ.0.0D0 ) RETURN
C
      MRIGHT = LSAME('R',SIDE)
C
      IF ( INCX.EQ.1  ) THEN
        IF ( .NOT. MRIGHT ) THEN
          DO 10 I = 1,N
C           Hamilton Product
            HTMPS  = ALPHA(1)*X(1,I) - ALPHA(2)*X(2,I)
            HTMPI  = ALPHA(1)*X(2,I) + ALPHA(2)*X(1,I)
            HTMPJ  = ALPHA(1)*X(3,I) - ALPHA(2)*X(4,I)
            X(4,I) = ALPHA(1)*X(4,I) + ALPHA(2)*X(3,I)
            X(1,I) = HTMPS
            X(2,I) = HTMPI
            X(3,I) = HTMPJ
  10      CONTINUE
        ELSE
          DO 20 I = 1,N
C           Hamilton Product
            HTMPS   =  X(1,I)*ALPHA(1) - X(2,I)*ALPHA(2)
            HTMPI   =  X(1,I)*ALPHA(2) + X(2,I)*ALPHA(1)
            HTMPJ   =  X(3,I)*ALPHA(1) + X(4,I)*ALPHA(2)
            X(4,I)  = -X(3,I)*ALPHA(2) + X(4,I)*ALPHA(1)
            X(1,I)  = HTMPS
            X(2,I)  = HTMPI
            X(3,I)  = HTMPJ
  20      CONTINUE
        ENDIF
      ELSE

        IX = 1
        IF ( INCX.LT.0 ) IX = (-N + 1)*INCX + 1

        IF ( .NOT. MRIGHT ) THEN
          DO 30 I = 1,N
C           Hamilton Product
            HTMPS   = ALPHA(1)*X(1,IX) - ALPHA(2)*X(2,IX)
            HTMPI   = ALPHA(1)*X(2,IX) + ALPHA(2)*X(1,IX)
            HTMPJ   = ALPHA(1)*X(3,IX) - ALPHA(2)*X(4,IX)
            X(4,IX) = ALPHA(1)*X(4,IX) + ALPHA(2)*X(3,IX)
            X(1,IX) = HTMPS
            X(2,IX) = HTMPI
            X(3,IX) = HTMPJ

            IX = IX + INCX
  30      CONTINUE
        ELSE
          DO 40 I = 1,N
C           Hamilton Product
            HTMPS   =  X(1,IX)*ALPHA(1) - X(2,IX)*ALPHA(2)
            HTMPI   =  X(1,IX)*ALPHA(2) + X(2,IX)*ALPHA(1)
            HTMPJ   =  X(3,IX)*ALPHA(1) + X(4,IX)*ALPHA(2)
            X(4,IX) = -X(3,IX)*ALPHA(2) + X(4,IX)*ALPHA(1)
            X(1,IX) = HTMPS
            X(2,IX) = HTMPI
            X(3,IX) = HTMPJ

            IX = IX + INCX
  40      CONTINUE
        ENDIF
      ENDIF
C
      END 
