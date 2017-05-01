C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HSCALC(SIDE, N, ALPHA, X, INCX)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      LOGICAL    MRIGHT
      CHARACTER  SIDE
      INTEGER*4  N, INCX, I, IX
      COMPLEX*16 ALPHA, X(2,*)
C
      COMPLEX*16 HTMPS, HTMPJ
C
      IF ( N.LE.0 ) RETURN
      IF ( ABS(ALPHA).EQ.0.0D0  ) RETURN
C
      MRIGHT = LSAME('R',SIDE)
C
      IF ( INCX.EQ.1  ) THEN
        IF ( .NOT. MRIGHT ) THEN
          DO 10 I = 1,N
C           Hamilton Product
            X(1,I) = ALPHA * X(1,I)
            X(2,I) = ALPHA * X(2,I)
  10      CONTINUE
        ELSE
          DO 20 I = 1,N
C           Hamilton Product
            X(1,I) = X(1,I) * ALPHA
            X(2,I) = X(2,I) * CONJG(ALPHA)
  20      CONTINUE
        ENDIF
      ELSE

        IX = 1
        IF ( INCX.LT.0 ) IX = (-N + 1)*INCX + 1

        IF ( .NOT. MRIGHT ) THEN
          DO 30 I = 1,N
C           Hamilton Product
            X(1,IX) = ALPHA * X(1,IX)
            X(2,IX) = ALPHA * X(2,IX)

            IX = IX + INCX
  30      CONTINUE
        ELSE
          DO 40 I = 1,N
C           Hamilton Product
            X(1,IX) = X(1,IX) * ALPHA
            X(2,IX) = X(2,IX) * CONJG(ALPHA)

            IX = IX + INCX
  40      CONTINUE
        ENDIF
      ENDIF
C
      END 
