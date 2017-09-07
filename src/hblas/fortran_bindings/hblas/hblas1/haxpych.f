C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HAXPYCH(SIDE, N, ALPHA, X, INCX, Y, INCY)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      LOGICAL   MRIGHT
      CHARACTER SIDE
      INTEGER*4 N, INCX, INCY, I, IX, IY
      REAL*8    ALPHA(2), X(4,*), Y(4,*)
C
      IF ( N.LE.0 ) RETURN
      IF ( ALPHA(1).EQ.0.0D0 .AND. ALPHA(2).EQ.0.0D0 ) RETURN
C
      MRIGHT = LSAME('R',SIDE)
C
      IF ( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
        IF ( .NOT. MRIGHT ) THEN
          DO 10 I = 1,N
C           Hamilton Product
            Y(1,I) = Y(1,I) + ALPHA(1)*X(1,I) - ALPHA(2)*X(2,I)
            Y(2,I) = Y(2,I) + ALPHA(1)*X(2,I) + ALPHA(2)*X(1,I)
            Y(3,I) = Y(3,I) + ALPHA(1)*X(3,I) - ALPHA(2)*X(4,I)
            Y(4,I) = Y(4,I) + ALPHA(1)*X(4,I) + ALPHA(2)*X(3,I)
  10      CONTINUE
        ELSE
          DO 20 I = 1,N
C           Hamilton Product
            Y(1,I) = Y(1,I) + X(1,I)*ALPHA(1) - X(2,I)*ALPHA(2)
            Y(2,I) = Y(2,I) + X(1,I)*ALPHA(2) + X(2,I)*ALPHA(1)
            Y(3,I) = Y(3,I) + X(3,I)*ALPHA(1) + X(4,I)*ALPHA(2)
            Y(4,I) = Y(4,I) - X(3,I)*ALPHA(2) + X(4,I)*ALPHA(1)
  20      CONTINUE
        ENDIF
      ELSE

        IX = 1
        IY = 1
        IF ( INCX.LT.0 ) IX = (-N + 1)*INCX + 1
        IF ( INCY.LT.0 ) IY = (-N + 1)*INCY + 1

        IF ( .NOT. MRIGHT ) THEN
          DO 30 I = 1,N
C           Hamilton Product
            Y(1,IY) = Y(1,IY) + ALPHA(1)*X(1,IX) - ALPHA(2)*X(2,IX)
            Y(2,IY) = Y(2,IY) + ALPHA(1)*X(2,IX) + ALPHA(2)*X(1,IX)
            Y(3,IY) = Y(3,IY) + ALPHA(1)*X(3,IX) - ALPHA(2)*X(4,IX)
            Y(4,IY) = Y(4,IY) + ALPHA(1)*X(4,IX) + ALPHA(2)*X(3,IX)

            IX = IX + INCX
            IY = IY + INCY
  30      CONTINUE
        ELSE
          DO 40 I = 1,N
C           Hamilton Product
            Y(1,IY) = Y(1,IY) + X(1,IX)*ALPHA(1) - X(2,IX)*ALPHA(2)
            Y(2,IY) = Y(2,IY) + X(1,IX)*ALPHA(2) + X(2,IX)*ALPHA(1)
            Y(3,IY) = Y(3,IY) + X(3,IX)*ALPHA(1) + X(4,IX)*ALPHA(2)
            Y(4,IY) = Y(4,IY) - X(3,IX)*ALPHA(2) + X(4,IX)*ALPHA(1)

            IX = IX + INCX
            IY = IY + INCY
  40      CONTINUE
        ENDIF
      ENDIF
C
      END 
