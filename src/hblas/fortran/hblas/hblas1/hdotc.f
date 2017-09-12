C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HDOTC(R, N, X, INCX, Y, INCY)
C
      INTEGER*4 N, INCX, INCY, I,IX,IY
      REAL*8    X(4,*), Y(4,*), R(4)
C
      R(1) = 0.0D+0
      R(2) = 0.0D+0
      R(3) = 0.0D+0
      R(4) = 0.0D+0
C
      IF ( N.LE.0 ) RETURN
      IF ( INCX.EQ.1 .AND. INCY.EQ.1 ) THEN
        DO 10 I = 1,N
C         Hamilton Product
          R(1) = R(1) + X(1,I)*Y(1,I) + X(2,I)*Y(2,I) +
     $                  X(3,I)*Y(3,I) + X(4,I)*Y(4,I)

          R(2) = R(2) + X(1,I)*Y(2,I) - X(2,I)*Y(1,I) -
     $                  X(3,I)*Y(4,I) + X(4,I)*Y(3,I)

          R(3) = R(3) + X(1,I)*Y(3,I) + X(2,I)*Y(4,I) -
     $                  X(3,I)*Y(1,I) - X(4,I)*Y(2,I)

          R(4) = R(4) + X(1,I)*Y(4,I) - X(2,I)*Y(3,I) +
     $                  X(3,I)*Y(2,I) - X(4,I)*Y(1,I)
  10    CONTINUE
      ELSE

        IX = 1
        IY = 1
        IF ( INCX.LT.0 ) IX = (-N + 1)*INCX + 1
        IF ( INCY.LT.0 ) IY = (-N + 1)*INCY + 1

        DO 20 I = 1,N
C         Hamilton Product
          R(1) = R(1) + X(1,IX)*Y(1,IY) + X(2,IX)*Y(2,IY) +
     $                  X(3,IX)*Y(3,IY) + X(4,IX)*Y(4,IY)

          R(2) = R(2) + X(1,IX)*Y(2,IY) - X(2,IX)*Y(1,IY) -
     $                  X(3,IX)*Y(4,IY) + X(4,IX)*Y(3,IY)

          R(3) = R(3) + X(1,IX)*Y(3,IY) + X(2,IX)*Y(4,IY) -
     $                  X(3,IX)*Y(1,IY) - X(4,IX)*Y(2,IY)

          R(4) = R(4) + X(1,IX)*Y(4,IY) - X(2,IX)*Y(3,IY) +
     $                  X(3,IX)*Y(2,IY) - X(4,IX)*Y(1,IY)
          IX = IX + INCX
          IY = IY + INCY
  20    CONTINUE
      ENDIF
C
      END 
