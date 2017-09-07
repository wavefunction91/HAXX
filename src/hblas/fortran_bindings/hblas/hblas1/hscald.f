C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HSCALD(SIDE, N, ALPHA, X, INCX)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
      CHARACTER SIDE
      INTEGER*4 N, INCX, I, IX
      REAL*8    ALPHA, X(4,*)
C
      IF ( N.LE.0 ) RETURN
      IF ( ALPHA.EQ.0.0D0 ) RETURN
C
      IF ( INCX.EQ.1 ) THEN
        DO 10 I = 1,N
          X(1,I) = ALPHA*X(1,I)
          X(2,I) = ALPHA*X(2,I)
          X(3,I) = ALPHA*X(3,I)
          X(4,I) = ALPHA*X(4,I)
  10    CONTINUE
      ELSE

        IX = 1
        IF ( INCX.LT.0 ) IX = (-N + 1)*INCX + 1

        DO 20 I = 1,N
          X(1,IX) = ALPHA*X(1,IX)
          X(2,IX) = ALPHA*X(2,IX)
          X(3,IX) = ALPHA*X(3,IX)
          X(4,IX) = ALPHA*X(4,IX)

          IX = IX + INCX
  20    CONTINUE
      ENDIF
C
      END
