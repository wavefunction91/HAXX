C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HDEXP(M, N, A, LDA, B, LDB)
C
      IMPLICIT REAL*8(A-H,O,Z)
      DIMENSION A(4,LDA,*), B(LDB,*)
C
      JB = 1
      Do 20 J = 1,N
        IB = 1
        Do 10 I = 1,M
          B(IB,JB)     =  A(1,I,J)
          B(IB,JB+1)   =  A(2,I,J)
          B(IB,JB+2)   =  A(3,I,J)
          B(IB,JB+3)   =  A(4,I,J)

          B(IB+1,JB)   =  -A(2,I,J)
          B(IB+1,JB+1) =   A(1,I,J)
          B(IB+1,JB+2) =  -A(4,I,J)
          B(IB+1,JB+3) =   A(3,I,J)

          B(IB+2,JB)   =  -A(3,I,J)
          B(IB+2,JB+1) =   A(4,I,J)
          B(IB+2,JB+2) =   A(1,I,J)
          B(IB+2,JB+3) =  -A(2,I,J)

          B(IB+3,JB)   =  -A(4,I,J)
          B(IB+3,JB+1) =  -A(3,I,J)
          B(IB+3,JB+2) =   A(2,I,J)
          B(IB+3,JB+3) =   A(1,I,J)

          IB = IB + 4
  10    CONTINUE
        JB = JB + 4
  20  CONTINUE
C
      END
