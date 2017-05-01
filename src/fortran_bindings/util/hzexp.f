C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HZEXP(M, N, A, LDA, B, LDB)
C
      IMPLICIT COMPLEX*16(A-H,O,Z)
      DIMENSION A(2,LDA,*), B(LDB,*)
C
      JB = 1
      Do 20 J = 1,N
        IB = 1
        Do 10 I = 1,M
          B(IB,JB)     =  A(1,I,J)
          B(IB+1,JB+1) =  CONJG(A(1,I,J))
          B(IB,JB+1)   =  A(2,I,J)
          B(IB+1,JB)   = -CONJG(A(2,I,J))

          IB = IB + 2
  10    CONTINUE
        JB = JB + 2
  20  CONTINUE
C
      END
