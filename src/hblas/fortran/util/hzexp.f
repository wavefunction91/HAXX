C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HZEXP1(M, N, A, LDA, B, LDB)
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
C
C
C
      Subroutine HZEXP2(M, N, A, LDA, B, LDB)
C
      IMPLICIT COMPLEX*16(A-H,O,Z)
      DIMENSION A(2,LDA,*), B(LDB,*)
C
      Do 10 J = 1,N
      Do 10 I = 1,M
        B(I,J)     =  A(1,I,J)
        B(I+M,J+N) =  CONJG(A(1,I,J))
        B(I,J+N)   =  A(2,I,J)
        B(I+M,J)   = -CONJG(A(2,I,J))
  10  CONTINUE
C
      END
