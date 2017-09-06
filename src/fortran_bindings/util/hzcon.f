C 
C     This file is a part of HAXX
C     
C     Copyright (c) 2017 David Williams-Young
C     All rights reserved.
C     
C     See LICENSE.txt 
C
      Subroutine HZCON1(UPLO, M, N, A, LDA, B, LDB)
C
      CHARACTER UPLO
      LOGICAL NOTL
      INTEGER*4 M,N,LDA,LDB
      COMPLEX*16 A(2,LDA,*), B(LDB,*)
C
      LOGICAL LSAME
      EXTERNAL LSAME
C
C
      NOTL = LSAME(UPLO,'U')
C
      JB = 1
      Do 20 J = 1,N
        IB = 1
        Do 10 I = 1,M
          A(1,I,J) = B(IB,JB)     

          IF( NOTL ) THEN
            A(2,I,J) = B(IB,JB+1)     
          ELSE
            A(2,I,J) = -CONJG(B(IB+1,JB))     
          ENDIF

          IB = IB + 2
  10    CONTINUE
        JB = JB + 2
  20  CONTINUE
C
      END
C
C
C
      Subroutine HZCON2(UPLO, M, N, A, LDA, B, LDB)
C
      CHARACTER UPLO
      LOGICAL NOTL
      INTEGER*4 M,N,LDA,LDB
C
      Complex*16 A(2,LDA,*), B(LDB,*)
C
C
      LOGICAL LSAME
      EXTERNAL LSAME
C

      NOTL = LSAME(UPLO,'U')

      Do 10 J = 1,N
      Do 10 I = 1,M

        A(1,I,J) = B(I,J)     
        IF( NOTL ) THEN
          A(2,I,J) = B(I,J+N)     
        ELSE
          A(2,I,J) = -CONJG(B(I+M,J))     
        ENDIF


  10  CONTINUE
C
      END
