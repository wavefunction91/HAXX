/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */
#ifndef __INCLUDED_SIMD_QOP_HPP__
#define __INCLUDED_SIMD_QOP_HPP__

#include <iostream>

namespace HAXX {

  // Single Quaternion Multiplication

  inline __m256d MULDQ_NN(__m256d &x, __m256d &y) {

    const __m256i maskScalar = _mm256_set_epi64x(0,0,0,0x8000000000000000);
    // Get the low and high-bits of x and y
    __m128d xLo = GET_LOD_256(x);
    __m128d yLo = GET_LOD_256(y);
    __m128d xHi = GET_HID_256(x);
    __m128d yHi = GET_HID_256(y);

    // 1123 Permutation of x
    __m256d x1123 = _mm256_permute_pd(x, 0xB); // 1011: 0123 -> 1123

    // 2231 Permutation of x
    __m128d rLo = _mm_shuffle_pd(xHi,xHi, 0x0); // 0000: 23;23 -> 22
    __m128d rHi = _mm_shuffle_pd(xHi,xLo, 0x3); // 0010: 23;01 -> 31
    __m256d x2231 = D256_FROM_D128(_mm_castpd_ps(rLo), _mm_castpd_ps(rHi));

    // 1000 Permutation of y
    rLo = _mm_shuffle_pd(yLo,yLo, 0x1); // 0001: 01;01 -> 11
    rHi = _mm_shuffle_pd(yLo,yLo, 0x0); // 0000: 01;01 -> 00
    __m256d y1000 = D256_FROM_D128(_mm_castpd_ps(rLo), _mm_castpd_ps(rHi));

    // 2312 Permutation of y
    rHi = _mm_shuffle_pd(yLo,yHi, 0x1); // 0001: 01;23 -> 12
    __m256d y2312 = D256_FROM_D128(_mm_castpd_ps(yHi), _mm_castpd_ps(rHi));


    __m256d t1  = MULD(x1123,y1000);
    __m256d t12 = FMAD(x2231,y2312,t1);

    // Flip scalar sign bit on t12
    t12 = _mm256_xor_pd(t12,_mm256_castsi256_pd(maskScalar));

    
    // 3312 Permutation of x
    rLo = _mm_shuffle_pd(xHi,xHi, 0x3); // 0010: 23;23 -> 33
    rHi = _mm_shuffle_pd(xLo,xHi, 0x1); // 0001: 01;23 -> 12
    __m256d x3312 = D256_FROM_D128(_mm_castpd_ps(rLo), _mm_castpd_ps(rHi));

    // 0000 Permutation of x
    rHi = _mm_shuffle_pd(xLo,xLo, 0x0); // 0000: 01;01 -> 00
    __m256d x0000 = D256_FROM_D128(_mm_castpd_ps(rHi), _mm_castpd_ps(rHi));


    // 3231 Permutation of y
    rLo = _mm_shuffle_pd(yHi,yHi, 0x1); // 0001: 23;23 -> 32
    rHi = _mm_shuffle_pd(yHi,yLo, 0x3); // 0010: 23;01 -> 31
    __m256d y3231 = D256_FROM_D128(_mm_castpd_ps(rLo), _mm_castpd_ps(rHi));

    

   __m256d t3  = MULD(x3312,y3231);
   __m256d t03 = FMSD(x0000,y,t3);

    return ADDD(t03,t12);

  };

  inline __m256d MULDQ_CN(__m256d &x, __m256d &y) {

    const __m256i maskVec    = _mm256_set_epi64x(0x8000000000000000,0x8000000000000000,0x8000000000000000,0);
    x = _mm256_xor_pd(x,_mm256_castsi256_pd(maskVec));
    return MULDQ_NN(x,y);

  }


  // Batched Quaternion Multiplication 
  
  /**
   *  \brief Increment each of the quaternion components
   *  of the quaternion product by their LHS scalar 
   *  contribution.
   *
   *  R(S) += X(S) * Y(S)
   *  R(I) += X(S) * Y(I)
   *  R(J) += X(S) * Y(J)
   *  R(K) += X(S) * Y(J)
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   */
  inline void INC_MULD4Q_SCALAR(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
  
    r1 = FMAD(x1,y1,r1);
    r2 = FMAD(x1,y2,r2);
    r3 = FMAD(x1,y3,r3);
    r4 = FMAD(x1,y4,r4);
  
  }
  
  /**
   *  \brief Initialize each of the quaternion components
   *  of the quaternion product by their LHS scalar 
   *  contribution.
   *
   *  R(S) = X(S) * Y(S)
   *  R(I) = X(S) * Y(I)
   *  R(J) = X(S) * Y(J)
   *  R(K) = X(S) * Y(J)
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   */
  inline void MULD4Q_SCALAR(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
  
    r1 = MULD(x1,y1);
    r2 = MULD(x1,y2);
    r3 = MULD(x1,y3);
    r4 = MULD(x1,y4);
  
  }
  
  /**
   *  \brief Increment the components of the unconjugated quaternion 
   *  product by their LHS vector contributions.
   *
   *
   *  R(S) -=  X(I)*Y(I) + X(J)*Y(J) + X(K)*Y(K)
   *  R(I) +=  X(I)*Y(S) + X(J)*Y(K) - X(K)*Y(J)
   *  R(J) += -X(I)*Y(K) + X(J)*Y(S) + X(K)*Y(I)
   *  R(K) +=  X(I)*Y(J) - X(J)*Y(I) - X(K)*Y(S)
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   *
   */
  inline void VEC_MULD4Q_NN(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
    r1 = FNMAD(x2,y2,r1); 
    r1 = FNMAD(x3,y3,r1); 
    r1 = FNMAD(x4,y4,r1); 
   
    r2 = FMAD(x2,y1,r2); 
    r2 = FMAD(x3,y4,r2); 
    r2 = FNMAD(x4,y3,r2); 
   
    r3 = FNMAD(x2,y4,r3); 
    r3 = FMAD(x3,y1,r3); 
    r3 = FMAD(x4,y2,r3); 
   
    r4 = FMAD(x2,y3,r4); 
    r4 = FNMAD(x3,y2,r4); 
    r4 = FMAD(x4,y1,r4);
  
  
  };
  
  
  /**
   *  \brief Increment the components of the left-conjugated quaternion 
   *  product by their LHS vector contributions.
   *
   *
   *  R(S) +=  X(I)*Y(I) + X(J)*Y(J) + X(K)*Y(K)
   *  R(I) += -X(I)*Y(S) - X(J)*Y(K) + X(K)*Y(J)
   *  R(J) +=  X(I)*Y(K) - X(J)*Y(S) - X(K)*Y(I)
   *  R(K) += -X(I)*Y(J) + X(J)*Y(I) + X(K)*Y(S)
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   *
   */
  inline void VEC_MULD4Q_CN(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
    r1 = FMAD(x2,y2,r1); 
    r1 = FMAD(x3,y3,r1); 
    r1 = FMAD(x4,y4,r1); 
   
    r2 = FNMAD(x2,y1,r2); 
    r2 = FNMAD(x3,y4,r2); 
    r2 = FMAD(x4,y3,r2); 
   
    r3 = FMAD(x2,y4,r3); 
    r3 = FNMAD(x3,y1,r3); 
    r3 = FNMAD(x4,y2,r3); 
   
    r4 = FNMAD(x2,y3,r4); 
    r4 = FMAD(x3,y2,r4); 
    r4 = FNMAD(x4,y1,r4);
  
  
  };
  
  
  /**
   *  Perform the packed unconjugated quaternion product of
   *  4 packed quaternions in place.
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   */
  inline void MULD4Q_NN(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
    MULD4Q_SCALAR(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
    VEC_MULD4Q_NN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
  
  };
  
  
  
  /**
   *  Perform the packed left-conjugated quaternion product of
   *  4 packed quaternions in place.
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   */
  inline void MULD4Q_CN(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
    MULD4Q_SCALAR(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
    VEC_MULD4Q_CN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
  
  };
  
  
  /**
   *  Perform the packed unconjugated quaternion product of
   *  4 packed quaternions while incrementing the result
   *  storage.
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   */
  inline void INC_MULD4Q_NN(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
    INC_MULD4Q_SCALAR(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
        VEC_MULD4Q_NN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
  
  };
  
  
  
  
  
  /**
   *  Perform the packed left-conjugated quaternion product of
   *  4 packed quaternions while incrementing the result
   *  storage.
   *
   *  \param[in]      x1  Packed X(S)
   *  \param[in]      x2  Packed X(I)
   *  \param[in]      x3  Packed X(J)
   *  \param[in]      x4  Packed X(K)
   *
   *  \param[in]      y1  Packed Y(S)
   *  \param[in]      y2  Packed Y(I)
   *  \param[in]      y3  Packed Y(J)
   *  \param[in]      y4  Packed Y(K)
   *
   *  \param[in/out]  x1  Packed R(S)
   *  \param[in/out]  x2  Packed R(I)
   *  \param[in/out]  x3  Packed R(J)
   *  \param[in/out]  x4  Packed R(K)
   */
  inline void INC_MULD4Q_CN(
    __m256d &x1, __m256d &x2, __m256d &x3, __m256d &x4,
    __m256d &y1, __m256d &y2, __m256d &y3, __m256d &y4,
    __m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4) {
  
    INC_MULD4Q_SCALAR(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
        VEC_MULD4Q_CN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
  
  };

};

#endif
