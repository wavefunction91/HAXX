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

  // Misc operations

  #define QCONJ_256D(X) \
    _mm256_xor_pd(X,\
      _mm256_castsi256_pd( \
        _mm256_set_epi64x(\
          0x8000000000000000,\
          0x8000000000000000,\
          0x8000000000000000,\
          0x0000000000000000\
        )\
      )\
    )


  // Single Quaternion Multiplication


  /**
   *  \brief Multiply two quaternions as 4 packed  64-bit floats
   *  int a 256-bit vector lane.
   *
   *  \param[in]  x  Packed LHS quaternion
   *  \param[in]  y  Packed RHS quaternion
   *
   *  \returns r = x * y
   */ 
  inline __m256d MULDQ_NN(const __m256d &x, const __m256d &y) {

    const __m256i maskScalar = _mm256_set_epi64x(0,0,0,0x8000000000000000);
    // Get the low and high-bits of x and y
    __m128d xLo = GET_LO_128D_256D(x);
    __m128d yLo = GET_LO_128D_256D(y);
    __m128d xHi = GET_HI_128D_256D(x);
    __m128d yHi = GET_HI_128D_256D(y);

    // 1123 Permutation of x
    __m256d x1123 = _mm256_permute_pd(x, 0xB); // 1011: 0123 -> 1123

    // 2231 Permutation of x
    __m128d rLo = _mm_shuffle_pd(xHi,xHi, 0x0); // 0000: 23;23 -> 22
    __m128d rHi = _mm_shuffle_pd(xHi,xLo, 0x3); // 0010: 23;01 -> 31
    __m256d x2231 = SET_256D_FROM_128D(rLo,rHi);

    // 1000 Permutation of y
    rLo = _mm_shuffle_pd(yLo,yLo, 0x1); // 0001: 01;01 -> 11
    rHi = _mm_shuffle_pd(yLo,yLo, 0x0); // 0000: 01;01 -> 00
    __m256d y1000 = SET_256D_FROM_128D(rLo,rHi);

    // 2312 Permutation of y
    rHi = _mm_shuffle_pd(yLo,yHi, 0x1); // 0001: 01;23 -> 12
    __m256d y2312 = SET_256D_FROM_128D(yHi,rHi);

    // t1 = x1123 * y1000
    __m256d t1  = _mm256_mul_pd(x1123,y1000);

    // t12 = t1 + x2231 * y2312 
    __m256d t12 = FMA_256D(x2231,y2312,t1);

    // Flip scalar sign bit on t12
    t12 = _mm256_xor_pd(t12,_mm256_castsi256_pd(maskScalar));

    
    // 3312 Permutation of x
    rLo = _mm_shuffle_pd(xHi,xHi, 0x3); // 0010: 23;23 -> 33
    rHi = _mm_shuffle_pd(xLo,xHi, 0x1); // 0001: 01;23 -> 12
    __m256d x3312 = SET_256D_FROM_128D(rLo,rHi);

    // 0000 Permutation of x
    rHi = _mm_shuffle_pd(xLo,xLo, 0x0); // 0000: 01;01 -> 00
    __m256d x0000 = SET_256D_FROM_128D(rHi,rHi);


    // 3231 Permutation of y
    rLo = _mm_shuffle_pd(yHi,yHi, 0x1); // 0001: 23;23 -> 32
    rHi = _mm_shuffle_pd(yHi,yLo, 0x3); // 0010: 23;01 -> 31
    __m256d y3231 = SET_256D_FROM_128D(rLo,rHi);

    
    // t3 = x3312 * y3231
    __m256d t3  = _mm256_mul_pd(x3312,y3231);

    // t03 = t3 + x000y * y
    __m256d t03 = FMS_256D(x0000,y,t3);

    // r = t03 + t12
    return _mm256_add_pd(t03,t12);

  };


  /**
   *  \brief Multiply two quaternions as 4 packed  64-bit floats
   *  int a 256-bit vector lane.
   *
   *  \param[in]  x  Packed LHS quaternion
   *  \param[in]  y  Packed RHS quaternion
   *
   *  \returns r = CONJ(x) * y
   */ 
  inline __m256d MULDQ_CN(__m256d &x, __m256d &y) {

/*
    const __m256i maskVec = 
      _mm256_set_epi64x(0x8000000000000000,0x8000000000000000,
                        0x8000000000000000,0);

    // Flip the sign bits on the vector part of x -> CONJ(x)
    x = _mm256_xor_pd(x,_mm256_castsi256_pd(maskVec));
*/

    // r = CONJ(x) * y
    return MULDQ_NN(QCONJ_256D(x),y);

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
  
  
    r1 = FMA_256D(x1,y1,r1);
    r2 = FMA_256D(x1,y2,r2);
    r3 = FMA_256D(x1,y3,r3);
    r4 = FMA_256D(x1,y4,r4);
  
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
  
  
    r1 = _mm256_mul_pd(x1,y1);
    r2 = _mm256_mul_pd(x1,y2);
    r3 = _mm256_mul_pd(x1,y3);
    r4 = _mm256_mul_pd(x1,y4);
  
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
  
    r1 = FNMA_256D(x2,y2,r1); 
    r1 = FNMA_256D(x3,y3,r1); 
    r1 = FNMA_256D(x4,y4,r1); 
   
    r2 = FMA_256D(x2,y1,r2); 
    r2 = FMA_256D(x3,y4,r2); 
    r2 = FNMA_256D(x4,y3,r2); 
   
    r3 = FNMA_256D(x2,y4,r3); 
    r3 = FMA_256D(x3,y1,r3); 
    r3 = FMA_256D(x4,y2,r3); 
   
    r4 = FMA_256D(x2,y3,r4); 
    r4 = FNMA_256D(x3,y2,r4); 
    r4 = FMA_256D(x4,y1,r4);
  
  
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
  
    r1 = FMA_256D(x2,y2,r1); 
    r1 = FMA_256D(x3,y3,r1); 
    r1 = FMA_256D(x4,y4,r1); 
   
    r2 = FNMA_256D(x2,y1,r2); 
    r2 = FNMA_256D(x3,y4,r2); 
    r2 = FMA_256D(x4,y3,r2); 
   
    r3 = FMA_256D(x2,y4,r3); 
    r3 = FNMA_256D(x3,y1,r3); 
    r3 = FNMA_256D(x4,y2,r3); 
   
    r4 = FNMA_256D(x2,y3,r4); 
    r4 = FMA_256D(x3,y2,r4); 
    r4 = FNMA_256D(x4,y1,r4);
  
  
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
