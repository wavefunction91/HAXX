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



// Quaternion Multiplication macros
  

// Evaluate 4 quaternion products simultaneously and in place
#define MULD4Q(CONJ,x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
\
    r1 = MULD(x1,y1);\
    r2 = MULD(x1,y2);\
    r3 = MULD(x1,y3);\
    r4 = MULD(x1,y4);\
\
    VEC_MULD4Q_##CONJ(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);


// Evaluate 4 quaternion products simultaneously and increment
// the result storage
#define INC_MULD4Q(CONJ,x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
\
    r1 = FMAD(x1,y1,r1);\
    r2 = FMAD(x1,y2,r2);\
    r3 = FMAD(x1,y3,r3);\
    r4 = FMAD(x1,y4,r4);\
\
    VEC_MULD4Q_##CONJ(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);



// Evaluate the vector parts of 4 quaternion products
// simultaneously and increment the result storage
//
// R = X * Y
#define VEC_MULD4Q_NN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
\
    r1 = FSMD(r1,x2,y2);\
    r1 = FSMD(r1,x3,y3);\
    r1 = FSMD(r1,x4,y4);\
\
    r2 = FMAD(x2,y1,r2);\
    r2 = FMAD(x3,y4,r2);\
    r2 = FSMD(r2,x4,y3);\
\
    r3 = FSMD(r3,x2,y4);\
    r3 = FMAD(x3,y1,r3);\
    r3 = FMAD(x4,y2,r3);\
\
    r4 = FMAD(x2,y3,r4);\
    r4 = FSMD(r4,x3,y2);\
    r4 = FMAD(x4,y1,r4);

// Evaluate the vector parts of 4 quaternion products
// (where the LHS is conjugated) simultaneously and increment 
// the result storage
//
// R = CONJ(X) * Y
#define VEC_MULD4Q_CN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
\
    r1 = FMAD(x2,y2,r1);\
    r1 = FMAD(x3,y3,r1);\
    r1 = FMAD(x4,y4,r1);\
\
    r2 = FSMD(r2,x2,y1);\
    r2 = FSMD(r2,x3,y4);\
    r2 = FMAD(x4,y3,r2);\
\
    r3 = FMAD(x2,y4,r3);\
    r3 = FSMD(r3,x3,y1);\
    r3 = FSMD(r3,x4,y2);\
\
    r4 = FSMD(r4,x2,y3);\
    r4 = FMAD(x3,y2,r4);\
    r4 = FSMD(r4,x4,y1);



// Simultaneously multiply 4 quaternions (R = X * Y)
#define MULD4Q_NN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
  MULD4Q(NN,x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);

// Simultaneously multiply (and increment) 4 quaternions
// (R += X * Y)
#define INC_MULD4Q_NN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
  INC_MULD4Q(NN,x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);
  
// Simultaneously multiply 4 quaternions (R = CONJ(X) * Y)
#define MULD4Q_CN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
  MULD4Q(CN,x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);

// Simultaneously multiply (and increment) 4 quaternions
// (R += CONJ(X) * Y)
#define INC_MULD4Q_CN(x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4) \
  INC_MULD4Q(CN,x1,x2,x3,x4,y1,y2,y3,y4,r1,r2,r3,r4);

#endif
