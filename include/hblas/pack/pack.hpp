/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */
#ifndef __INCLUDED_HBLAS_PACK_PACK_HPP
#define __INCLUDED_HBLAS_PACK_PACK_HPP

#include "haxx.hpp"
#include "util/simd.hpp"
#include "util/boilerplate.hpp"

#include "hblas/pack/typewrapper.hpp"
#include "hblas/pack/packops.hpp"


namespace HAXX {

  /**
   *  \brief Packing function for packing along the fastest running index
   *  of a matrix.
   */ 
  template <size_t PK, typename T, typename PackOps, typename... Args>
  void TPACK(const HAXX_INT M, const HAXX_INT N, T *X,
    const HAXX_INT LDX, T *Xp, Args... args) {
  
  
    // Sanity Check
    static_assert(
      (sizeof(typename PackOps::load_t) % sizeof(T)) == 0,
      "The size of the loaded variable must be an integer multiple of T" 
    );
  
  
  
    HAXX_INT i,j;
  
  
    T *_x,*_x1,*_x2;
    T *_xp;
    
  
    _x  = X;
    _xp = Xp;
  
    auto alpha = PackOps::TypeWrapper::cacheScalar(args...);
  
  
    auto load_preOP = 
      [&](auto &t){ 
        auto x = PackOps::load(t);        
        return PackOps::preOP(x,alpha);
      };
  
    auto ptrInc = [](auto &t){ return t + 1; };
  
    auto store = [](auto *&p, auto &t){ PackOps::store(p,t);    };
  
    j = N / PK;
    if( j > 0 )
    do {
  
      auto xCol = ptrSeq<PK>(_x,LDX);
  
      _x += PK*LDX;
  
      for( i = 0; i < M; i++ ) {
  
        auto x   = apply( load_preOP, xCol );
        auto x_t = PackOps::OP(x);
  
        arrExc( store, _xp, x_t );
  
        xCol = apply( ptrInc, xCol );
        _xp += PK;
  
      }
      
      j--;
  
    } while(j > 0);
  
    if( N % PK ) {
  
      #ifdef _NDEBUG
        #undef _NDEBUG
        assert(false /* This has yet to be worked out */);
        #define _NDEBUG
      #endif
  
  #if 0
      for( i = 0; i < M; i++ ) {
  
        auto x1 = PackOps::load(_x);
        auto x2 = PackOps::load();
  
        x1 = PackOps::preOP(x1,alpha);
        x2 = PackOps::preOP(x2,alpha);
  
        auto x1_t = PackOps::OP1(x1, x2);
        auto x2_t = PackOps::OP2(x1, x2);
  
        PackOps::store(_xp,   x1_t);
        PackOps::store(_xp+1, x2_t);
  
        _xp +=2; _x++;
  
      }
  #endif
  
    }
  
  }
  
  /**
   *  \brief Packing function for packing along the slowest running index
   *  of a matrix.
   */ 
  template <size_t PK, typename T, typename PackOps, typename... Args>
  inline void NPACK(const HAXX_INT M, const HAXX_INT N, T *X, 
    const HAXX_INT LDX, T *Xp, Args... args) {
  
    // Sanity Check
    static_assert(
      (sizeof(typename PackOps::load_t) % sizeof(T)) == 0,
      "The size of the loaded variable must be an integer multiple of T" 
    );
  
    HAXX_INT i,j;
    T *_x,*_x1,*_x2;
    T *_xp;
  
    _x = X;
    _xp = Xp;
  
    auto alpha = PackOps::TypeWrapper::cacheScalar(args...);
  
    auto load_preOP = 
      [&](auto &t){ 
        auto x = PackOps::load(t);        
        return PackOps::preOP(x,alpha);
      };
  
    auto ptrInc = [&](auto &t){ return t + LDX; };
  
    auto store = [](auto *&p, auto &t){ PackOps::store(p,t);    };
  
    i = M / PK;
    if( i > 0 )
    do {
  
      auto xCol = ptrSeq<PK>(_x,1);
  
      _x += PK;
  
      for( j = 0; j < N; j++ ) {
  
        auto x   = apply( load_preOP, xCol );
        auto x_t = PackOps::OP(x);
        arrExc( store, _xp, x_t );
  
        xCol = apply( ptrInc, xCol );
        _xp += PK;
  
      }
  
      i--;
    } while( i > 0 );
  
    if( M % PK ) {
  
      #ifdef _NDEBUG
        #undef _NDEBUG
        assert(false /* This has yet to be worked out */);
        #define _NDEBUG
      #endif
  
  #if 0
      for( j = 0; j < N; j++ ) {
  
        auto x1 = PackOps::load(_x);
        auto x2 = PackOps::load();
  
        x1 = PackOps::preOP(x1,alpha);
        x2 = PackOps::preOP(x2,alpha);
  
        auto x1_t = PackOps::OP1(x1, x2);
        auto x2_t = PackOps::OP2(x1, x2);
  
        PackOps::store(_xp,   x1_t);
        PackOps::store(_xp+1, x2_t);
  
        _x += LDX; _xp += 2;
  
      }
  #endif
  
    }
  
  };



}; // namespace HAXX

#endif
