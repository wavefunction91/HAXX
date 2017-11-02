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

#include <cassert>

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
  
  
    // Sanity Checks
      

    static_assert( PK > 0, "Packing dimension must be non-zero" );      

    static_assert(
      (sizeof(typename PackOps::load_t) % sizeof(T)) == 0,
      "The size of the loaded variable must be an integer multiple of T" 
    );

    constexpr size_t NLOAD = sizeof(typename PackOps::load_t) / sizeof(T);

    static_assert(
      ( PK % NLOAD ) == 0,
      "The packing dimension must be an integer multiple of NLOAD"
    );

  
  
    HAXX_INT i,j;
  
  
    T *_x,*_x1,*_x2;
    T *_xp;
    
  
    _x  = X;
    _xp = Xp;
  
    auto alpha = PackOps::TypeWrapper::cacheScalar(args...);
  
  
    auto load_preOP = 
      overload(
        [&](auto &t){ 
          auto x = PackOps::load(t);        
          return PackOps::preOP(x,alpha);
        },
        [&](){ 
          auto x = PackOps::load();        
          return PackOps::preOP(x,alpha);
        }
      );
  
    auto ptrInc = [](auto &t){ return t + 1; };
  
    auto store = [](auto *&p, auto &t){ PackOps::store(p,t);    };

    auto loopBody = [&](auto &xCol, auto *&xp) {

      for( i = 0; i < M; i++ ) {
  
        auto x   = apply_n<PK>( load_preOP, xCol );
        auto x_t = PackOps::OP(x);
  
        arrExc( store, xp, x_t );
  
        xCol = apply( ptrInc, xCol );
        xp += PK;
  
      }

    };
  
    j = N / PK;
    if( j > 0 )
    do {
  
      auto xCol = ptrSeq<PK>(_x,LDX);
      _x += PK*LDX;

      loopBody(xCol,_xp);
      
      j--;
  
    } while(j > 0);
  

    j = N % PK;
    if( j ) {
  
      #ifdef _NDEBUG
        #undef _NDEBUG
        assert(j == 1 /* This has yet to be worked out */);
        #define _NDEBUG
      #else
        assert(j == 1 /* This has yet to be worked out */);
      #endif
  
      auto xCol = std::make_tuple(_x);
      loopBody(xCol,_xp);

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
      overload(
        [&](auto &t){ 
          auto x = PackOps::load(t);        
          return PackOps::preOP(x,alpha);
        },
        [&](){ 
          auto x = PackOps::load();        
          return PackOps::preOP(x,alpha);
        }
      );
  
    auto ptrInc = [&](auto &t){ return t + LDX; };
  
    auto store = [](auto *&p, auto &t){ PackOps::store(p,t);    };


    auto loopBody = [&](auto &xCol, auto *&xp) {

      for( j = 0; j < N; j++ ) {
  
        auto x   = apply_n<PK>( load_preOP, xCol );
        auto x_t = PackOps::OP(x);
        arrExc( store, xp, x_t );
  
        xCol = apply( ptrInc, xCol );
        xp += PK;
  
      }

    };

  
    i = M / PK;
    if( i > 0 )
    do {
  
      auto xCol = ptrSeq<PK>(_x,1);
      _x += PK;

      loopBody(xCol,_xp); 
  
      i--;
    } while( i > 0 );
  

    i = M % PK;
    if( i ) {
  
      #ifdef _NDEBUG
        #undef _NDEBUG
        assert(i==1 /* This has yet to be worked out */);
        #define _NDEBUG
      #else
        assert(i==1 /* This has yet to be worked out */);
      #endif

      auto xCol = std::make_tuple(_x);
      loopBody(xCol,_xp);
  
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
