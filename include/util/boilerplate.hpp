/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */
#ifndef __INCLUDED_UTIL_BOILERPLATE_HPP__
#define __INCLUDED_UTIL_BOILERPLATE_HPP__

#include <tuple>

namespace HAXX {


  // Generate a sequence of pointers

  /**
   *  \brief Generate a tuple sequence of integers (implementation)
   *
   *  Implementation function for ptrSeq.
   */
  template <typename T, size_t... Is>
  constexpr auto pointer_sequence_impl(T seed, size_t INC, 
      std::index_sequence<Is...>) {
  
    return std::make_tuple( (seed + Is*INC)... );
  
  };
  
  /**
   *  \brief Generate a tuple sequence of integers
   *
   *  Generate a sequence of integers and store it in a tuple.
   */
  template <size_t N, typename T>
  constexpr auto pointer_sequence(T seed, size_t INC) {
  
    return pointer_sequence_impl(seed, INC, 
        std::make_index_sequence< N >{});
  
  };



  // Apply a function to a tuple (FIXME: deprecate in C++17)


  /**
   *  \brief Apply a function to a tuple and store the results in a
   *  tuple. (implementation)
   *
   *  Implementation function for apply
   */
  template < typename F, class Tuple, size_t... Is>
  constexpr auto apply_impl( const F &op, Tuple &t, 
      std::index_sequence<Is...>) {
  
    return std::make_tuple( op(std::get<Is>(t))... );
  
  };
  
  
  
  /**
   *  \brief Apply a function to a tuple and store the results in a
   *  tuple. 
   */
  template < typename F, class Tuple>
  constexpr auto apply( const F &op, Tuple &t) {
  
    return apply_impl(op,t,
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
  
  };



  /**
   *  \brief Apply a function N times to a tuple of size M (N >= M).
   *  (implementation case N > M).
   */ 
  template < size_t N, typename F, class Tuple, 
    typename std::enable_if_t< (N > 0), int> = 0 >
  constexpr auto apply_n_impl( const F &op, const Tuple &prev ) {
  
    return std::tuple_cat( apply_n_impl<N-1>(op,prev), std::make_tuple(op()) );
  
  };
  
  /**
   *  \brief Apply a function N times to a tuple of size M (N >= M).
   *  (implementation case N == M).
   *
   *  Delegates to apply
   */ 
  template < size_t N, typename F, class Tuple, 
    typename std::enable_if_t< (N == 0), int> = 0 >
  constexpr auto apply_n_impl( const F &op, const Tuple &prev ) {
  
    return apply(op,prev);
  
  };
  
  /**
   *  \brief Apply a function N times to a tuple of size M (N >= M).
   *
   *  For the remaining (N-M) executations of the function, call
   *  the function with no arguments (requires proper overload).
   */ 
  template < size_t N, typename F, class Tuple >
  constexpr auto apply_n( const F &op, const Tuple &t ) {
  
    static_assert( 
      (N >= std::tuple_size<Tuple>::value), 
      "N must be >= sizeof(Tuple)"
    );
  
    return apply_n_impl<N - std::tuple_size<Tuple>::value>(op,t); 
  
  };




  /**
   *  \brief Apply a binary function to a commonly indexed tuple
   *  and array. (implementation case sizeof(array) > 1).
   */ 
  template <typename F, typename T, typename U, typename... Args>
  constexpr inline std::enable_if_t<(sizeof...(Args) > 0),void> 
    array_tuple_execute_impl( const F &op, T *ptr, U &param, Args... args) {
  
    op(ptr,param); array_tuple_execute_impl(op,++ptr,args...);
  
  };

  /**
   *  \brief Apply a binary function to a commonly indexed tuple
   *  and array. (implementation case sizeof(array) == 1).
   */ 
  template <typename F, typename T, typename U, typename... Args>
  constexpr inline void array_tuple_execute_impl( const F &op, T *ptr, 
      U &param) {
  
    op(ptr,param); 
  
  };
  
  /**
   *  \brief Apply a binary function to a commonly indexed tuple
   *  and array. (expansion of tuple).
   */ 
  template <typename F, typename T, class Tuple, size_t... Is>
  constexpr inline void array_tuple_execute_impl( const F &op, T* ptr, 
    Tuple &t, std::index_sequence<Is...>) {
  
    array_tuple_execute_impl(op,ptr,(std::get<Is>(t))...);
  
  };
  
  /**
   *  \brief Apply a binary function to a commonly indexed tuple
   *  and array.
   *
   *  Compiles to
   *  op(ptr+0,std::get<0>(t));
   *  op(ptr+1,std::get<1>(t));
   *  op(ptr+2,std::get<2>(t));
   *  ...
   */ 
  template <typename F, typename T, class Tuple>
  constexpr inline void array_tuple_execute(const F &op, T* ptr, Tuple &t) {
  
    array_tuple_execute_impl(op,ptr,t,
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
  
  }



  template <class... Fs> struct overload_set;
  
  template <class F1, class... Fs>
  struct overload_set<F1,Fs...> : F1, overload_set<Fs...>::type {
  
    typedef overload_set type;
    overload_set(F1 head, Fs... tail) :
      F1(head), overload_set<Fs...>::type(tail...){}
  
    using F1::operator();
    using overload_set<Fs...>::type::operator();
  
  };
  
  template< class F >
  struct overload_set<F> : F {
  
    typedef F type;
    using F::operator();
  
  };
  
  
  template <class... Fs>
  typename overload_set<Fs...>::type overload(Fs... x) {
  
    return overload_set<Fs...>(x...);
  
  };


}; // namespace HAXX



#endif
