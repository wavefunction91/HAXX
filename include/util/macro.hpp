/*
 *  This file is a part of HAXX
 *  
 *  Copyright (c) 2017 David Williams-Young
 *  All rights reserved.
 *  
 *  See LICENSE.txt 
 */

#ifndef __INCLUDED_UTIL_MACRO_HPP
#define __INCLUDED_UTIL_MACRO_HPP

// Misc macros

// Compute rank-2 index
#define RANK2_INDX(i,j,N) ( (i) + (j)*(N) )

// Alignment checking
#define IS_ALIGNED(X,B) ( ((unsigned long)(X) & (B-1)) == 0 )

// Fix a number for mod arithmitic
#define FixMod(X,N) (( (X) % (N) ) ? (X) + (N) - ((X) % (N)) : (X))

#endif
