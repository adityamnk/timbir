/*************************************************************
*  
*  FILE:     typeutil.h
*  AUTHOR:   Ben Brame
*  DATE:     Jan. 24, 2009
*  PURPOSE:  Define architecture independent integer types
* 
*************************************************************/

#ifndef _TYPEUTIL_H_
#define _TYPEUTIL_H_

#ifdef _WIN64
#define __WINDOWS__
#endif

#ifdef _WIN32
#define __WINDOWS__
#endif

#ifdef __WINDOWS__

typedef char int8_t;
typedef short int16_t;
typedef int int32_t;
typedef long long int64_t;

typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;

#else

#include <inttypes.h>
#include <stdint.h>

#endif

#endif /* _TYPEUTIL_H_ */

