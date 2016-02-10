#ifndef __dbg_h__
#define __dbg_h__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef NO_DEBUG
#define debug(ptr,M, ...)
#else
#define debug(ptr,M, ...) fprintf(ptr, "[DEBUG] (%s:%d): " M, __FILE__,__LINE__, ##__VA_ARGS__)
#endif

#define clean_errno() (errno == 0 ? "None" : strerror(errno))

#define log_err(ptr,M, ...) fprintf(ptr, "[ERROR] (%s:%d: errno: %s) " M, __FILE__,__LINE__, clean_errno(), ##__VA_ARGS__)

#define log_warn(ptr,M, ...) fprintf(ptr, "[WARN] (%s:%d: errno: %s) " M, __FILE__,__LINE__, clean_errno(), ##__VA_ARGS__)

#define log_info(ptr,M, ...) fprintf(ptr, "[INFO] (%s:%d) " M, __FILE__,__LINE__, ##__VA_ARGS__)

#define check_error(A,print,ptr,M, ...) if(A) {if(print) log_err(ptr,M,##__VA_ARGS__); errno=0; goto error;}

#define check_warn(print,ptr,M, ...) {if(print) log_warn(ptr,M,##__VA_ARGS__);}

#define check_info(print,ptr,M, ...) {if(print) log_info(ptr,M,##__VA_ARGS__);}

#define check_debug(print,ptr,M, ...) {if(print) debug(ptr,M,##__VA_ARGS__);}

#define sentinel(print,ptr,M, ...)  {if(print) log_err(ptr,M,##__VA_ARGS__); errno=0; goto error;}

#define check_mem(A,print,ptr) check_error(!(A),print,ptr,"Out of memory.")

#endif

