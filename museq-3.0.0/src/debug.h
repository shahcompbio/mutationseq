#ifndef __debug_h
#define __debug_h

#include <stdio.h>
#include <errno.h>
#include <string.h>
#define trace(M, ...) fprintf(stderr, "\tTRACE %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#ifdef NDEBUG
#define debug(M, ...)
#else
#define trace(M, ...) fprintf(stderr, "\tTRACE %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#endif
#endif
