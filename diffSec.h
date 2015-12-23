#ifndef __diffSec_H__
#define __diffSec_H__

#include <sys/time.h>

/**
 * compute difference of struct timeval and return it in seconds
 */

double diffSec(const struct timeval t0, 
	       const struct timeval t1){
  time_t diff_sec = (t1.tv_sec - t0.tv_sec);
  suseconds_t diff_usec = (t1.tv_usec - t0.tv_usec);
  return (double)diff_sec + ((double)diff_usec * 1e-6 );
}

/**
 * copy functino for struct timeval
 */

void cpTimeval(const struct timeval source,
	       struct timeval *target){
  (*target).tv_sec = source.tv_sec;
  (*target).tv_usec = source.tv_usec;
  return;
}

#endif
