#ifndef __adaboost_H__
#define __adaboost_H__ 

#include <math.h>
#include "calloc_errchk.h"

int adaboost_apply(const unsigned long *lernerAxis,
		   const int *lernerPred,
		   const double *beta,
		   const unsigned long T,
		   const unsigned long N,
		   const unsigned int **x,
		   int **pred){

  unsigned int i, t;
  double threshold = 0, sum;
  {
    for(t = 0; t < T; t++){
      threshold -= log(beta[t]);
    }
    threshold /= 2;
  }
  for(i = 0; i < N; i++){
    sum = 0;
    for(t = 0; t < T; t++){
      if(lernerPred[t] != get_bit(x[i], lernerAxis[t])){
	sum -= log(beta[t]);
      }
    }
    (*pred)[i] = ((sum >= threshold) ? 1 : 0);
  }
  return 0;
}

#endif
