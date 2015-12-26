#ifndef __THRESHOLD_H__
#define __THRESHOLD_H__

#include <stdlib.h>
#include <stdio.h>
#include "constant.h"
#include "calloc_errchk.h"

int double_comp(const void *cmp1,
		const void *cmp2){
  double d1 = *((double *)cmp1);
  double d2 = *((double *)cmp2);
  if(d1 < d2){
    return -1;
  }else if(d1 > d2){
    return 1;
  }else{
    return 0;
  } 
}

double get_threshold(const double *mij,
		     const double percentile,
		     const unsigned long num){
  unsigned long i;
  double threshold = 0;
  int index = (int)(num * percentile), int_percentile = (int)(100 * percentile);
  double *cpy = calloc_errchk(num, sizeof(double), "calloc: cpy");
  for(i = 0; i < num; i++){
    cpy[i] = mij[i];
  }
  
  qsort(cpy, num, sizeof(double), double_comp);

  threshold = cpy[index];

  fprintf(stderr, "info: %d %% percentile = %f\n",
	  int_percentile, threshold);

  free(cpy);
  return threshold;
}

#endif
