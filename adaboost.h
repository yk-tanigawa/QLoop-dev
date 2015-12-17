#ifndef __adaboost_H__
#define __adaboost_H__ 

#include <math.h>
#include "calloc_errchk.h"

inline void adaDataOnTheFly(const int h_i,
			 const int h_j,
			 const int **kmerFreq,
			 const int k,
			 short *x){
  const unsigned long nkmers = 1 << (2 * k);
  unsigned long l, m;
  for(l = 0; l < nkmers; l++){
    for(m = 0; m < nkmers; m++){
      x[l * nkmers + m] = (((kmerFreq[h_i][l] * kmerFreq[h_j][m]) > 0) ? 1 : 0);
    }
  }   
  return;
}

inline short adaDataOnTheFlyAxis(const int h_i,
			      const int h_j,
			      const int **kmerFreq,
			      const int k,
			      const int axis){
  const unsigned long nkmers = 1 << (2 * k);
  return (((kmerFreq[h_i][axis / nkmers] * kmerFreq[h_j][axis % nkmers]) > 0) ? 1 : 0);
}


int adaboostLearn(const int *y,
		  const int *h_i,
		  const int *h_j,
		  const int **kmerFreq,
		  const int k,
		  const unsigned long T,
		  const unsigned long N,
		  const unsigned long dim,
		  unsigned long **adaAxis,
		  int **adaSign,
		  double **adaBeta){
  short *selected, *x, xaxis;
  double *w, *p, *err, wsum, epsilon, min, max;
  unsigned long t, i, d, argmind, argmaxd;

  *adaAxis = calloc_errchk(T, sizeof(unsigned long), "calloc adaAxis");
  *adaSign = calloc_errchk(T, sizeof(int), "calloc adaSign");
  *adaBeta = calloc_errchk(T, sizeof(double), "calloc adaBeta");

  selected = calloc_errchk(dim, sizeof(short), "calloc selected");
  x = calloc_errchk(dim, sizeof(short), "calloc x");
  w = calloc_errchk(N, sizeof(double), "calloc w");
  p = calloc_errchk(N, sizeof(double), "calloc p");
  err = calloc_errchk(dim, sizeof(double), "calloc count");

  for(i = 0; i < N; i++){
    w[i] = 1.0 / N;
  }

  for(t = 0; t < T; t++){
    /* step 1 : compute normalized weights p[] */
    {
      wsum = 0;
      for(i = 0; i < N; i++){
	wsum += w[i];
      }
      for(i = 0; i < N; i++){
	p[i] = w[i] / wsum;
      }
    }

    /* step 2 : find the most appropriate axis (weak lerner) */
    {
      for(d = 0; d < dim; d++){
	err[d] = 0;
      }
      for(i = 0; i < N; i++){
	/* compute data point x */
	adaDataOnTheFly(h_i[i], h_j[i], kmerFreq, k, x);
	for(d = 0; d < dim; d++){	  
	  if(x[d]!= y[i]){
	    err[d] += p[i];
	  }
	}
      }

      {
	d = 0;
	while(selected[d] != 0){
	  d++;
	}
	max = min = err[d];
	argmaxd = argmind = d;
	for(d++; d < dim; d++){
	  if(selected[d] == 0){
	    if(err[d] < min){
	      min = err[d];
	      argmind = d;
	    }else if(err[d] > max){
	      max = err[d];
	      argmaxd = d;
	    }
	  }
	}
      }

      {
	if(max + min > 1.0){
	  /** 
	   * min > 1 - max 
	   *  argmaxd is the best axis
	   */
	  selected[argmaxd]++;
	  (*adaAxis)[t] = argmaxd;
	  (*adaSign)[t] = 1;
	  epsilon = 1 - max;
	}else{
	  /*  argmind is the best axis */
	  selected[argmind]++;
	  (*adaAxis)[t] = argmind;
	  (*adaSign)[t] = 0;
	  epsilon = min;       
	}      
      }
    }

    /* step 3: compute new weithgts */
    {
      (*adaBeta)[t] = epsilon / (1 - epsilon);
      for(i = 0; i < N; i++){
	/* compute data */
	xaxis = adaDataOnTheFlyAxis(h_i[i], h_j[i], kmerFreq, k, (*adaAxis)[t]);
	if(((*adaSign)[t] == 0 && xaxis == y[i]) ||
	   ((*adaSign)[t] == 1 && xaxis != y[i])){
	  w[i] *= (*adaBeta)[t];
	}
      }
    }
    
    fprintf(stderr, "t = %d\tbeta[t] = %e\tlerner = (%ld, %d)\n",
	    t, (*adaBeta)[t], (*adaAxis)[t], (*adaSign)[t]);	       
  }

  free(selected);
  free(x);
  free(w);
  free(p);
  free(err);
  return 0;
}


#if 0
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
      if(lernerPred[t] != x[i][lernerAxis[t]]){
	sum -= log(beta[t]);
      }
    }
    (*pred)[i] = ((sum >= threshold) ? 1 : 0);
  }
  return 0;
}
#endif

#endif
