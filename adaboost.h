#ifndef __adaboost_H__
#define __adaboost_H__ 

#include <sys/time.h>
#include <math.h>
#include "calloc_errchk.h"
#include "diffSec.h"
#include "bit_op.h"

inline void adaDataOnTheFly(const int h_i,
			    const int h_j,
			    const int **kmerFreq,
			    const int k,
			    unsigned short *x){
  const unsigned long nkmers = 1 << (2 * k);
  unsigned long l, m;
  for(l = 0; l < nkmers; l++){
    for(m = 0; m < nkmers; m++){
      x[l * nkmers + m] = (((kmerFreq[h_i][l] * kmerFreq[h_j][m]) > 0) ? 1 : 0);
    }
  }   
  return;
}

inline unsigned short adaDataOnTheFlyAxis(const int h_i,
					  const int h_j,
					  const int **kmerFreq,
					  const int k,
					  const int axis){
  const unsigned long nkmers = 1 << (2 * k);
  return (((kmerFreq[h_i][axis / nkmers] * kmerFreq[h_j][axis % nkmers]) > 0) ? 1 : 0);
}


int adaboostLearn(const unsigned int *y,
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
  short *selected;
  unsigned short *x, xaxis;
  double *w, *p, *err, wsum, epsilon, min, max;
  unsigned long t, i, d, argmind, argmaxd;
  struct timeval timeStart, timePrev, time;

  *adaAxis = calloc_errchk(T, sizeof(unsigned long), "calloc adaAxis");
  *adaSign = calloc_errchk(T, sizeof(int), "calloc adaSign");
  *adaBeta = calloc_errchk(T, sizeof(double), "calloc adaBeta");

  selected = calloc_errchk(dim, sizeof(short), "calloc selected");
  x = calloc_errchk(dim, sizeof(unsigned short), "calloc x");
  w = calloc_errchk(N, sizeof(double), "calloc w");
  p = calloc_errchk(N, sizeof(double), "calloc p");
  err = calloc_errchk(dim, sizeof(double), "calloc count");

  for(i = 0; i < N; i++){
    w[i] = 1.0 / N;
  }

  gettimeofday(&timeStart, NULL);
  cpTimeval(timeStart, &timePrev);
  
  fprintf(stderr, "t\tbeta\t\tweak lerner\ttime(itration)\ttime(total) [s]\n");

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
    
  gettimeofday(&time, NULL);
    
    fprintf(stderr, "%ld\t%e\t(%ld, %d)\t%e\t%e\n",
	    t, (*adaBeta)[t], (*adaAxis)[t], (*adaSign)[t],
	    diffSec(timePrev, time), diffSec(timeStart, time));	       

    cpTimeval(time, &timePrev);
  }

  free(selected);
  free(x);
  free(w);
  free(p);
  free(err);
  return 0;
}

#if 0
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
#endif

/* ada boost : bit mode */

int constructBitTable(const int *h_i, 
		      const int *h_j, 
		      const int **kmerFreq,
		      const unsigned long nHic, 
		      const int k,
		      unsigned int ***x){
  const unsigned long nkmers = 1 << (2 * k);
  const unsigned long nkmerpairs = 1 << (4 * k);
  const size_t intBits = 8 * sizeof(unsigned int);
  unsigned long n, l, m, lm;
  unsigned int buf = 0;
  {
    *x = calloc_errchk(nHic, sizeof(unsigned int *), "calloc x");
    for(n = 0; n < nHic; n++){
      (*x)[n] = calloc_errchk(((nkmerpairs + intBits - 1) / intBits), sizeof(unsigned int), "calloc x[n]");
    }
  }

  for(n = 0; n < nHic; n++){
    lm = 0;
    buf = 0;
    for(l = 0; l < nkmers; l++){
      for(m = 0; m < nkmers; m++){
	buf = buf << 1;
	buf += (((kmerFreq[h_i[n]][l] * kmerFreq[h_j[n]][m]) > 0) ? 1 : 0);
	if((lm % intBits) == (intBits - 1)){
	  (*x)[n][lm / intBits] = buf;
	  buf = 0;
	}
	lm++;
      }
    }   
  }
  

  return 0;
}

int adaboostBitLearn(const unsigned int *y,
		     const unsigned int **x,
		     const unsigned long T,
		     const unsigned long N,
		     const unsigned long dim,
		     unsigned long **adaAxis,
		     int **adaSign,
		     double **adaBeta){
  short *selected;
  double *w, *p, *err, wsum, epsilon, min, max;
  unsigned long t, i, d, argmind, argmaxd;
  struct timeval timeStart, timePrev, time;

  *adaAxis = calloc_errchk(T, sizeof(unsigned long), "calloc adaAxis");
  *adaSign = calloc_errchk(T, sizeof(int), "calloc adaSign");
  *adaBeta = calloc_errchk(T, sizeof(double), "calloc adaBeta");

  selected = calloc_errchk(dim, sizeof(short), "calloc selected");
  w = calloc_errchk(N, sizeof(double), "calloc w");
  p = calloc_errchk(N, sizeof(double), "calloc p");
  err = calloc_errchk(dim, sizeof(double), "calloc count");

  for(i = 0; i < N; i++){
    w[i] = 1.0 / N;
  }

  gettimeofday(&timeStart, NULL);
  cpTimeval(timeStart, &timePrev);
  
  fprintf(stderr, "t\tbeta\t\tweak lerner\ttime(itration)\ttime(total) [s]\n");

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
	for(d = 0; d < dim; d++){	  
	  if(get_bit(x[i], d) != y[i]){
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
	if(((*adaSign)[t] == 0 && get_bit(x[i], (*adaAxis)[t]) == y[i]) ||
	   ((*adaSign)[t] == 1 && get_bit(x[i], (*adaAxis)[t]) != y[i])){
	  w[i] *= (*adaBeta)[t];
	}
      }
    }
    
  gettimeofday(&time, NULL);
    
    fprintf(stderr, "%ld\t%e\t(%ld, %d)\t%e\t%e\n",
	    t, (*adaBeta)[t], (*adaAxis)[t], (*adaSign)[t],
	    diffSec(timePrev, time), diffSec(timeStart, time));	       

    cpTimeval(time, &timePrev);
  }

  free(selected);
  free(w);
  free(p);
  free(err);
  return 0;
}


#if 0

/** 
 *   The following code is copied from another project, 
 * and hence does not cabaple with this program 
 */
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
