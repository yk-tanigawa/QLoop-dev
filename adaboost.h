#ifndef __adaboost_H__
#define __adaboost_H__ 

#include <sys/time.h>
#include <math.h>
#include "calloc_errchk.h"
#include "diffSec.h"
#include "io.h"

typedef struct _adaboost{
  unsigned long *axis;
  double *beta;
  unsigned int *sign;
} adaboost;

/* arguments for function adaboost_comp_err */
typedef struct _adaboost_comp_err_args{
  /* thread specific info */
  int thread_id;
  unsigned long begin;
  unsigned long end;
  /* shared param(s) */
  unsigned long N;
  /* shared data */
  unsigned int **kmer_freq;
  unsigned int *h_i;
  unsigned int *h_j;
  /* array with 0.5 * 16^k elements */
  unsigned int *marked;
  unsigned int *l1;
  unsigned int *m1;
  unsigned int *l2;
  unsigned int *m2;
  double *err;
  /* array with N elements */
  double *p;
  unsigned int *y;
} adaboost_comp_err_args;

void *adaboost_comp_err(void *args){
  adaboost_comp_err_args *params = (adaboost_comp_err_args *)args;
  unsigned int kmerpair = 0, x = 0, pred;
 
#if 0
  fprintf(stderr, "thread %d: start [%d : %d]\n",
	  params->thread_id, params->begin, params->end);  
#endif

  for(kmerpair = params->begin; kmerpair <= params->end; kmerpair++){
    (params->err)[kmerpair] = 0;
  }
  for(kmerpair = params->begin; kmerpair <= params->end; kmerpair++){
    if(params->marked[kmerpair] == 0){
      for(x = 0; x < params->N; x++){
	pred = 
	  (params->kmer_freq)[(params->h_i)[x]][(params->l1)[kmerpair]] * 
	  (params->kmer_freq)[(params->h_j)[x]][(params->m1)[kmerpair]] +
	  (params->kmer_freq)[(params->h_i)[x]][(params->l2)[kmerpair]] * 
	  (params->kmer_freq)[(params->h_j)[x]][(params->m2)[kmerpair]];
	if((params->y)[x] != (pred > 0 ? 1 : 0)){
	  (params->err)[kmerpair] += (params->p)[x];
	}
      }
    }
  }  
  return NULL;
}

int set_kmer_pairs(unsigned int k,
		   unsigned int **l1,
		   unsigned int **m1,
		   unsigned int **l2,
		   unsigned int **m2){
  const unsigned int kmerpair_num = 1 << (4 * k - 1);  
  {
    *l1 = calloc_errchk(kmerpair_num, sizeof(unsigned int), "calloc l1");
    *m1 = calloc_errchk(kmerpair_num, sizeof(unsigned int), "calloc m1");
    *l2 = calloc_errchk(kmerpair_num, sizeof(unsigned int), "calloc l2");
    *m2 = calloc_errchk(kmerpair_num, sizeof(unsigned int), "calloc m2");
  }

  /**
   * need to fill four arrays
   */

  return 0;
}

int adaboost_learn(const command_line_arguements *cmd_args,
		   unsigned int **kmer_freq,
		   hic *hic,
		   adaboost **model){
  const unsigned long kmerpair_num = 1 << (4 * (cmd_args->k) - 1);  
  unsigned int *marked, *l1, *l2, *m1, *m2, *y;
  double *err, *w, *p, wsum, epsilon, min, max;

  {
    *model = calloc_errchk(1, sizeof(adaboost), "calloc adaboost");
    marked = calloc_errchk(kmerpair_num, sizeof(unsigned int), "calloc: marked");
    err = calloc_errchk(kmerpair_num, sizeof(double), "calloc: err");
    y = calloc_errchk(hic->nrow, sizeof(unsigned int), "calloc: y");
    p = calloc_errchk(hic->nrow, sizeof(double), "calloc: p");
    set_kmer_pairs(cmd_args->k, &l1, &m1, &l2, &m2);
  }

  if(cmd_args->exec_thread_num >= 1){
    int i = 0;
    //    unsigned int t;
    adaboost_comp_err_args *params;
    pthread_t *threads = NULL;

    {
      params = calloc_errchk(cmd_args->exec_thread_num,			   
			     sizeof(adaboost_comp_err_args),
			     "calloc: adaboost_comp_err_args");
      threads = calloc_errchk(cmd_args->exec_thread_num,			   
			      sizeof(pthread_t),
			      "calloc: threads");        
      /* set variables */
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	params[i].thread_id = i;
	params[i].begin = ((i == 0) ? 0 : params[i - 1].end + 1);
	params[i].end = ((i == (cmd_args->exec_thread_num - 1)) ?
			 kmerpair_num :
			 ((kmerpair_num / cmd_args->exec_thread_num) * (i + 1) - 1));
	params[i].N = hic->nrow;
	params[i].kmer_freq = kmer_freq;
	params[i].h_i = hic->i;
	params[i].h_j = hic->j;
	params[i].marked = marked;
	params[i].l1 = l1;
	params[i].m1 = m1;
	params[i].l2 = l2;
	params[i].m2 = m2;
	params[i].err = err;
	params[i].p = p;
	params[i].y = y;
      }
    }

    /* pthread create */
    for(i = 0; i < cmd_args->exec_thread_num; i++){
      pthread_create(&threads[i], NULL, 
		     adaboost_comp_err,
		     (void*)&params[i]);
    }
    
    /* pthread join */
    for(i = 0; i < cmd_args->exec_thread_num; i++){
      pthread_join(threads[i], NULL);
    }
  }

  return 0;
}



#if 0

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


int dump_results(FILE *fp,
		 const unsigned long *adaAxis,
		 const int *adaSign,
		 const double *adaBeta,
		 const unsigned long T,
		 const int k){
  const unsigned long nkmers = 1 << (2 * k);
  unsigned long t;
  char **kmerStrings;

  setKmerStrings(k, &kmerStrings);
  
  for(t = 0; t < T; t++){
    fprintf(fp, "%e\t%s\t%s\t%ld\t%d\n",
	    adaBeta[t],
	    kmerStrings[adaAxis[t] / nkmers],
	    kmerStrings[adaAxis[t] % nkmers],
	    adaAxis[t],
	    adaSign[t]);
  }

  for(t = 0; t < nkmers; t++){
    free(kmerStrings[t]);
  }
  free(kmerStrings);

  return 0;
}



int adaboostLearnOnTheFly(const unsigned int *y,
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


#endif
#endif
