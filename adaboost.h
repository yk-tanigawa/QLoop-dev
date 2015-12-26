#ifndef __adaboost_H__
#define __adaboost_H__ 

#include <sys/time.h>
#include <math.h>
#include "calloc_errchk.h"
#include "diffSec.h"
#include "io.h"
#include "kmer.h"


/* adaboost results*/
typedef struct _adaboost{
  unsigned long T;
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
  const unsigned int **kmer_freq;
  const unsigned int *h_i;
  const unsigned int *h_j;
  /* array with 2^(4k-1) + 2^(2k-1) elements */
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

int adaboost_show_itr(FILE *fp, 
		      const adaboost *model,
		      const char **kmer_strings,
		      const unsigned int *l1,
		      const unsigned int *m1,
		      const unsigned int *l2,
		      const unsigned int *m2,
		      const unsigned long t){
  fprintf(fp, "%ld\t%e\t%d\t%ld\t%s\t%s\t%s\t%s\n",
	  t, 
	  (model->beta)[t],
	  (model->sign)[t],
	  (model->axis)[t],
	  kmer_strings[l1[(model->axis)[t]]],
	  kmer_strings[m1[(model->axis)[t]]],
	  kmer_strings[l2[(model->axis)[t]]],
	  kmer_strings[m2[(model->axis)[t]]]);
  return 0;
}
		      
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

int set_kmer_pairs(const unsigned int k,
		   unsigned int **l1,
		   unsigned int **m1,
		   unsigned int **l2,
		   unsigned int **m2){
  const unsigned long kmer_num = 1 << (2 * k);
  const unsigned long kmer_pair_num = 1 << (4 * k);
  const unsigned long canonical_kmer_pair_num = (1 << (4 * k - 1)) + (1 << (2 * k - 1));  

  {
    *l1 = calloc_errchk(canonical_kmer_pair_num, sizeof(unsigned int), "calloc l1");
    *m1 = calloc_errchk(canonical_kmer_pair_num, sizeof(unsigned int), "calloc m1");
    *l2 = calloc_errchk(canonical_kmer_pair_num, sizeof(unsigned int), "calloc l2");
    *m2 = calloc_errchk(canonical_kmer_pair_num, sizeof(unsigned int), "calloc m2");
  }

  {
    unsigned long lm = 0, revcomp_lm = 0, next = 0;
    for(lm = 0; lm < kmer_pair_num; lm++){
      revcomp_lm = rev_comp(lm, 2 * k);
      if(lm <= revcomp_lm){
	/**
	 * where, l2 = rev_comp(m1, k)
	 *        m2 = rev_comp(l1, k)
	 * note:
	 *        concat(l2 + m2) = rev_comp(concat(l1 + m1))
	 */	
	(*l1)[next] = lm / kmer_num;
	(*m1)[next] = lm % kmer_num;
	(*l2)[next] = revcomp_lm / kmer_num;
	(*m2)[next] = revcomp_lm % kmer_num;
	next++;
      }
    }
  }
  return 0;
}

int adaboost_learn(const command_line_arguements *cmd_args,
		   const unsigned int **kmer_freq,
		   hic *hic,
		   adaboost **model){
  const unsigned long canonical_kmer_pair_num = 
    (1 << (4 * (cmd_args->k) - 1)) + (1 << (2 * (cmd_args->k) - 1));  
  unsigned long n, lm, argmin_lm, argmax_lm;
  unsigned int *marked, *l1, *l2, *m1, *m2, *y, pred;
  double *err, *w, *p, wsum, epsilon, min, max;
  char **kmer_strings;

  /* allocate memory */
  {
    *model = calloc_errchk(1, sizeof(adaboost), "calloc adaboost");
    marked = calloc_errchk(canonical_kmer_pair_num, sizeof(unsigned int), "calloc: marked");
    set_kmer_pairs(cmd_args->k, &l1, &m1, &l2, &m2);
    y = calloc_errchk(hic->nrow, sizeof(unsigned int), "calloc: y");
    err = calloc_errchk(canonical_kmer_pair_num, sizeof(double), "calloc: err");
    w = calloc_errchk(hic->nrow, sizeof(double), "calloc: p");
    p = calloc_errchk(hic->nrow, sizeof(double), "calloc: p");
    for(n = 0; n < hic->nrow; n++){
      w[n] = 1.0 / (hic->nrow);
    }
    set_kmer_strings(cmd_args->k, &kmer_strings);
  }

  if(cmd_args->exec_thread_num >= 1){
    int i = 0;
    unsigned long t;
    adaboost_comp_err_args *params;
    pthread_t *threads = NULL;

    /* prepare for thread programming */
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
	params[i].end =
	  ((i == (cmd_args->exec_thread_num - 1)) ?
	   canonical_kmer_pair_num :
	   ((canonical_kmer_pair_num / cmd_args->exec_thread_num) * (i + 1) - 1));
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

    /* AdaBoost iterations */
    for(t = 0; t < cmd_args->iteration_num; t++){
      /* step 1 : compute normalized weights p[] */
      {
	wsum = 0;
	for(n = 0; n < (hic->nrow); n++){
	  wsum += w[n];
	}
	for(n = 0; n < (hic->nrow); n++){
	  p[n] = w[n] / wsum;
	}
      }

      /* step 2 : find the most appropriate axis (weak lerner) */
      {
	/* compute err for each kmer pair using pthread */
	{
	  /* pthread create */
	  for(i = 0; i < cmd_args->exec_thread_num; i++){
	    pthread_create(&threads[i], NULL, adaboost_comp_err, (void*)&params[i]);	
	  }      
	  /* pthread join */
	  for(i = 0; i < cmd_args->exec_thread_num; i++){
	    pthread_join(threads[i], NULL);
	  }
	}
	/* find best stamp */
	{
	  lm = 0;
	  /* skip arleady selected kmer pairs */
	  while(marked[lm] != 0){
	    lm++;
	  }
	  /* find max and min*/
	  max = min = err[lm];
	  argmax_lm = argmin_lm = lm;
	  for(lm++; lm < canonical_kmer_pair_num; lm++){
	    if(marked[lm] == 0){
	      if(err[lm] < min){
		min = err[lm];
		argmin_lm = lm;
	      }else if(err[lm] > max){
		max = err[lm];
		argmax_lm = lm;
	      }
	    }
	  }
	  /* compare max and min */
	  {
	    if(max + min > 1.0){
	      /** 
	       * min > 1 - max 
	       *  argmaxd is the best axis
	       */
	      marked[argmax_lm]++;
	      ((*model)->axis)[t] = argmax_lm;
	      ((*model)->sign)[t] = 1;
	      epsilon = 1 - max;
	    }else{
	      /*  argmind is the best axis */
	      marked[argmin_lm]++;
	      ((*model)->axis)[t] = argmin_lm;
	      ((*model)->sign)[t] = 0;
	      epsilon = min;       
	    }      	    
	  }
	}
      }
      /* step 3 : compute new weights */
      {
	((*model)->beta)[t] = epsilon / (1 - epsilon);
	for(n = 0; n < hic->nrow; n++){
	  pred = 
	    ((kmer_freq[hic->i[n]][l1[((*model)->axis)[t]]] * 
	      kmer_freq[hic->j[n]][m1[((*model)->axis)[t]]] +
	      kmer_freq[hic->i[n]][l2[((*model)->axis)[t]]] * 
	      kmer_freq[hic->j[n]][m2[((*model)->axis)[t]]]) > 0) ? 1 : 0;
	  if(((((*model)->sign)[t] == 0) && pred == y[n]) ||
	     ((((*model)->sign)[t] == 1) && pred != y[n]))
	    w[n] *= ((*model)->beta)[t];
	}
      }
      adaboost_show_itr(stderr, *model, (const char**)kmer_strings, l1, m1, l2, m2, t);
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
