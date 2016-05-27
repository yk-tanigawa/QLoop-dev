#ifndef __l2boost_H__
#define __l2boost_H__ 

#include <sys/time.h>
#include <math.h>
#include <pthread.h>
#include "calloc_errchk.h"
#include "diffSec.h"
#include "kmer.h"
#include "hic.h"


/* L2 boost results*/
typedef struct _l2boost{
  double *res_sq;
  double *beta;
} l2boost;

typedef struct _cmpUdX_args{
  /* thread specific info */
  int thread_id;
  unsigned long begin;
  unsigned long end;
  /* shared param(s) */
  unsigned long n;
  /* shared data */
  const double **feature;
  const hic *data;
  const canonical_kp *ckps;
  /* shared data */
  double *U;
  double *UdX;
  double *Xnorm;
} cmpUdX_args;

void *cmpUdX(void *args){
  /* unstack parameters */
  const cmpUdX_args *params = (cmpUdX_args *)args;
  const double **feature = params->feature;
  const unsigned int *h_i = params->data->i;
  const unsigned int *h_j = params->data->j;
  const unsigned int *kmer1 = params->ckps->kmer1;
  const unsigned int *kmer2 = params->ckps->kmer2;
  const unsigned int *revcmp1 = params->ckps->revcmp1;
  const unsigned int *revcmp2 = params->ckps->revcmp2;

  /* compute the dot product between U and X^{(j)} */
  unsigned int i, j;
  double sum;
  for(j = params->begin; j < params->end; j++){
    sum = 0;
    for(i = 0; i < params->n; i++){
      sum += (params->U)[i] * ((feature[h_i[i]][kmer1[j]] *
				feature[h_j[i]][kmer2[j]]) +
			       (feature[h_i[i]][revcmp1[j]] *
				feature[h_j[i]][revcmp2[j]]));
    }
    (params->UdX)[j] = sum;
  }
  return NULL;
}

void *cmpXnorm(void *args){
  /* unstack parameters */
  const cmpUdX_args *params = (cmpUdX_args *)args;
  const double **feature = params->feature;
  const unsigned int *h_i = params->data->i;
  const unsigned int *h_j = params->data->j;
  const unsigned int *kmer1 = params->ckps->kmer1;
  const unsigned int *kmer2 = params->ckps->kmer2;
  const unsigned int *revcmp1 = params->ckps->revcmp1;
  const unsigned int *revcmp2 = params->ckps->revcmp2;

  /* compute the dot product between U and X^{(j)} */
  unsigned int i, j;
  double sum, pf;
  for(j = params->begin; j < params->end; j++){
    sum = 0;
    for(i = 0; i < params->n; i++){
      pf =  ((feature[h_i[i]][kmer1[j]] *
	      feature[h_j[i]][kmer2[j]]) +
	     (feature[h_i[i]][revcmp1[j]] *
	      feature[h_j[i]][revcmp2[j]]));
      sum += pf * pf;     
    }
    (params->Xnorm)[j] = sum;
  }
  return NULL;
}

int pthread_prep(const int thread_num,
		 const unsigned long n,
		 const unsigned long p,		 
		 const double **feature,
		 const hic *data,
		 const canonical_kp *ckps,
		 double *U,
		 double *UdX,
		 double *Xnorm,
		 cmpUdX_args **params,
		 pthread_t **threads){
  int i = 0;

  *params = calloc_errchk(thread_num,			   
			  sizeof(cmpUdX_args),
			  "calloc: cmpUdX_args[]");
  *threads = calloc_errchk(thread_num,			   
			   sizeof(pthread_t),
			   "calloc: threads[]");        
  /* set variables */
  for(i = 0; i < thread_num; i++){
    (*params)[i].thread_id = i;
    (*params)[i].begin = ((i == 0) ? 0 : (*params)[i - 1].end);
    (*params)[i].end =
      ((i == (thread_num - 1)) ? p : (p / thread_num) * (i + 1));
    (*params)[i].n = n;
    (*params)[i].feature = feature;
    (*params)[i].data = data;
    (*params)[i].ckps = ckps;
    (*params)[i].U     = U;
    (*params)[i].UdX   = UdX;
    (*params)[i].Xnorm = Xnorm;
  }
  return 0;
}

int l2boost_train(const cmd_args *args,
		  const double **feature,
		  const hic *data,
		  const canonical_kp *ckps,
		  const unsigned int iternum,
		  l2boost **model,
		  FILE *fp){
  
  const unsigned long n = data->nrow;
  const unsigned long p = ckps->num;
  const int thread_num = args->thread_num;
  unsigned long s = 0;
  double gamma = 0;
  double *U, *UdX, *Xnorm;
  unsigned int m = 0;

  /* allocate memory */
  {
    *model = calloc_errchk(1, sizeof(l2boost), "calloc l2boost");
    (*model)->res_sq = calloc_errchk(iternum, sizeof(double),
				     "calloc l2boost -> res_sq");
    (*model)->beta   = calloc_errchk(p, sizeof(double),
				     "calloc l2boost -> beta");
    U     = calloc_errchk(n, sizeof(double), "calloc U[]");
    UdX   = calloc_errchk(p, sizeof(double), "calloc UdX[]");
    Xnorm = calloc_errchk(p, sizeof(double), "calloc Xnorm[]");
  }

  if(thread_num >= 1){
    int t = 0;
    cmpUdX_args *params;
    pthread_t *threads;

    /* prepare for thread programming */
    pthread_prep(thread_num, n, p, 
		 feature, data, ckps,
		 U, UdX, Xnorm,
		 &params, &threads);

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start computation of Xnorm with %d threads\n",
	    thread_num);

    /* compute Xnorm ||X^{(j)}|| */
    {
      for(t = 0; t < thread_num; t++){
	pthread_create(&threads[t], NULL, 
		       cmpXnorm, (void*)&params[t]);		       
      } 
      for(t = 0; t < thread_num; t++){
	pthread_join(threads[t], NULL);
      }
    }

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "Xnorm finished\n");


  }
#if 0
    gettimeofday(&t0, NULL);

    /* AdaBoost iterations */
    for(t = 0; t < cmd_args->iteration_num; t++){
      /* step 1 : compute normalized weights p[] */
      {
	wsum = 0;
	for(n = 0; n < (hic->nrow); n++){
	  wsum += w[n];
	}
	for(n = 0; n < (hic->nrow); n++){
	  p[n] = 1.0 * w[n] / wsum;
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
	    ((kmer_freq[hic->i[n]][kp->l1[((*model)->axis)[t]]] * 
	      kmer_freq[hic->j[n]][kp->m1[((*model)->axis)[t]]] +
	      kmer_freq[hic->i[n]][kp->l2[((*model)->axis)[t]]] * 
	      kmer_freq[hic->j[n]][kp->m2[((*model)->axis)[t]]]) > 0) ? 1 : 0;
	  if(((((*model)->sign)[t] == 0) && pred == y[n]) ||
	     ((((*model)->sign)[t] == 1) && pred != y[n])){
	    w[n] *= ((*model)->beta)[t];
	  }
	}
      }
      gettimeofday(&time, NULL);
      adaboost_show_itr(stderr, 
			*model, (const char**)kmer_strings, kp, 
			t, diffSec(t0, time));
    }
  }
  
  /* write to file OR stderr */
  {
    if(output_file == NULL){
      adaboost_show_all(stderr, *model, (const char**)kmer_strings, kp);
    }else{
      FILE *fp;
      if((fp = fopen(output_file, "w")) == NULL){
	fprintf(stderr, "error: fopen %s\n%s\n",
		output_file, strerror(errno));
	exit(EXIT_FAILURE);
      }
      fprintf(stderr, "%s: info: AdaBoost: writing results to file: %s\n",
	      cmd_args->prog_name, output_file);
      adaboost_show_all(fp, *model, (const char**)kmer_strings, kp);
      fclose(fp);
    }
  }
#endif
  return 0;
}




#if 0

int adaboost_show_itr(FILE *fp, 
		      const adaboost *model,
		      const char **kmer_strings,
		      const canonical_kp *kp,
		      const unsigned long t,
		      double time);

int adaboost_show_all(FILE *fp, 
		      const adaboost *model,
		      const char **kmer_strings,
		      const canonical_kp *kp);

void *adaboost_comp_err(void *args);

int adaboost_set_y(hic *hic,
		   const double threshold,
		   unsigned int **y);


int adaboost_learn(const command_line_arguements *cmd_args,
		   const unsigned int **kmer_freq,
		   hic *hic,
		   const double threshold,
		   const canonical_kp *kp,
		   adaboost **model,
		   const char *output_file);
#endif

#if 0
int adaboost_show_itr(FILE *fp, 
		      const adaboost *model,
		      const char **kmer_strings,
		      const canonical_kp *kp,
		      const unsigned long t,
		      double time){
  fprintf(fp, "%ld\t%e\t%d\t%ld\t%s\t%s\t%s\t%s\t%f\t%f\n",
	  t, 
	  (model->beta)[t],
	  (model->sign)[t],
	  (model->axis)[t],
	  kmer_strings[kp->l1[(model->axis)[t]]],
	  kmer_strings[kp->m1[(model->axis)[t]]],
	  kmer_strings[kp->l2[(model->axis)[t]]],
	  kmer_strings[kp->m2[(model->axis)[t]]],
	  time,
	  time / (t + 1));
  return 0;
}

int adaboost_show_all(FILE *fp, 
		      const adaboost *model,
		      const char **kmer_strings,
		      const canonical_kp *kp){
  unsigned long t;
  for(t = 0; t < model->T; t++){
    fprintf(fp, "%ld\t%e\t%d\t%ld\t%s\t%s\t%s\t%s\n",
	    t, 
	    (model->beta)[t],
	    (model->sign)[t],
	    (model->axis)[t],
	    kmer_strings[kp->l1[(model->axis)[t]]],
	    kmer_strings[kp->m1[(model->axis)[t]]],
	    kmer_strings[kp->l2[(model->axis)[t]]],
	    kmer_strings[kp->m2[(model->axis)[t]]]);
  }
  return 0;
}
#endif		      

#if 0
void *adaboost_comp_err(void *args){
  adaboost_comp_err_args *params = (adaboost_comp_err_args *)args;
  unsigned int kmerpair = 0, x = 0, pred;
 
  for(kmerpair = params->begin; kmerpair <= params->end; kmerpair++){
    (*(params->err))[kmerpair] = 0;
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
	  (*(params->err))[kmerpair] += (params->p)[x];
	}
      }
    }
  }  
  return NULL;
}

int adaboost_set_y(hic *hic,
		   const double threshold,
		   unsigned int **y){
  unsigned int x;
  *y = calloc_errchk(hic->nrow, sizeof(unsigned int), "calloc: y");
  for(x = 0; x < hic->nrow; x++){
    (*y)[x] = (hic->mij)[x] > threshold ? 1 : 0; 
  }
  return 0;
}

int adaboost_learn(const command_line_arguements *cmd_args,
		   const unsigned int **kmer_freq,
		   hic *hic,
		   const double threshold,
		   const canonical_kp *kp,
		   adaboost **model,
		   const char *output_file){
  const unsigned long canonical_kmer_pair_num = 
    (1 << (4 * (cmd_args->k) - 1)) + (1 << (2 * (cmd_args->k) - 1));  
  unsigned long n, lm, argmin_lm, argmax_lm;
  unsigned int *marked, *y, pred;
  double *err, *w, *p, wsum, epsilon, min, max;
  char **kmer_strings;
  struct timeval t0, time;

  /* allocate memory */
  {
    *model = calloc_errchk(1, sizeof(adaboost), "calloc adaboost");
    (*model)->axis = calloc_errchk(cmd_args->iteration_num, sizeof(unsigned long), "calloc adaboost -> axis");
    (*model)->beta = calloc_errchk(cmd_args->iteration_num, sizeof(double), "calloc adaboost -> beta");
    (*model)->sign = calloc_errchk(cmd_args->iteration_num, sizeof(unsigned int), "calloc adaboost -> sign");
    (*model)->T = cmd_args->iteration_num;
    marked = calloc_errchk(canonical_kmer_pair_num, sizeof(unsigned int), "calloc: marked");
    err = calloc_errchk(canonical_kmer_pair_num, sizeof(double), "calloc: err");
    w = calloc_errchk(hic->nrow, sizeof(double), "calloc: p");
    p = calloc_errchk(hic->nrow, sizeof(double), "calloc: p");
    for(n = 0; n < hic->nrow; n++){
      w[n] = 1.0 / (hic->nrow);
    }
    adaboost_set_y(hic, threshold, &y);
    set_kmer_strings(cmd_args->k, &kmer_strings);
  }

  /* mark k-mer pairs containing forbidden k-mers */
  {
    unsigned int mask, m;
    unsigned int forbidden_kmer[1] = {141}; /* GATC */
    unsigned int forbidden_kmer_len[1] = {4};
    unsigned int forbidden_kmer_num = 1;
    for(lm = 0; lm < canonical_kmer_pair_num; lm++){
      for(m = 0; m < forbidden_kmer_num; m++){
	if(marked[lm] != 0){
	  break;
	}else{
	  mask = (1 << (2 * forbidden_kmer_len[m])) - 1;
	  for(n = 0; 
	      n < cmd_args->k - forbidden_kmer_len[m] + 1;
	      n++){
	    if(      (((kp->l1[lm]) >> (2 * n)) & mask) == forbidden_kmer[m]){
	      marked[lm] = 1;
	      break;
	    }else if((((kp->m1[lm]) >> (2 * n)) & mask) == forbidden_kmer[m]){
	      marked[lm] = 1;
	      break;
	    }else if((((kp->l2[lm]) >> (2 * n)) & mask) == forbidden_kmer[m]){
	      marked[lm] = 1;
	      break;
	    }else if((((kp->m2[lm]) >> (2 * n)) & mask) == forbidden_kmer[m]){	
	      marked[lm] = 1;
	      break;
	    }
	  }
	}
      }
    }
  }

  {
    n = 0;
    for(lm = 0; lm < canonical_kmer_pair_num; lm++){
      if(marked[lm] != 0){
	n++;
      }
    }
    fprintf(stderr, "%s: info: AdaBoost: %ld out of %ld k-mer pairs are filtered out\n",
	    cmd_args->prog_name, n, canonical_kmer_pair_num);
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
	params[i].l1 = kp->l1;
	params[i].m1 = kp->m1;
	params[i].l2 = kp->l2;
	params[i].m2 = kp->m2;
	params[i].err = &err;
	params[i].p = p;
	params[i].y = y;
      }
    }

    gettimeofday(&t0, NULL);

    /* AdaBoost iterations */
    for(t = 0; t < cmd_args->iteration_num; t++){
      /* step 1 : compute normalized weights p[] */
      {
	wsum = 0;
	for(n = 0; n < (hic->nrow); n++){
	  wsum += w[n];
	}
	for(n = 0; n < (hic->nrow); n++){
	  p[n] = 1.0 * w[n] / wsum;
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
	    ((kmer_freq[hic->i[n]][kp->l1[((*model)->axis)[t]]] * 
	      kmer_freq[hic->j[n]][kp->m1[((*model)->axis)[t]]] +
	      kmer_freq[hic->i[n]][kp->l2[((*model)->axis)[t]]] * 
	      kmer_freq[hic->j[n]][kp->m2[((*model)->axis)[t]]]) > 0) ? 1 : 0;
	  if(((((*model)->sign)[t] == 0) && pred == y[n]) ||
	     ((((*model)->sign)[t] == 1) && pred != y[n])){
	    w[n] *= ((*model)->beta)[t];
	  }
	}
      }
      gettimeofday(&time, NULL);
      adaboost_show_itr(stderr, 
			*model, (const char**)kmer_strings, kp, 
			t, diffSec(t0, time));
    }
  }
  
  /* write to file OR stderr */
  {
    if(output_file == NULL){
      adaboost_show_all(stderr, *model, (const char**)kmer_strings, kp);
    }else{
      FILE *fp;
      if((fp = fopen(output_file, "w")) == NULL){
	fprintf(stderr, "error: fopen %s\n%s\n",
		output_file, strerror(errno));
	exit(EXIT_FAILURE);
      }
      fprintf(stderr, "%s: info: AdaBoost: writing results to file: %s\n",
	      cmd_args->prog_name, output_file);
      adaboost_show_all(fp, *model, (const char**)kmer_strings, kp);
      fclose(fp);
    }

  }
  return 0;
}
#endif
#endif
