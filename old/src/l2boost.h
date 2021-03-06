#ifndef __l2boost_H__
#define __l2boost_H__ 

#include <sys/time.h>
#include <math.h>
#include <pthread.h>
#include "calloc_errchk.h"
#include "diffSec.h"
#include "kmer.h"
#include "hic.h"

/* boost results */
typedef struct _boost{
  double *res_sq;
  double *beta;
  unsigned int nextiter;
  unsigned int iternum;
} boost;

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
  const kmer *kmers;
  /* shared data */
  double *U;
  double *UdX;
  double *Xnormsq;
} cmpUdX_args;

void *boost_cmpXnormsq(void *args);
int boost_dump_beta(const boost *model, 
		    const unsigned long p);
int boost_pthread_prep(const int thread_num,
		       const unsigned long n,
		       const unsigned long p,		 
		       const double **feature,
		       const hic *data,
		       const canonical_kp *ckps,
		       const kmer *kmers,
		       double *U,
		       double *UdX,
		       double *Xnormsq,
		       cmpUdX_args **params,
		       pthread_t **threads);
unsigned long boost_select_axis(const double *UdX, 
				const double *Xnormsq,
				const unsigned long p);
int boost_step_dump_head(FILE *fp);
int boost_step_dump(const boost *model,
		    const unsigned int m,
		    const unsigned long s,
		    const double v_gamma,
		    const double sec_per_step,
		    const double sec_total,
		    FILE *fp);
int boost_init(const cmd_args *args,
	       const canonical_kp *ckps,
	       const kmer *kmers,
	       const unsigned int iternum,
	       const char *file,
	       boost **model,
	       FILE *fp_out);

void *l2_cmpUdX(void *args);
int l2_update_U(double *U, 
		double *residual_square,
		const unsigned int m,
		const unsigned long n,
		const unsigned long s,
		const double **feature,
		const hic *data,
		const canonical_kp *ckps,
		const double gamma, 
		const double v);
int l2_train(const cmd_args *args,
	     const double **feature,
	     const hic *data,
	     const canonical_kp *ckps,
	     const double v,
	     boost **model,
	     FILE *fp);


/**
 * Boosting common functions
 **/

void *boost_cmpXnormsq(void *args){
  /* unstack parameters */
  const cmpUdX_args *params = (cmpUdX_args *)args;
  const double **feature = params->feature;
  const unsigned int *h_i = params->data->i;
  const unsigned int *h_j = params->data->j;
  const unsigned int *kmer1 = params->ckps->kmer1;
  const unsigned int *kmer2 = params->ckps->kmer2;
  const unsigned int *revcmp1 = params->ckps->revcmp1;
  const unsigned int *revcmp2 = params->ckps->revcmp2;

  /* compute ||X^{(j)}||^2 */
  unsigned int i, j;
  double sum, pf;
  for(j = params->begin; j < params->end; j++){
    sum = 0;
    for(i = 0; i < params->n; i++){
      /* pf : pairwise feature */
      pf =  ((feature[h_i[i]][kmer1[j]] *
	      feature[h_j[i]][kmer2[j]]) +
	     (feature[h_i[i]][revcmp1[j]] *
	      feature[h_j[i]][revcmp2[j]]));
      sum += pf * pf;     
    }
    (params->Xnormsq)[j] = sum;
  }
  return NULL;
}

int boost_dump_beta(const boost *model, 
	       const unsigned long p){
  unsigned long j;
  for(j = 0; j < p; j++){
    if((model->beta)[j] != 0){
      fprintf(stderr, "%ld\t%e\n",
	      j, (model->beta)[j]);
    }
  }
  return 0;
}

int boost_pthread_prep(const int thread_num,
		 const unsigned long n,
		 const unsigned long p,		 
		 const double **feature,
		 const hic *data,
		 const canonical_kp *ckps,
		 const kmer *kmers,
		 double *U,
		 double *UdX,
		 double *Xnormsq,
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
    (*params)[i].data    = data;
    (*params)[i].ckps    = ckps;
    (*params)[i].kmers   = kmers;
    (*params)[i].U       = U;
    (*params)[i].UdX     = UdX;
    (*params)[i].Xnormsq = Xnormsq;
  }
  return 0;
}

unsigned long boost_select_axis(const double *UdX, 
			     const double *Xnormsq,
			     const unsigned long p){
  unsigned long argmax = 0;
  double max = UdX[argmax] * UdX[argmax] / Xnormsq[argmax];
  unsigned long j;
  for(j = 1; j < p; j++){
    if(max < UdX[j] * UdX[j] / Xnormsq[j]){
      argmax = j;
      max = UdX[argmax] * UdX[argmax] / Xnormsq[argmax];
    }
  }
  return argmax;
}

int boost_step_dump_head(FILE *fp){
  fprintf(fp, "iter \t axis \t gamma \t residuals \t step t \t total t\n");
  return 0;
}

int boost_step_dump(const boost *model,
		 const unsigned int m,
		 const unsigned long s,
		 const double v_gamma,
		 const double sec_per_step,
		 const double sec_total,
		 FILE *fp){
  fprintf(fp, "%d\t%ld\t%e\t%e\t%f\t%f\n",
	  m, s, v_gamma, (model->res_sq)[m], sec_per_step, sec_total);
  fflush(fp);
  return 0;
}

int boost_init(const cmd_args *args,
	       const canonical_kp *ckps,
	       const kmer *kmers,
	       const unsigned int iternum,
	       const char *file,
	       boost **model,
	       FILE *fp_out){  
  unsigned long p;
  if(ckps != NULL){
    p = ckps->num;
  }else if(kmers != NULL){
    p = kmers->num;
  }else{
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "ckps and kmers are NULL\n");
    exit(EXIT_FAILURE);   
  }

  /* allocate memory */
  {
    *model = calloc_errchk(1, sizeof(boost), "calloc boost");
    (*model)->res_sq = calloc_errchk(iternum + 1, sizeof(double),
				     "calloc boost -> res_sq");
    (*model)->beta   = calloc_errchk(p, sizeof(double),
				     "calloc boost -> beta");
    (*model)->iternum = iternum;
    (*model)->nextiter = 1;
  }

  boost_step_dump_head(fp_out);

  if(file == NULL){
    fprintf(fp_out, "%d\t%ld\t%e\t%e\t%f\t%f\n",
	    0, (long int)0, 0.0, 1.0, 0.0, 0.0);
  }else{
    /* read from a file */
    {
      FILE *fp;
      char buf[BUF_SIZE];
      char gamma_str[BUF_SIZE], residuals_str[BUF_SIZE];
      char step_t[BUF_SIZE], total_t[BUF_SIZE];
      unsigned int iter;
      unsigned long axis;
      unsigned int m, model_len;
      
      model_len = mywc(file) - 2;
      
      if(model_len >= iternum){
	fprintf(stderr, "%s [ERROR] ", args->prog_name);
	fprintf(stderr, "iternum (%d) is smaller than saved file (%d)\n",
		iternum, model_len);
	exit(EXIT_FAILURE);
      }
      
      fprintf(stderr, "%s [INFO] ", args->prog_name);
      fprintf(stderr, "start reading model from %s\n", file);
      
      if((fp = fopen(file, "r")) == NULL){
	fprintf(stderr, "error: fopen %s\n%s\n",
		file, strerror(errno));
	exit(EXIT_FAILURE);
      }
      
      fprintf(stderr, "%s [INFO] \t ", args->prog_name);
      boost_step_dump_head(stderr);
      
      /* skip a header line */
      fgets(buf, BUF_SIZE, fp);
      
      m = 0;
      while(fgets(buf, BUF_SIZE, fp) && m <= model_len){
	/* parse the input line */
	sscanf(buf, "%d\t%ld\t%s\t%s\t%s\t%s", 
	       &iter, &axis, gamma_str, residuals_str,
	       step_t, total_t);
	
	if(iter != m){
	  fprintf(stderr, "%s [ERROR] ", args->prog_name);
	  fprintf(stderr, "iteration step %d is missing\n", m);
	  exit(EXIT_FAILURE);
	}
	
	(*model)->beta[axis] += strtod(gamma_str, NULL);
	(*model)->res_sq[m] = strtod(residuals_str, NULL);      
	
	boost_step_dump(*model, m, axis, 
		     strtod(gamma_str, NULL), 0, 0,
		     fp_out);
	fprintf(stderr, "%s [INFO] \t ", args->prog_name);
	boost_step_dump(*model, m, axis, 
		     strtod(gamma_str, NULL), 0, 0,
		     stderr);
	
	
	(*model)->nextiter = ++m;      
      }
      
      fclose(fp);
      
      fprintf(stderr, "%s [INFO] ", args->prog_name);
      fprintf(stderr, "%d iterations has loaded from file\n", m - 1);
    }
  }
  return 0;
}


/**
 * L2 Boosting 
 **/

void *l2_cmpUdX(void *args){
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

int l2_update_U(double *U, 
		double *residual_square,
		const unsigned int m,
		const unsigned long n,
		const unsigned long s,
		const double **feature,
		const hic *data,
		const canonical_kp *ckps,
		const double gamma, 
		const double v){
  const unsigned int *h_i = data->i;
  const unsigned int *h_j = data->j;
  const unsigned int *kmer1 = ckps->kmer1;
  const unsigned int *kmer2 = ckps->kmer2;
  const unsigned int *revcmp1 = ckps->revcmp1;
  const unsigned int *revcmp2 = ckps->revcmp2;
  unsigned long i;
  double sum = 0, pf = 0;
  for(i = 0; i < n; i++){
    /* pf : pairwise feature */
    pf =  ((feature[h_i[i]][kmer1[s]] *
	    feature[h_j[i]][kmer2[s]]) +
	   (feature[h_i[i]][revcmp1[s]] *
	    feature[h_j[i]][revcmp2[s]]));
    U[i] -= v * gamma * pf;
    sum += 1.0 * U[i] * U[i] / n;
  }
  residual_square[m] = sum;
  return 0;
}
			  
int l2_train(const cmd_args *args,
		  const double **feature,
		  const hic *data,
		  const canonical_kp *ckps,
		  const double v,
		  boost **model,
		  FILE *fp){  
  const unsigned long n = data->nrow;
  const unsigned long p = ckps->num;
  const int thread_num = args->thread_num;
  unsigned long s = 0;
  double gamma = 0;
  double *U, *UdX, *Xnormsq;
  unsigned int m = 0;
  struct timeval time_start, time_prev, time;

  /* allocate memory */
  {
    U       = calloc_errchk(n, sizeof(double), "calloc U[]");
    UdX     = calloc_errchk(p, sizeof(double), "calloc UdX[]");
    Xnormsq = calloc_errchk(p, sizeof(double), "calloc Xnormsq[]");
  }

  /* initialize residuals U[] := Y[] and 
   * compute \sum_i U[i]^2                */
  {
    const double *Y = data->mij;
    unsigned long i;
    ((*model)->res_sq)[0] = 0;
    for(i = 0; i < n; i++){
      U[i] = Y[i];
      ((*model)->res_sq)[0] += U[i] * U[i] / n;
    }
  }

  /* If we load some model from a file, update residuals */
  if(((*model)->nextiter) > 1){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start computation of residuals\n");
    const unsigned int *h_i = data->i;
    const unsigned int *h_j = data->j;
    const unsigned int *kmer1 = ckps->kmer1;
    const unsigned int *kmer2 = ckps->kmer2;
    const unsigned int *revcmp1 = ckps->revcmp1;
    const unsigned int *revcmp2 = ckps->revcmp2;
    unsigned long i, j;
    double pf;
    /* pf : pairwise feature */
    for(j = 0; j < p; j++){
      if(((*model)->beta[j]) != 0){
	for(i = 0; i < n; i++){
	  pf =  ((feature[h_i[i]][kmer1[j]] *
		  feature[h_j[i]][kmer2[j]]) +
		 (feature[h_i[i]][revcmp1[j]] *
		  feature[h_j[i]][revcmp2[j]]));
	  U[i] -= ((*model)->beta[j]) * pf;
	}
      }
    }
  }

  if(thread_num >= 1){
    int t = 0;
    cmpUdX_args *params;
    pthread_t *threads;

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start computation of Xnormsq with %d threads\n",
	    thread_num);

    /* compute Xnormsq ||X^{(j)}||^2 */
    {
      gettimeofday(&time_prev, NULL);

      /* prepare for thread programming */
      boost_pthread_prep(thread_num, n, p, 
			 feature, data, ckps, NULL,
			 U, UdX, Xnormsq, 
			 &params, &threads);
		   
      for(t = 0; t < thread_num; t++){
	pthread_create(&threads[t], NULL, 
		       boost_cmpXnormsq, (void*)&params[t]);		       
      } 
      for(t = 0; t < thread_num; t++){
	pthread_join(threads[t], NULL);
      }
      gettimeofday(&time, NULL);
    }

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "Xnormsq finished in %f sec.\n", diffSec(time_prev, time));

#if 0
    {
      unsigned int tmp;
      for(tmp = 0; tmp < p; tmp++){
	fprintf(stderr, "%e\t", Xnormsq[tmp]);
      }
      fprintf(stderr, "\n");
    }
#endif

    cpTimeval(time, &time_prev);
    cpTimeval(time, &time_start);
    
    for(m = (*model)->nextiter; m <= (*model)->iternum; m++){
      /* compute inner product $U \cdot X^{(j)}$ */
      {
	/* prepare for thread programming */
	boost_pthread_prep(thread_num, n, p, 
			   feature, data, ckps, NULL,
			   U, UdX, Xnormsq, 
			   &params, &threads);
		   
	for(t = 0; t < thread_num; t++){
	  pthread_create(&threads[t], NULL, 
			 l2_cmpUdX, (void*)&params[t]);		       
	} 
	for(t = 0; t < thread_num; t++){
	  pthread_join(threads[t], NULL);
	}
      }

      /* select axis */
      s = boost_select_axis((const double *)UdX, (const double *)Xnormsq, p);			 
      gamma = UdX[s] / Xnormsq[s];
      
      ((*model)->beta)[s] += v * gamma;

      /* Update U[] and sum of residual square */
      l2_update_U(U, (*model)->res_sq, 
		  (const unsigned int)m, n, s, 
		  feature, data, ckps,
		  (const double)gamma, v);
      
      gettimeofday(&time, NULL);
      boost_step_dump(*model,
			(const unsigned int)m, 
			(const unsigned long)s,
			v * (const double)gamma,
			(const double)diffSec(time_prev, time),
			(const double)diffSec(time_start, time),
			fp);
      fprintf(stderr, "%s [INFO] \t ", args->prog_name);
      boost_step_dump(*model,
			(const unsigned int)m, 
			(const unsigned long)s,
			v * (const double)gamma,
			(const double)diffSec(time_prev, time),
			(const double)diffSec(time_start, time),
			stderr);
      cpTimeval(time, &time_prev);
    }

  }

  {
    boost_dump_beta((const boost *)*model, p);
  }
  return 0;
}


/**
 * AdaBoost
 **/

void *ada_cmpUdX(void *args){
  /* unstack parameters */
  const cmpUdX_args *params = (cmpUdX_args *)args;
  const double **feature = params->feature;
  const unsigned int *h_i = params->data->i;
  const unsigned int *h_j = params->data->j;
  const double *Y = params->data->mij;
  const double *beta_x = params->U;
  const unsigned int *kmer = params->kmers->kmer1;

  /* compute the dot product between U and X^{(j)} */
  unsigned int i, j;
  unsigned long TT, TF, FT, FF;
  for(j = params->begin; j < params->end; j++){
    TT = TF = FT = FF = 0;
    for(i = 0; i < params->n; i++){
      /* (i, m^{i,j}) */
      if(((beta_x[2 * i] * Y[i]) > 0) || 
	 ((beta_x[2 * i] == 0) && (Y[i] >= 0))){
	if(((feature[h_i[i]][kmer[j]] * Y[i]) > 0) || 
	   ((feature[h_i[i]][kmer[j]] == 0) && (Y[i] >= 0))){
	  TT += 1;
	}else{
	  TF += 1;
	}
      }else{
	if(((feature[h_i[i]][kmer[j]] * Y[i]) > 0) || 
	   ((feature[h_i[i]][kmer[j]] == 0) && (Y[i] >= 0))){
	  FT += 1;
	}else{
	  FF += 1;
	}
      }      
      
      /* (j, m^{i,j}) */
      if(((beta_x[2 * i + 1] * Y[i]) > 0) || 
	 ((beta_x[2 * i + 1] == 0) && (Y[i] >= 0))){
	if(((feature[h_j[i]][kmer[j]] * Y[i]) > 0) || 
	   ((feature[h_j[i]][kmer[j]] == 0) && (Y[i] >= 0))){
	  TT += 1;
	}else{
	  TF += 1;
	}
      }else{
	if(((feature[h_j[i]][kmer[j]] * Y[i]) > 0) || 
	   ((feature[h_j[i]][kmer[j]] == 0) && (Y[i] >= 0))){
	  FT += 1;
	}else{
	  FF += 1;
	}
      }      
    }
    (params->UdX)[j] = 1.0 * (TT - TF) * exp(-1) + 1.0 * (FT- FF) * exp(1);
  }
  return NULL;
}

int ada_update_beta_x(double *beta_x, 
		      double *residual_square,
		      const unsigned int m,
		      const unsigned long n,
		      const unsigned long s,
		      const double **feature,
		      const hic *data,
		      const kmer *kmers,
		      const double gamma, 
		      const double v){
  const unsigned int *h_i = data->i;
  const unsigned int *h_j = data->j;
  const double *Y = data->mij;
  const unsigned int *kmer = kmers->kmer1;

  unsigned long i = 0;
  unsigned long err = 0;

  /* update beta_x */
  for(i = 0; i < n; i++){
    /* (i, m^{i,j}) */
    if((feature[h_i[i]][kmer[s]]) >= 0){
      beta_x[2 * i] += v * gamma;
    }else{
      beta_x[2 * i] -= v * gamma;
    }
    /* (j, m^{i,j}) */
    if((feature[h_j[i]][kmer[s]]) >= 0){
      beta_x[2 * i + 1] += v * gamma;
    }else{
      beta_x[2 * i + 1] -= v * gamma;
    }
  }

  /* count error rate */
  for(i = 0; i < n; i++){
    if(((beta_x[2 * i] * Y[i]) < 0) ||
       ((beta_x[2 * i] == 0) && (Y[i] < 0))){
      err += 1;
    }
    if(((beta_x[2 * i + 1] * Y[i]) < 0) ||
       ((beta_x[2 * i + 1] == 0) && (Y[i] < 0))){
      err += 1;
    }
  }

  residual_square[m] = 0.5 * err / n;
  return 0;
}

#if 1
int ada_train(const cmd_args *args,
	      const double **feature,
	      const hic *data,
	      const kmer *kmers,
	      const double v,
	      boost **model,
	      FILE *fp){  
  const unsigned long n = data->nrow;
  const unsigned long p = kmers->num;
  const int thread_num = args->thread_num;
  unsigned long s = 0;
  double gamma = 0;
  double *beta_x, *UdX, *Xnormsq;
  unsigned int m = 0;
  struct timeval time_start, time_prev, time;

  /* allocate memory */
  {
    beta_x  = calloc_errchk(2 * n, sizeof(double), "calloc beta_x[]");
    UdX     = calloc_errchk(p, sizeof(double), "calloc UdX[]");
    Xnormsq = calloc_errchk(p, sizeof(double), "calloc Xnormsq[]");
  }

  /* count error rate */
  {
    const double *Y = data->mij;
    unsigned long i = 0, err = 0;    
    for(i = 0; i < n; i++){
      if(((beta_x[2 * i] * Y[i]) < 0) ||
	 ((beta_x[2 * i] == 0) && (Y[i] < 0))){
	err += 1;
      }
      if(((beta_x[2 * i + 1] * Y[i]) < 0) ||
	 ((beta_x[2 * i + 1] == 0) && (Y[i] < 0))){
	err += 1;
      }
    }
  
    ((*model)->res_sq)[(*model)->nextiter - 1] = 0.5 * err / n;
  }
  
  /* If we load some model from a file, update residuals */
  if(((*model)->nextiter) > 1){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start computation of residuals\n");
    const unsigned int *h_i = data->i;
    const unsigned int *h_j = data->j;
    const unsigned int *kmer = kmers->kmer1;
    unsigned long i, j;
    for(j = 0; j < p; j++){
      if(((*model)->beta[j]) != 0){
	for(i = 0; i < n; i++){
	  beta_x[2 * i]     = feature[h_i[i]][kmer[j]];
	  beta_x[2 * i + 1] = feature[h_j[i]][kmer[j]];
	}
      }
    }
  }

  /* prepare Xnormsq ||X^{(j)}||^2 */
  {
    unsigned long j;
    for(j = 0; j < p; j++){
      Xnormsq[j] = 2.0 * n;
    }

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "Xnormsq[] is uniform\n");
  }

#if 1
  if(thread_num >= 1){
    int t = 0;
    cmpUdX_args *params;
    pthread_t *threads;

    gettimeofday(&time, NULL);
    cpTimeval(time, &time_prev);
    cpTimeval(time, &time_start);
#if 1
    for(m = (*model)->nextiter; m <= (*model)->iternum; m++){
      /* compute inner product $U \cdot X^{(j)}$ */
      {
	/* prepare for thread programming */
	boost_pthread_prep(thread_num, n, p, 
			   feature, data, NULL, kmers,
			   beta_x, UdX, Xnormsq, 
			   &params, &threads);
		   
	for(t = 0; t < thread_num; t++){
	  pthread_create(&threads[t], NULL, 
			 ada_cmpUdX, (void*)&params[t]);		       
	} 
	for(t = 0; t < thread_num; t++){
	  pthread_join(threads[t], NULL);
	}
      }

      /* select axis */
      s = boost_select_axis((const double *)UdX, (const double *)Xnormsq, p);			 
      gamma = UdX[s] / Xnormsq[s];
      
      ((*model)->beta)[s] += v * gamma;
#if 1
      /* Update U[] and sum of residual square */
      ada_update_beta_x(beta_x, (*model)->res_sq, 
			(const unsigned int)m, n, s, 
			feature, data, kmers,
			(const double)gamma, v);
      
      gettimeofday(&time, NULL);
      boost_step_dump(*model,
		      (const unsigned int)m, 
		      (const unsigned long)s,
		      v * (const double)gamma,
		      (const double)diffSec(time_prev, time),
		      (const double)diffSec(time_start, time),
		      fp);
      fprintf(stderr, "%s [INFO] \t ", args->prog_name);
      boost_step_dump(*model,
		      (const unsigned int)m, 
		      (const unsigned long)s,
		      v * (const double)gamma,
		      (const double)diffSec(time_prev, time),
		      (const double)diffSec(time_start, time),
		      stderr);
      cpTimeval(time, &time_prev);
#endif
    }

#endif
  }

#endif
  {
    boost_dump_beta((const boost *)*model, p);
  }

  return 0;
}

#endif
#endif
