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
  unsigned int nextiter;
  unsigned int iternum;
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
      /* pf : pairwise feature */
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

int dump_model(const l2boost *model, 
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

unsigned long l2_select_axis(const double *UdX, 
			     const double *Xnorm,
			     const unsigned long p){
  unsigned long argmax = 0;
  double max = UdX[argmax] / Xnorm[argmax];
  unsigned long j;
  for(j = 1; j < p; j++){
    if(max < UdX[j] / Xnorm[j]){
      argmax = j;
      max = UdX[argmax] / Xnorm[argmax];
    }
  }
  return argmax;
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
  residual_square[m] = 0;
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


int l2boost_step_dump(const l2boost *model,
		      const unsigned int m,
		      const unsigned long s,
		      const double gamma,
		      const double sec_per_step,
		      const double sec_total,
		      const char *prog_name,
		      FILE *fp){
  fprintf(fp, "%s [INFO] \t ", prog_name);
  fprintf(fp, "%d\t%ld\t%e\t%e\t%f\t%f\n",
	  m, s, gamma, (model->res_sq)[m], sec_per_step, sec_total);
  return 0;
}

int l2boost_init(const canonical_kp *ckps,		 
		 const unsigned int iternum,
		 l2boost **model){  
  const unsigned long p = ckps->num;

  /* allocate memory */
  {
    *model = calloc_errchk(1, sizeof(l2boost), "calloc l2boost");
    (*model)->res_sq = calloc_errchk(iternum + 1, sizeof(double),
				     "calloc l2boost -> res_sq");
    (*model)->beta   = calloc_errchk(p, sizeof(double),
				     "calloc l2boost -> beta");
    (*model)->iternum = iternum;
    (*model)->nextiter = 1;
  }

  return 0;
}

int l2boost_train(const cmd_args *args,
		  const double **feature,
		  const hic *data,
		  const canonical_kp *ckps,
		  const double v,
		  l2boost **model,
		  FILE *fp){  
  const unsigned long n = data->nrow;
  const unsigned long p = ckps->num;
  const int thread_num = args->thread_num;
  unsigned long s = 0;
  double gamma = 0;
  double *U, *UdX, *Xnorm;
  unsigned int m = 0;
  struct timeval time_start, time_prev, time;

  /* allocate memory */
  {
    U     = calloc_errchk(n, sizeof(double), "calloc U[]");
    UdX   = calloc_errchk(p, sizeof(double), "calloc UdX[]");
    Xnorm = calloc_errchk(p, sizeof(double), "calloc Xnorm[]");
  }

  /* initialize residuals U[] := Y[] and 
   * compute \sum_i U[i]^2                */
  {
    const double *Y = data->mij;
    unsigned long i;
    ((*model)->res_sq)[0] = 0;
    for(i = 0; i < n; i++){
      U[i] = 1.0 * Y[i];
      ((*model)->res_sq)[0] += 1.0 * U[i] * U[i] / n;
    }
  }

  if(thread_num >= 1){
    int t = 0;
    cmpUdX_args *params;
    pthread_t *threads;

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start computation of Xnorm with %d threads\n",
	    thread_num);

    /* compute Xnorm ||X^{(j)}|| */
    {
      gettimeofday(&time_prev, NULL);

      /* prepare for thread programming */
      pthread_prep(thread_num, n, p, 
		   feature, data, ckps,
		   U, UdX, Xnorm, 
		   &params, &threads);
		   
      for(t = 0; t < thread_num; t++){
	pthread_create(&threads[t], NULL, 
		       cmpXnorm, (void*)&params[t]);		       
      } 
      for(t = 0; t < thread_num; t++){
	pthread_join(threads[t], NULL);
      }
      gettimeofday(&time, NULL);
    }

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "Xnorm finished in %f sec.\n", diffSec(time_prev, time));
    cpTimeval(time, &time_prev);
    cpTimeval(time, &time_start);

    fprintf(fp, "%s [INFO] \t ", args->prog_name);
    fprintf(fp, "iter \t axis \t gamma \t residuals \t step t \t total t\n");

    l2boost_step_dump(*model,
		      (const unsigned int)m, 
		      (const unsigned long)s,
		      (const double)gamma,
		      (const double)diffSec(time_prev, time),
		      (const double)diffSec(time_start, time),
		      (const char *)args->prog_name,
		      fp);
    

    for(m = (*model)->nextiter; m <= (*model)->iternum; m++){
      /* compute inner product $U \cdot X^{(j)}$ */
      {
	/* prepare for thread programming */
	pthread_prep(thread_num, n, p, 
		     feature, data, ckps,
		     U, UdX, Xnorm, 
		     &params, &threads);
		   
	for(t = 0; t < thread_num; t++){
	  pthread_create(&threads[t], NULL, 
			 cmpUdX, (void*)&params[t]);		       
	} 
	for(t = 0; t < thread_num; t++){
	  pthread_join(threads[t], NULL);
	}
      }

      /* select axis */
      s = l2_select_axis((const double *)UdX, (const double *)Xnorm, p);			 
      gamma = UdX[s] / Xnorm[s];
      
      ((*model)->beta)[s] += v * gamma;

      /* Update U[] and sum of residual square */
      l2_update_U(U, (*model)->res_sq, 
		  (const unsigned int)m, n, s, 
		  feature, data, ckps,
		  (const double)gamma, v);
      
      gettimeofday(&time, NULL);
      l2boost_step_dump(*model,
			(const unsigned int)m, 
			(const unsigned long)s,
			(const double)gamma,
			(const double)diffSec(time_prev, time),
			(const double)diffSec(time_start, time),
			(const char *)args->prog_name,
			fp);
      cpTimeval(time, &time_prev);
    }

  }

  {
    dump_model((const l2boost *)*model, p);
  }
  return 0;
}

#endif
