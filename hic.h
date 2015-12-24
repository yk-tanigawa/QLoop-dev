#ifndef __HIC_H__
#define __HIC_H__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "constant.h"
#include "cmd_args.h"
#include "mywc.h"
#include "calloc_errchk.h"
#include "io.h"
#include "diffSec.h"


/* normalized O/E converted Hi-C data */
typedef struct _hic {
  unsigned long nrow;
  unsigned int res;
  unsigned int *invalid;
  unsigned int *i;
  unsigned int *j;
  double *mij;
} hic;

/* Hi-C raw data */
typedef struct _hic_raw{
  hic *hic;
  double *norm;
  double *exp;
  unsigned long norm_len;
  unsigned long exp_len;
} hic_raw;

/* arguments for function hic_prep_thread */
typedef struct _hic_prep_thread_args{
  int thread_id;
  unsigned long begin;
  unsigned long end;
  unsigned long min_dist;
  unsigned long max_dist;
  unsigned int *h_invalid;
  unsigned int *h_i;
  unsigned int *h_j;
  double *h_mij;
  double *norm;
  double *exp;
} hic_prep_thread_args;

/**
 * read Hi-C data from a file 
 */

int hic_read(const char *hic_file,
	     const unsigned int res,
	     hic **data){

  /* allocate memory */
  *data = calloc_errchk(1, sizeof(hic), "calloc hic");
  {
    (*data)->nrow = mywc(hic_file);  
    (*data)->invalid = calloc_errchk((*data)->nrow, sizeof(unsigned int),
				     "calloc hic (*data)->invalid");
    (*data)->i = calloc_errchk((*data)->nrow, sizeof(unsigned int),
			       "calloc hic (*data)->i");
    (*data)->j = calloc_errchk((*data)->nrow, sizeof(unsigned int),
			       "calloc hic (*data)->j");
    (*data)->mij = calloc_errchk((*data)->nrow, sizeof(double), 
				 "calloc hic (*data)->mij");
    (*data)->res = res;
  }

  /* read from a file */
  {
    FILE *fp;
    char buf[256], tmp_mij_str[64];
    unsigned int tmp_i = 0, tmp_j = 0;
    unsigned long row = 0;

    if((fp = fopen(hic_file, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      hic_file, strerror(errno));
      exit(EXIT_FAILURE);
    }
    
    while(fgets(buf, BUF_SIZE, fp) && row < (*data)->nrow){
      sscanf(buf, "%d\t%d\t%s", &tmp_i, &tmp_j, (char *)(&tmp_mij_str));	 
      ((*data)->mij)[row] = strtod(tmp_mij_str, NULL);	  
      if(tmp_i <= tmp_j){
	((*data)->i)[row] = tmp_i / res;
	((*data)->j)[row] = tmp_j / res;
      }else{
	((*data)->i)[row] = tmp_j / res;
	((*data)->j)[row] = tmp_i / res;
      }
      row++;
    }
  
    fclose(fp);
  }

  return 0;
}


/**
 * read Hi-C raw data 
 *  - one raw data matrix
 *  - two vectors for normalization and O/E conversion
 */

inline char *res2str(const unsigned int res){
  switch(res){
    case 1000:
      return "1kb";
    default:
      fprintf(stderr, "resolution size %d is not supported\n", res);
     exit(EXIT_FAILURE);
  }
}

/* set appropriate file names */
inline void set_hic_file_names(const char *hicDir,
			       const unsigned int res,
			       const unsigned int chr,
			       const char *norm,
			       const char *exp,
			       char **hic_raw_file,
			       char **hic_norm_file,
			       char **hic_exp_file){
  {
    char file_head[F_NAME_LEN]; 
    sprintf(file_head, "%s/%s_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%s",
	    hicDir, res2str(res), chr, chr, res2str(res));
    
    *hic_raw_file = calloc_errchk(F_NAME_LEN, sizeof(char), "calloc: hic_raw_file");
    sprintf(*hic_raw_file, "%s.%s", file_head, "RAWobserved");

    if(norm != NULL){
      *hic_norm_file = calloc_errchk(F_NAME_LEN, sizeof(char), "calloc: hic_norm_file");
      sprintf(*hic_norm_file, "%s.%s%s", file_head, norm, "norm");
    }else{
      *hic_norm_file = NULL;
    }

    if(exp != NULL){
      *hic_exp_file = calloc_errchk(F_NAME_LEN, sizeof(char), "calloc: hic_exp_file");
      sprintf(*hic_exp_file, "%s.%s%s", file_head, exp, "expected");
    }else{
      *hic_exp_file = NULL;
    }
  }


  return;
}

/* read one matrix and two vectors */
int hic_raw_read(const char *hic_raw_dir,
		 const unsigned int res,
		 const unsigned int chr,
		 const char *norm,
		 const char *exp,
		 hic_raw **raw){

  char *hic_raw_file, *hic_norm_file, *hic_exp_file;    

  set_hic_file_names(hic_raw_dir, res, chr, norm, exp,
		     &hic_raw_file, 
		     &hic_norm_file,
		     &hic_exp_file);

  *raw = calloc_errchk(1, sizeof(hic_raw), "calloc hic_raw raw");

  if(hic_norm_file != NULL){
    read_double(hic_norm_file, &((*raw)->norm), &((*raw)->norm_len));
  }
  if(hic_exp_file != NULL){
    read_double(hic_exp_file, &((*raw)->exp), &((*raw)->exp_len));
  }

  hic_read(hic_raw_file, res, &((*raw)->hic));

  return 0;
}


/**
 * - normalization and/or O/E conversion
 * - size selection
 */

/* normalization and O/E conversion */
void *hic_prep_thread_norm_exp(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row, ij_dist;
 
  for(row = params-> begin; row <= params -> end; row++){
    ij_dist = (params->h_j)[row] - (params->h_i)[row];

    if(params->min_dist <= ij_dist && ij_dist <= params->max_dist){

      (params->h_mij)[row] /= ((params->norm)[(params->h_i)[row]] *
			       (params->norm)[(params->h_j)[row]] *
			       (params->exp)[(params->h_j)[row] - (params->h_i)[row]]);

      if(isnan((params->h_mij)[row]) || isinf((params->h_mij)[row])){
	(params->h_invalid)[row] = 1;
      }
    }else{
      (params->h_invalid)[row] = 1;
    }
  }

  return NULL;
}

/* normalization only */
void *hic_prep_thread_norm(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row, ij_dist;

  for(row = params-> begin; row <= params -> end; row++){
    ij_dist = (params->h_j)[row] - (params->h_i)[row];

    if(params->min_dist <= ij_dist && ij_dist <= params->max_dist){

      (params->h_mij)[row] /= ((params->norm)[(params->h_i)[row]] *
			       (params->norm)[(params->h_j)[row]]);			 

      if(isnan((params->h_mij)[row]) || isinf((params->h_mij)[row])){
	(params->h_invalid)[row] = 1;
      }
    }else{
      (params->h_invalid)[row] = 1;
    }
  }

  return NULL;
}

/* O/E conversion only */
void *hic_prep_thread_exp(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row, ij_dist;

  for(row = params-> begin; row <= params -> end; row++){
    ij_dist = (params->h_j)[row] - (params->h_i)[row];

    if(params->min_dist <= ij_dist && ij_dist <= params->max_dist){

      (params->h_mij)[row] /= (params->exp)[(params->h_j)[row] - (params->h_i)[row]];

      if(isnan((params->h_mij)[row]) || isinf((params->h_mij)[row])){
	(params->h_invalid)[row] = 1;
      }
    }else{
      (params->h_invalid)[row] = 1;
    }
  }
  return NULL;
}

/* O/E conversion only */
void *hic_prep_thread(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row, ij_dist;
  for(row = params-> begin; row <= params -> end; row++){
    ij_dist = (params->h_j)[row] - (params->h_i)[row];

    if(params->max_dist < ij_dist ||
       ij_dist < params->min_dist ||
       isnan((params->h_mij)[row]) ||
       isinf((params->h_mij)[row])){

      (params->h_invalid)[row] = 1;
    }
  }
  return NULL;
}


int hic_prep(const command_line_arguements *cmd_args,
	     hic **hic){
  struct timeval ts, tg;


  hic_raw *raw;
  hic_raw_read(cmd_args->hicRaw_dir,
	       cmd_args->res,  cmd_args->chr,	       
	       cmd_args->norm, cmd_args->exp,	       
	       &raw);

  gettimeofday(&ts, NULL);

  if(cmd_args->exec_thread_num >= 1){
    int i = 0;
    hic_prep_thread_args *params;
    pthread_t *threads = NULL;

    params = calloc_errchk(cmd_args->exec_thread_num,			   
			   sizeof(hic_prep_thread_args),
			   "calloc: hic_prep_thread_args");
    threads = calloc_errchk(cmd_args->exec_thread_num,			   
			    sizeof(pthread_t),
			    "calloc: threads");
        
    /* set variables */
    for(i = 0; i < cmd_args->exec_thread_num; i++){
      params[i].thread_id = i;
      params[i].begin = ((i == 0) ? 0 : params[i - 1].end + 1);
      params[i].end = ((i == (cmd_args->exec_thread_num - 1)) ?
		       raw->hic->nrow - 1 :
		       ((raw->hic->nrow / cmd_args->exec_thread_num) * (i + 1) - 1));
      params[i].h_invalid = raw->hic->invalid;
      params[i].h_i = raw->hic->i;
      params[i].h_j = raw->hic->j;
      params[i].h_mij = raw->hic->mij;
      params[i].min_dist = cmd_args->min_size / cmd_args->res;
      params[i].max_dist = cmd_args->max_size / cmd_args->res;
      params[i].norm = raw->norm;
      params[i].exp = raw->exp;
    }

    if(cmd_args->norm != NULL && cmd_args->exp != NULL){
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_create(&threads[i], NULL, 
		       hic_prep_thread_norm_exp,
		       (void*)&params[i]);
      }
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_join(threads[i], NULL);
      }
      free(raw->norm);
      free(raw->exp);
    }else if(cmd_args->norm != NULL){
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_create(&threads[i], NULL, 
		       hic_prep_thread_norm,
		       (void*)&params[i]);
      }
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_join(threads[i], NULL);
      }      
      free(raw->norm);
    }else if(cmd_args->exp != NULL){
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_create(&threads[i], NULL, 
		       hic_prep_thread_exp,
		       (void*)&params[i]);
      }
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_join(threads[i], NULL);
      }      
      free(raw->exp);
    }else{
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_create(&threads[i], NULL, 
		       hic_prep_thread,
		       (void*)&params[i]);
      }
      for(i = 0; i < cmd_args->exec_thread_num; i++){
	pthread_join(threads[i], NULL);
      }      
    }
    free(threads);
    free(params);
  }
  *hic = raw->hic;
  free(raw);

  gettimeofday(&tg, NULL);
  printf("%d\t%f\n", cmd_args->exec_thread_num, diffSec(ts, tg));

  return 0;
}


#if 0
int HicRead(const char *hic_file,
	    const unsigned int res,
	    int **h_i,
	    int **h_j,
	    double **h_mij,
	    unsigned long *nrow){
  FILE *fp;
  char buf[BUF_SIZE], tmp_mij_str[64];
  int tmp_i, tmp_j;
  unsigned long row = 0;

  *nrow = mywc(hic_file);

  /* calloc */
  {  
    *h_i = calloc_errchk(*nrow, sizeof(int), "calloc h_i");
    *h_j = calloc_errchk(*nrow, sizeof(int), "calloc h_j");
    *h_mij = calloc_errchk(*nrow, sizeof(double), "calloc h_mij");
  }

  if((fp = fopen(hic_file, "r")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    hic_file, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fp)){
    sscanf(buf, "%d\t%d\t%s", &tmp_i, &tmp_j, (char *)(&tmp_mij_str));
    (*h_i)[row] = tmp_i / res;
    (*h_j)[row] = tmp_j / res;
    (*h_mij)[row] = strtod(tmp_mij_str, NULL);
    row++;
  }
  
  fclose(fp);

  return 0;
}
#endif
#endif
