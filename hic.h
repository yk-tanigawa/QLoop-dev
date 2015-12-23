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


/* normalized O/E converted Hi-C data */
typedef struct _hic {
  unsigned long nrow;
  unsigned int res;
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
  int *h_i;
  int *h_j;
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
    (*data)->i = calloc_errchk((*data)->nrow, sizeof(unsigned int), "calloc hic (*data)->i");
    (*data)->j = calloc_errchk((*data)->nrow, sizeof(unsigned int), "calloc hic (*data)->j");
    (*data)->mij = calloc_errchk((*data)->nrow, sizeof(double), "calloc hic (*data)->mij");
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
      ((*data)->i)[row] = tmp_i / res;
      ((*data)->j)[row] = tmp_j / res;
      ((*data)->mij)[row] = strtod(tmp_mij_str, NULL);
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
int set_hic_file_names(const char *hicDir,
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


  return 0;
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
 * normalization and or O/E conversion
 */

/* normalization and O/E conversion */
void *hic_prep_thread_norm_exp(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row;

  for(row = params-> begin; row <= params -> end; row++){
    (params->h_mij)[row] /= ((params->norm)[(params->h_i)[row]] *
			     (params->norm)[(params->h_j)[row]] *
			     (params->exp)[(params->h_j)[row] - (params->h_i)[row]]);
  }
  return NULL;
}

/* normalization only */
void *hic_prep_thread_norm(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row;

  for(row = params-> begin; row <= params -> end; row++){
    (params->h_mij)[row] /= ((params->norm)[(params->h_i)[row]] *
			     (params->norm)[(params->h_j)[row]]);			 
  }
  return NULL;
}

/* O/E conversion only */
void *hic_prep_thread_exp(void *args){
  hic_prep_thread_args *params = (hic_prep_thread_args *)args;
  unsigned long row;

  for(row = params-> begin; row <= params -> end; row++){
    (params->h_mij)[row] /= (params->exp)[(params->h_j)[row] - (params->h_i)[row]];
  }
  return NULL;
}

int hic_prep(const command_line_arguements *cmd_args,
	     hic **hic){
  hic_raw *raw;
  hic_raw_read(cmd_args->hicRaw_dir,
	       cmd_args->res,  cmd_args->chr,	       
	       cmd_args->norm, cmd_args->exp,	       
	       &raw);
  *hic = raw->hic;
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
