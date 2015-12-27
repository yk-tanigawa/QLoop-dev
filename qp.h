#ifndef __QP_H__
#define __QP_H__

#include <stdlib.h>
#include <stdio.h>
#include "constant.h"
#include "calloc_errchk.h"
#include "hic.h"
#include "kmer.h"
#include "adaboost.h"

int qp_show_P(FILE *fp, const unsigned int dim,
	      double **matrix){
  unsigned int i, j;
  for(i = 0; i < dim; i++){
    for(j = 0; j < dim; j++){
      fprintf(fp, "%e\t", matrix[i][j]);
    }
    fprintf(fp, "\n");
  }
  return 0;
}

int qp_show_q(FILE *fp, const unsigned int dim,
	      double *vector){
  unsigned int i;
  for(i = 0; i < dim; i++){
    fprintf(fp, "%e\n", vector[i]);
  }
  return 0;
}


int qp_prep(const command_line_arguements *cmd_args,
	    const unsigned int **kmer_freq,
	    hic *data,
	    const canonical_kp *kp,
	    adaboost *model,
	    double ***P,
	    double **q,
	    const char *qp_file_P,
	    const char *qp_file_q){
  unsigned long stamp, x, i, j;
  unsigned int *pair_freq;
  
  fprintf(stderr, 
	  "%s: info: QP: QP preparation with %ld variables(k-mer pair stamps)\n",
	  cmd_args->prog_name, model->T);

  /* allocate memory */
  {
    *P = calloc_errchk(model->T, sizeof(double *), "calloc: P");
    for(stamp = 0; stamp < model->T; stamp++){
      (*P)[stamp] = calloc_errchk(model->T, sizeof(double), "calloc: P[]");
    }
    *q = calloc_errchk(model->T, sizeof(double), "calloc: q");		       
    pair_freq = calloc_errchk(model->T,
			      sizeof(unsigned int), "calloc: pair_freq");
  }

  for(x = 0; x < data->nrow; x++){
    for(stamp = 0; stamp < model->T; stamp++){
      pair_freq[stamp] =
	kmer_freq[data->i[x]][kp->l1[model->axis[stamp]]] * 
	kmer_freq[data->j[x]][kp->m1[model->axis[stamp]]] +
	kmer_freq[data->i[x]][kp->l2[model->axis[stamp]]] * 
	kmer_freq[data->j[x]][kp->m2[model->axis[stamp]]];
    }
    for(i = 0; i < model->T; i++){
      for(j = 0; j < model->T; j++){
	(*P)[i][j] += pair_freq[i] * pair_freq[j];
      }
    }
    for(i = 0; i < model->T; i++){
      (*q)[i] -= data->mij[x] * pair_freq[i];
    }
  }

  {
    for(i = 0; i < model->T; i++){
      for(j = 0; j < model->T; j++){
	(*P)[i][j] /= ((data->nrow) * (data->nrow));
      }
    }
    for(i = 0; i < model->T; i++){
      (*q)[i] /= (data->nrow);
    }
  }


  /* write to file OR stderr */
  {
    if(qp_file_P == NULL || qp_file_q == NULL){
      qp_show_P(stderr, model->T, *P);
      qp_show_q(stderr, model->T, *q);
    }else{
      FILE *fp_P, *fp_q;
      {
	if((fp_P = fopen(qp_file_P, "w")) == NULL){
	  fprintf(stderr, "error: fopen %s\n%s\n",
		  qp_file_P, strerror(errno));
	  exit(EXIT_FAILURE);
	}
	fprintf(stderr, "%s: info: QP: writing matrix P to file: %s\n",
		cmd_args->prog_name, qp_file_P);
	qp_show_P(fp_P, model->T, *P);
	fclose(fp_P);
      }
      {
	if((fp_q = fopen(qp_file_q, "w")) == NULL){
	  fprintf(stderr, "error: fopen %s\n%s\n",
		  qp_file_q, strerror(errno));
	  exit(EXIT_FAILURE);
	}
	fprintf(stderr, "%s: info: QP: writing vector q to file: %s\n",
		cmd_args->prog_name, qp_file_q);
	qp_show_q(fp_q, model->T, *q);
	fclose(fp_q);
      }
    }
  }


  return 0;
}


#endif
