#ifndef __PRED_H__
#define __PRED_H__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "constant.h"
#include "cmd_args.h"

#include "calloc_errchk.h"
#include "hic.h"
#include "kmer.h"
#include "l2boost.h"

int predict(const cmd_args *,		
	    const double **,
	    const hic *,
	    const canonical_kp *,
	    const boost *,	 
	    double **,
	    FILE *);
int pred_cmp_file(const cmd_args *,
		  const hic *,
		  const double *,
		  FILE *);

int predict(const cmd_args *args,		
	    const double **feature,
	    const hic *data,
	    const canonical_kp *ckps,
	    const boost *model,	 
	    double **pred,
	    FILE *fp){
  const unsigned long n = data->nrow;
  const unsigned long p = ckps->num;
  const unsigned int *h_i = data->i;
  const unsigned int *h_j = data->j;
  const unsigned int *kmer1 = ckps->kmer1;
  const unsigned int *kmer2 = ckps->kmer2;
  const unsigned int *revcmp1 = ckps->revcmp1;
  const unsigned int *revcmp2 = ckps->revcmp2;
  unsigned long i, j;
  double pf;

  fprintf(fp, "%s [INFO] ", args->prog_name);
  fprintf(fp, "start prediction of interatcion intensities \n");

  /* allocate memory */
  *pred = calloc_errchk(n, sizeof(double),
			"calloc pred[]");

  /* pf : pairwise feature */
  for(j = 0; j < p; j++){
    if((model->beta[j]) != 0){
      for(i = 0; i < n; i++){
	pf =  ((feature[h_i[i]][kmer1[j]] *
		feature[h_j[i]][kmer2[j]]) +
	       (feature[h_i[i]][revcmp1[j]] *
		feature[h_j[i]][revcmp2[j]]));
	(*pred)[i] += (model->beta)[j] * pf;
      }
    }
  }

  return 0;
}

int pred_cmp_file(const cmd_args *args,
		  const hic *data,
		  const double *pred,
		  FILE *fp){
  const unsigned int *h_i = data->i;
  const unsigned int *h_j = data->j;
  const double *mij = data->mij;
  const unsigned long n = data->nrow;
  unsigned long i;
  FILE *fp_file;
  char out_file_name[BUF_SIZE];
  sprintf(out_file_name, 
	  "%s.cmp", args->out_file);

  if((fp_file = fopen(out_file_name, "w")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    out_file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "%s [INFO] ", args->prog_name);
  fprintf(fp, "start writing cmp file to : %s\n",
	  out_file_name);

  fprintf(fp_file, "i\tj\tobs\tpred\n");

  for(i = 0; i < n; i++){
    fprintf(fp_file, "%d\t%d\t%e\t%e\n",
	    h_i[i], h_j[i], mij[i], pred[i]);
  }

  fclose(fp_file);


  return 0;
}

#endif

