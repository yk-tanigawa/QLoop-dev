#ifndef __QP_H__
#define __QP_H__

#include <stdlib.h>
#include <stdio.h>
#include "constant.h"
#include "calloc_errchk.h"
#include "hic.h"
#include "kmer.h"

int qp_prep(const command_line_arguements *cmd_args,
	    const unsigned int **kmer_freq,
	    hic *data,
	    const canonical_kp *kp,
	    double ***P,
	    double **q,
	    const char *qp_file){
  const unsigned long canonical_kmer_pair_num = 
    (1 << (4 * (cmd_args->k) - 1)) + (1 << (2 * (cmd_args->k) - 1));  
  unsigned long lm, x, i, j;
  unsigned int *pair_freq;
  
  /* allocate memory */
  {
    *P = calloc_errchk(canonical_kmer_pair_num,
		       sizeof(double *), "calloc: P");
    for(lm = 0; lm < canonical_kmer_pair_num; lm++){
      (*P)[lm] = calloc_errchk(canonical_kmer_pair_num,
			       sizeof(double), "calloc: P[]");
    }
    *q = calloc_errchk(canonical_kmer_pair_num,
		       sizeof(double), "calloc: q");
    pair_freq = calloc_errchk(canonical_kmer_pair_num,
			      sizeof(unsigned int), "calloc: pair_freq");

  }
#if 1
  for(x = 0; x < data->nrow; x++){
    for(lm = 0; lm < kp->num; lm++){
      pair_freq[lm] =
	kmer_freq[data->i[x]][kp->l1[lm]] * 
	kmer_freq[data->j[x]][kp->m1[lm]] +
	kmer_freq[data->i[x]][kp->l2[lm]] * 
	kmer_freq[data->j[x]][kp->m2[lm]];
    }
    for(i = 0; i < kp->num; i++){
      for(j = 0; j < kp->num; j++){
	(*P)[i][j] += pair_freq[i] * pair_freq[j];
      }
    }
    for(i = 0; i < kp->num; i++){
      (*q)[lm] -= data->mij[x] * pair_freq[i];
    }
  }
#endif

  return 0;
}


#endif
