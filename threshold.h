#ifndef __THRESHOLD_H__
#define __THRESHOLD_H__

#include <stdlib.h>
#include <stdio.h>
#include "constant.h"
#include "calloc_errchk.h"

typedef struct _thresholds {
  unsigned int nclass;
  double *representatives;
} thresholds;


int dump_thresholds(FILE *fp,
		    thresholds *th){
  unsigned int i;
  for(i = 0; i <= th->nclass; i++){
    fprintf(fp, "%d\t%f\n", i, th->representatives[i]);
  }
  return 0;
}

int write_histo(const command_line_arguements *cmd_args,
		thresholds *th,
		const char *out_file){
  FILE *fp;
  if((fp = fopen(out_file, "w")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    out_file, strerror(errno));
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "%s: info: histo: writing histogram to file: %s\n",
	  cmd_args->prog_name, out_file);
  dump_thresholds(fp, th);
  fclose(fp);
  return 0;
}

int double_comp(const void *cmp1,
		const void *cmp2){
  double d1 = *((double *)cmp1);
  double d2 = *((double *)cmp2);
  if(d1 < d2){
    return -1;
  }else if(d1 > d2){
    return 1;
  }else{
    return 0;
  } 
}

int set_thresholds(const double *mij,
		   const unsigned int nclass,
		   const unsigned long num,
		   thresholds **t){
  double *cpy;
  /* allocate memory */
  {
    *t = calloc_errchk(1, sizeof(thresholds), "calloc thresholds");
    (*t)->representatives = calloc_errchk(nclass + 1, sizeof(double),
					  "calloc: thresholds->representatives");
    (*t)->nclass = nclass;
  }
  /* copy and sort */
  {
    unsigned long i;
    cpy = calloc_errchk(num, sizeof(double), "calloc: cpy");
    for(i = 0; i < num; i++){
      cpy[i] = mij[i];
    }   
    qsort(cpy, num, sizeof(double), double_comp);
  }
  {
    unsigned int class;
    unsigned long index;
    for(class = 0; class < nclass; class++){
      index = (int)(1.0 * num * class / nclass + 0.5);
      (*t)->representatives[class] = cpy[index];
    }
    (*t)->representatives[nclass] = cpy[num - 1];
  }
  free(cpy);
  return 0;
}

double get_threshold(const command_line_arguements *cmd_args,
		     thresholds *th,
		     double percentile){
  int index = (int)(th->nclass * percentile + 0.5);

  fprintf(stderr, "%s: info: threshold: %e percentile = %f\n",
	  cmd_args->prog_name, 
	  percentile, (th->representatives)[index]);

  return (th->representatives)[index];  
}

#if 0
double get_threshold(const double *mij,
		     const double percentile,
		     const unsigned long num){
  unsigned long i;
  double threshold = 0;
  int index = (int)(num * percentile), int_percentile = (int)(100 * percentile);
  double *cpy = calloc_errchk(num, sizeof(double), "calloc: cpy");
  for(i = 0; i < num; i++){
    cpy[i] = mij[i];
  }
  
  qsort(cpy, num, sizeof(double), double_comp);

  threshold = cpy[index];

  fprintf(stderr, "info: %d %% percentile = %f\n",
	  int_percentile, threshold);

  free(cpy);
  return threshold;
}
#endif
#endif
