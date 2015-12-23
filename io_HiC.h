#ifndef __IO_HIC_H__
#define __IO_HIC_H__

#include <stdlib.h>
#include <stdio.h>
#include "mywc.h"

#define BUF_SIZE 4096

int HicRead(const char *file,
	    const int res,
	    int **h_i,
	    int **h_j,
	    double **h_mij,
	    unsigned long *nrow){
  FILE *fp;
  char buf[BUF_SIZE], tmp_mij_str[64];
  int tmp_i, tmp_j;
  unsigned long row = 0;

  *nrow = mywc(file);
  
  *h_i = calloc_errchk(*nrow, sizeof(int), "calloc h_i");
  *h_j = calloc_errchk(*nrow, sizeof(int), "calloc h_j");
  *h_mij = calloc_errchk(*nrow, sizeof(double), "calloc h_mij");

  if((fp = fopen(file, "r")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    file, strerror(errno));
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
