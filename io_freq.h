#ifndef __IO_FREQ_H_
#define __IO_FREQ_H_

#include <stdlib.h>
#include <stdio.h>
#include "mywc.h"

int kmerFreqRead(const char *file,
		 const char *dlim,
		 const unsigned long ncol,
		 int ***table,
		 unsigned long *nrow){
  FILE *fp;
  char buf[BUF_SIZE];
  char *tok;
  unsigned long col = 0, row = 0;

  *nrow = mywc(file);
  
  *table = calloc_errchk(sizeof(int *), *nrow, "calloc table");

  if((fp = fopen(file, "r")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    file, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fp)){
    if(buf[0] == '*'){
      (*table)[row] = NULL;
    }else{
      (*table)[row] = calloc_errchk(sizeof(int), ncol, "calloc table[l]");
      col = 0;
      tok = strtok(buf, dlim);
      while(tok != NULL && col < ncol){
	(*table)[row][col++] = atoi(tok);
	tok = strtok(NULL, dlim);
      }
    }
    row++;
  }
  
  fclose(fp);

  return 0;
}

#endif
