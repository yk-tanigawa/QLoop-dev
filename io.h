#ifndef __io_H__
#define __io_H__

#include "hic.h"
#include "mywc.h"
#include "calloc_errchk.h"

int read_double(const char *fileName,
		double **array,
		unsigned long *len){
  {
    *len = mywc(fileName);
    *array = calloc_errchk(*len, sizeof(double), "error(calloc) readDouble");
  }

  {
    FILE *fp;
    unsigned long i = 0;
    char buf[BUF_SIZE];

    if((fp = fopen(fileName, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      fileName, strerror(errno));
      exit(EXIT_FAILURE);
    }

    while(fgets(buf, BUF_SIZE, fp) != NULL && i < *len){
      (*array)[i++] = strtod(buf, NULL);
    }

    fclose(fp);
  }

  return 0;
}


#endif
