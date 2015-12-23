#ifndef __MY_WC_H__
#define __MY_WC_H__

#include <stdlib.h>
#include <stdio.h>

#define BUF_SIZE 4096

/**
 * wc -l
 */
long mywc(const char *fName){
  long lines = 0;
  FILE* fp;
  char buf[BUF_SIZE];
  
  if((fp = fopen(fName, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    fName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fp)){
    lines++;
  }

  fclose(fp);

  return lines;
}


#endif
