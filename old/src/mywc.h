#ifndef __MY_WC_H__
#define __MY_WC_H__

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "constant.h"

/**
 * wc -l
 */
unsigned long mywc(const char *file_name){
  unsigned long lines = 0;
  FILE* fp;
  char buf[MYWC_BUF_SIZE];
  
  if((fp = fopen(file_name, "r")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fp)){
    lines++;
  }

  fclose(fp);

  return lines;
}

/**
 * wc -b
 */
unsigned long mywc_b(const char *file_name){
  int fd;
  struct stat stbuf;

  if((fd = open(file_name, O_RDONLY)) == -1){
    fprintf(stderr, "error: open %s\n%s\n", 
	    file_name, strerror(errno));      
    exit(EXIT_FAILURE);
  }

  if(fstat(fd, &stbuf) == -1){
    fprintf(stderr, "error: fstat %s\n%s\n",
	    file_name, strerror(errno));
    exit(EXIT_FAILURE);
  }

  close(fd);

  return stbuf.st_size;
}


#endif
