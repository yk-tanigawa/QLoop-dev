#ifndef __io_H__
#define __io_H__

#include <stdlib.h>
#include <stdio.h>

#define BUF_SIZE 4096


long wc(const char *fName){
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

int readTableInt(const char *file,
		 const char *dlim,
		 const unsigned long ncol,
		 int ***table,
		 unsigned long *nrow){
  FILE *fp;
  char buf[BUF_SIZE];
  char *tok;
  unsigned long col = 0, row = 0;

  *nrow = wc(file);
  
  *table = calloc_errchk(sizeof(int *), *nrow, "calloc table");

  if((fp = fopen(file, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
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


inline char Binary2char(const long binaryNum){
  if(binaryNum == 0){
    return 'A'; 
  }else if(binaryNum == 1){
    return 'C';
  }else if(binaryNum == 2){
    return 'G'; 
  }else if(binaryNum == 3){
    return 'T';
  }else{
    return '?'; 
  }
}

int setKmerStrings(const int k,
		   char ***kmerStrings){
  const unsigned long kmerNum = 1 << (2 * k);
  int i;
  unsigned long l, m;
  *kmerStrings = calloc_errchk(sizeof(int *), kmerNum, "calloc kmerStrings");
  for(l = 0; l < kmerNum; l++){
    (*kmerStrings)[l] = calloc_errchk(sizeof(char), k + 1, "calloc kmerStrings[l]");
    m = l;
    for(i = k - 1; i >= 0; i--){
      (*kmerStrings)[l][i] = Binary2char((m & 3));
      m = (m >> 2);
    }
  }
  return 0;
}

#endif
