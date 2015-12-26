#ifndef __KMER_H__
#define __KMER_H__

#include <stdlib.h>
#include <stdio.h>
#include "calloc_errchk.h"

/* compute reverse complement of a given k-mer */
unsigned long rev_comp(unsigned long kmer, 
		       const unsigned int k){
  unsigned int i = 0;
  unsigned long revComp = 0;
  kmer = ((~kmer) & ((1 << (2 * k)) - 1));
  for(i = 0; i <  k; i++){
    revComp <<= 2;
    revComp += (kmer & 3);
    kmer >>= 2;
  }
  return revComp;
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

int set_kmer_strings(const int k,
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
