#ifndef __KMER_H__
#define __KMER_H__

#include <stdlib.h>
#include <stdio.h>
#include "calloc_errchk.h"

typedef struct _canonical_kp{
  unsigned int *l1;
  unsigned int *m1;
  unsigned int *l2;
  unsigned int *m2;
  unsigned long kmer_num;
  unsigned long kmer_pair_num;
  unsigned long num;
} canonical_kp;

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

int set_canonical_kmer_pairs(const unsigned int k,
			     canonical_kp **kp){
  {
    *kp = calloc_errchk(1, sizeof(canonical_kp), "canonical_kp");
    (*kp)->kmer_num = 1 << (2 * k);
    (*kp)->kmer_pair_num = 1 << (4 * k);
    (*kp)->num = (1 << (4 * k - 1)) + (1 << (2 * k - 1));
    (*kp)->l1 = calloc_errchk((*kp)->num, sizeof(unsigned int), "canonical_kp l1");
    (*kp)->m1 = calloc_errchk((*kp)->num, sizeof(unsigned int), "canonical_kp m1");
    (*kp)->l2 = calloc_errchk((*kp)->num, sizeof(unsigned int), "canonical_kp l2");
    (*kp)->m2 = calloc_errchk((*kp)->num, sizeof(unsigned int), "canonical_kp m2");
  }

  {
    unsigned long lm = 0, revcomp_lm = 0, next = 0;
    for(lm = 0; lm < (*kp)->kmer_pair_num; lm++){
      revcomp_lm = rev_comp(lm, 2 * k);
      if(lm <= revcomp_lm){
	/**
	 * where, l2 = rev_comp(m1, k)
	 *        m2 = rev_comp(l1, k)
	 * note:
	 *        concat(l2 + m2) = rev_comp(concat(l1 + m1))
	 */	
	((*kp)->l1)[next] = lm / ((*kp)->kmer_num);
	((*kp)->m1)[next] = lm % ((*kp)->kmer_num);
	((*kp)->l2)[next] = revcomp_lm / ((*kp)->kmer_num);
	((*kp)->m2)[next] = revcomp_lm % ((*kp)->kmer_num);
	next++;
      }
    }
  }

  return 0;
}



#endif
