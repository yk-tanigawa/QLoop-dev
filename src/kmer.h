#ifndef __KMER_H__
#define __KMER_H__

#include <stdlib.h>
#include <stdio.h>
#include "calloc_errchk.h"
#include "mywc.h"

typedef struct _canonical_kp{
  unsigned int *kmer1;
  unsigned int *kmer2;
  unsigned int *revcmp1; 
  unsigned int *revcmp2;
  unsigned long num;
} canonical_kp;

typedef struct _kmer{
  unsigned int *kmer1;
  int k;
  unsigned long num;
} kmer;

int canonical_kp_read(const cmd_args *, canonical_kp **);

int kmer_read(const cmd_args *, kmer **);

int canonical_kp_read(const cmd_args *args,
		      canonical_kp **ckps){
  {
    /* allocate memory */
    *ckps        = calloc_errchk(1, sizeof(canonical_kp), "calloc ckps");   
    (*ckps)->num = mywc(args->kmer_pair);
    (*ckps)->kmer1    = calloc_errchk((*ckps)->num, sizeof(unsigned int),
				      "calloc ckps (*ckps)->kmer1");
    (*ckps)->kmer2    = calloc_errchk((*ckps)->num, sizeof(unsigned int),
				      "calloc ckps (*ckps)->kmer2");
    (*ckps)->revcmp1  = calloc_errchk((*ckps)->num, sizeof(unsigned int),
				      "calloc ckps (*ckps)->revcmp1");
    (*ckps)->revcmp2  = calloc_errchk((*ckps)->num, sizeof(unsigned int),
				      "calloc ckps (*ckps)->revcmp2");
  }

  /* read from a file */
  {
    FILE *fp;
    char buf[BUF_SIZE];
    char kmer1str[BUF_SIZE], kmer2str[BUF_SIZE];
    char revcmp1str[BUF_SIZE], revcmp2str[BUF_SIZE];
    unsigned int kmer1, kmer2, revcmp1, revcmp2;
    unsigned long row = 0;

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start reading canonical k-mer pair file from %s\n",
	    args->kmer_pair);

    if((fp = fopen(args->kmer_pair, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      args->kmer_pair, strerror(errno));
      exit(EXIT_FAILURE);
    }
    
    while(fgets(buf, BUF_SIZE, fp) && row < (*ckps)->num){
      sscanf(buf, "%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s",
	     &kmer1, &kmer2, &revcmp1, &revcmp2, 
	     kmer1str, kmer2str, revcmp1str, revcmp2str);
      ((*ckps)->kmer1)[row]   = kmer1;
      ((*ckps)->kmer2)[row]   = kmer2;
      ((*ckps)->revcmp1)[row] = revcmp1;
      ((*ckps)->revcmp2)[row] = revcmp2;
      row++;
    }
  
    fclose(fp);

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "# of canonical k-mer pairs = %ld\n",
	    (*ckps)->num);
  }

  return 0;
}

int kmer_read(const cmd_args *args,
	      kmer **kmers){
  {
    /* allocate memory */
    *kmers        = calloc_errchk(1, sizeof(kmer), "calloc kmers");   
    (*kmers)->num = mywc(args->kmer);
    (*kmers)->kmer1 = calloc_errchk((*kmers)->num, sizeof(unsigned int),
				    "calloc kmers (*kmers)->kmer1");
  }

  /* read from a file */
  {
    FILE *fp;
    char buf[BUF_SIZE];
    char kmer1str[BUF_SIZE];
    unsigned int kmer1;
    unsigned long row = 0;

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start reading k-mer file from %s\n",
	    args->kmer);

    if((fp = fopen(args->kmer, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      args->kmer, strerror(errno));
      exit(EXIT_FAILURE);
    }
    
    while(fgets(buf, BUF_SIZE, fp) && row < (*kmers)->num){
      sscanf(buf, "%d\t%s",
	     &kmer1, kmer1str);
      ((*kmers)->kmer1)[row]   = kmer1;
      row++;
    }
  
    fclose(fp);

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "# of k-mers = %ld\n",
	    (*kmers)->num);
  }

  return 0;
}


#endif
