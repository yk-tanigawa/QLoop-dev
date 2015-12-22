#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <libgen.h>

#include "adaboost.h"
#include "calloc_errchk.h"
#include "io.h"

#define F_NAME_LEN 128
#define BUF_SIZE 4096

#if 0
int binarization(const double *source,
		 const double threshold,
		 const unsigned long length,
		 unsigned int **target){
  unsigned long i;
  *target = calloc_errchk(length, sizeof(int), "calloc target");
  for(i = 0; i < length; i++){
    (*target)[i] = ((source[i] > threshold) ? 1 : 0);
  }
  return 0;
}
#endif

int binarization(const double *source,
		 const double threshold,
		 const unsigned long length,
		 unsigned int **target){
  const size_t intBits = 8 * sizeof(unsigned int);
  unsigned long i;
  unsigned int buf = 0;
  *target = calloc_errchk((length + intBits - 1) / intBits, 
			  sizeof(unsigned int), "calloc target");
  {
    for(i = 0; i < length; i++){
      buf = (buf << 1);
      buf += ((source[i] > threshold) ? 1 : 0);
      if(((i + 1) % intBits) == 0){
	(*target)[i / intBits] = buf;
	buf = 0;
      }
    }
    if((length % intBits) != 0){
      buf = (buf << (intBits - (length % intBits)));
      (*target)[length / intBits] = buf;
    }
  }
  return 0;
}


int main_sub(const char *freqFile,
	     const char *hicFile,
	     const int k,
	     const int res,
	     const double threshold,
	     const unsigned long T,
	     const char *outFile,
	     const int bitmode){

  int **kmerFreq;
  unsigned long nBin, b;

  int *h_i, *h_j;
  double *h_mij;
  unsigned int **x;
  unsigned int *y;
  unsigned long nHic;

  unsigned long *adaAxis;
  int *adaSign;
  double *adaBeta;

  FILE *fp;

  if(bitmode == 0){
    /* generate data on the fly mode */
    /* load data */
    {
      readTableInt(freqFile, "\t", 1 << (2 * k), &kmerFreq, &nBin);
      readHic(hicFile, res, &h_i, &h_j, &h_mij, &nHic);
      binarization(h_mij, threshold, nHic, &y);
      free(h_mij);
    }

    /* execute adaboost */
    {
      adaboostLearn((const unsigned int *)y, 
		    (const int *)h_i, (const int*)h_j, 
		    (const int **)kmerFreq, 
		    k, T, (const unsigned long)nHic, 
		    1 << (4 * k),
		    &adaAxis, &adaSign, &adaBeta);
    }
  }else{
    /* bitmode */
    /* load data */
    {
      readTableInt(freqFile, "\t", 1 << (2 * k), &kmerFreq, &nBin);
      readHic(hicFile, res, &h_i, &h_j, &h_mij, &nHic);

      binarization(h_mij, threshold, nHic, &y);
      constructBitTable((const int *)h_i, (const int*)h_j, 			
			(const int **)kmerFreq, 
			(const unsigned long)nHic, k,
			&x);

      for(b = 0; b < nBin; b++){
	free(kmerFreq[b]);
      }
      free(kmerFreq);
      free(h_i);
      free(h_j);
      free(h_mij);
    }

    /* execute adaboost */
    {
      adaboostBitLearn((const unsigned int *)y, 
		       (const unsigned int **)x,
		       T, (const unsigned long)nHic, 
		       1 << (4 * k),
		       &adaAxis, &adaSign, &adaBeta);
    }
  }

  /* write or show the results */
  {
    if(outFile == NULL){
      dump_results(stderr, adaAxis, adaSign, adaBeta, T, k);
    }else{
      if((fp = fopen(outFile, "w")) == NULL){
	fprintf(stderr, "error: fdopen %s\n%s\n",
		outFile, strerror(errno));
	exit(EXIT_FAILURE);
      }
      
      fprintf(stderr, "Writing results to : %s\n", outFile);
      
      dump_results(fp, adaAxis, adaSign, adaBeta, T, k);
      
      fclose(fp);
    }
  }

  /* free memory */
  {
    if(bitmode == 0){
      for(b = 0; b < nBin; b++){
	free(kmerFreq[b]);
      }
      free(kmerFreq);
      free(h_i);
      free(h_j);
    }else{
      //for(b = 0; b < nHic; b++){
      for(b = 0; b < (unsigned long)(1 << (4 * k)); b++){
	free(x[b]);
      }
      free(x);
    }
    free(y);
    free(adaAxis);
    free(adaSign);
    free(adaBeta);
  }
  return 0;
}

int show_usage(const char *progName){
  printf("usage:\n\n");
  printf("example: %s -f %s -H %s -k %d -r %d, -t %f -T %ld -o %s\n", 
	 progName,
	 "./data/hoge.dat",
	 "./data/hic.dat",
	 5,
	 1000,
	 10.0,
	 10L,
	 "./data/");
  return 0;
}

int check_params(const char *freqFile,
		 const char *hicFile,
		 const int k,
		 const int res,
		 const double threshold,
		 const unsigned long T,
		 const char *outDir,
		 const int bitmode,
		 const char *progName){
  int errflag = 0;

  if(freqFile == NULL){
    fprintf(stderr, "[ERROR] input frequency file is not specified\n");
    errflag++;
  }else if(errflag == 0){
    printf("  k-mer frequency file:    %s\n", freqFile);
  }

  if(hicFile == NULL){
    fprintf(stderr, "[ERROR] Hi-C data directory is not specified\n");
    errflag++;
  }else if(errflag == 0){
    printf("  Hi-C data dir: %s\n", hicFile);
  }

  if(k <= 0){
    fprintf(stderr, "[ERROR] k is not specified\n");
    errflag++;
  }else if(errflag == 0){
    printf("  k: %d\n", k);
  }

  if(res <= 0){
    fprintf(stderr, "[ERROR] resolution is not specified\n");
    errflag++;
  }else if(errflag == 0){
    printf("  resolution: %d\n", res);
  }

  if(threshold <= 0){
    fprintf(stderr, "[ERROR] threshold is not specified\n");
    errflag++;
  }else if(errflag == 0){
    printf("  threshold: %e\n", threshold);
  }

  if(T <= 0){
    fprintf(stderr, "[ERROR] number of weak lerners is not specified\n");
    errflag++;
  }else if((unsigned long)(1 << (4 * k)) < T){
    fprintf(stderr, "[ERROR] number of weak lerners(T = %ld) exceeds 16^k (%d)\n", T, 1 << (4 * k));
    errflag++;
  }else if(errflag == 0){
    printf("  num. of weak lerners: %ld\n", T);
  }

  if(outDir == NULL){
    fprintf(stderr, "[WARNING] output directory is not specified\n");
    fprintf(stderr, "          results will be written to stderr\n");
  }else if(errflag == 0){
    printf("  Output dir:    %s\n", outDir);
  }

  if(bitmode != 0 && errflag == 0){
    printf("  [INFO] bitmode\n");
  }

  if(errflag > 0){
    show_usage(progName);
    exit(EXIT_FAILURE);
  }

  return 0;
}


#if 0
int main(int argc, char **argv){
  char *freqFile = NULL, *hicFile = NULL, *outDir = NULL, *outFile = NULL;
  int k = 0, res = 0, bitmode = 0;
  unsigned long T = 0;
  double threshold = 0;

  int opt = 0, opt_idx = 0;
  struct option long_opts[] = {
    {"help",    no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    /* options */
    {"frequency", required_argument, NULL, 'f'},
    {"hic",       required_argument, NULL, 'H'},
    {"k",         required_argument, NULL, 'k'},
    {"res",       required_argument, NULL, 'r'},
    {"threshold", required_argument, NULL, 't'},
    {"T",         required_argument, NULL, 'T'},
    {"out",       required_argument, NULL, 'o'},
    {"bitmode",   no_argument,       NULL, 'b'},
    {0, 0, 0, 0}
  };

  while((opt = getopt_long(argc, argv, "hvf:H:k:r:t:T:o:b",
			   long_opts, &opt_idx)) != -1){
    switch (opt){
      case 'h': /* help */
	show_usage(argv[0]);
	exit(EXIT_SUCCESS);
      case 'v': /* version*/
	printf("version: 0.10\n");
	exit(EXIT_SUCCESS);
      case 'f': /* frequency */
	freqFile = optarg;
	break;
      case 'H': /* hic */
	hicFile = optarg;
	break;
      case 'k': /* k */
	k = atoi(optarg);
	break;
      case 'r': /* res */
	res = atoi(optarg);
	break;
      case 't': /* threshold */
	threshold = atof(optarg);
	break;
      case 'T': /* T */
	T = atol(optarg);
	break;
      case 'o': /* out */
	outDir = optarg;
	break;
      case 'b': /* bitmode */
	bitmode = 1;
	break;
    }
  }

  check_params(freqFile, hicFile, k, res, threshold, T, outDir, bitmode, argv[0]);

  if(outDir != NULL){
    outFile = calloc_errchk(F_NAME_LEN, sizeof(char), "calloc outFile");
    if(bitmode == 0){
      sprintf(outFile, "%s%s.k%d.t%f.T%ld.stamps",
	      outDir, basename(hicFile), k, threshold, T);
    }else{
      sprintf(outFile, "%s%s.k%d.t%f.T%ld.bit5.stamps",
	      outDir, basename(hicFile), k, threshold, T);
    }
  }

  main_sub(freqFile, hicFile, k, res, threshold, T, outFile, bitmode);
 
  return 0;
}
#endif


long revComp(long kmer, 
	     const int k){
  long revComp = 0;
  kmer = ((~kmer) & ((1 << (2 * k)) - 1));
  int i;
  for(i = 0; i <  k; i++){
    revComp <<= 2;
    revComp += (kmer & 3);
    kmer >>= 2;
  }
  return revComp;
}


int main(void){
  int k = 3;
  long i = 0;
  char **kmerStrings;

  setKmerStrings(k, &kmerStrings);
  
  for(i = 0; i < (1 << (2 * k)); i++){
    printf("%s\t", kmerStrings[i]);
    printf("%s\n", kmerStrings[revComp(i, k)]);
  }

  return 0;
}
