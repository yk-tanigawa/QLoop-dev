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

#include "calloc_errchk.h"
#include "io.h"

#define F_NAME_LEN 128
#define BUF_SIZE 4096

int readVectorInt(const char *file,
		 int **vector,
		 unsigned long *nrow){
  FILE *fp;
  char buf[BUF_SIZE];
  unsigned long row = 0;

  *nrow = wc(file);
  
  *vector = calloc_errchk(sizeof(int), *nrow, "calloc vector");

  {
    if((fp = fopen(file, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      file, strerror(errno));
      exit(EXIT_FAILURE);
    }
    
    while(fgets(buf, BUF_SIZE, fp)){
      (*vector)[row++] = atoi(buf);
    }
  
    fclose(fp);
  }

  return 0;
}

int set_lm(const int *lm,
	   const unsigned long num,
	   const int k,
	   int **l,
	   int **m){
  const int nkmer = 1 << (2 * k);
  unsigned long i;
  
  *l = calloc_errchk(num, sizeof(int), "calloc l");
  *m = calloc_errchk(num, sizeof(int), "calloc m");
  
  for(i = 0; i < num; i++){
    (*l)[i] = lm[i] % nkmer;
    (*m)[i] = lm[i] / nkmer;
  }
  return 0;
}

int dumpResultsQP_q(FILE *fp,
		    const double *q, 
		    const unsigned long nkp){
  unsigned int i;
  for(i = 0; i < nkp; i++){
    fprintf(fp, "%e\n", q[i]);
  }
  return 0;
}

int dumpResultsQP_P(FILE *fp,
		    const double **P, 
		    const unsigned long nkp){
  unsigned int i, j;
  for(i = 0; i < nkp; i++){
    for(j = 0; j < nkp; j++){
      fprintf(fp, "%e\t", P[i][j]);
    }
    fprintf(fp, "\n");
  }
  return 0;
}


int main_sub(const char *freqFile,
	     const char *hicFile,
	     const char *kmerPairFile,
	     const int k,
	     const int res,
	     const char *outDir){

  int **kmerFreq;
  unsigned long nBin;

  int *h_i, *h_j;
  double *h_mij;
  unsigned long nHic;
  
  int *kp_lm, *kp_l, *kp_m;
  unsigned long nkp;

  double **P, *q;
  int *pairFreq;

  unsigned long b, n, lm, l, m;

  char outFileP[F_NAME_LEN], outFileq[F_NAME_LEN];

  FILE *fpP, *fpq;


  /* generate data on the fly mode */
  /* load data */
  {
    readTableInt(freqFile, "\t", 1 << (2 * k), &kmerFreq, &nBin);
    readHic(hicFile, res, &h_i, &h_j, &h_mij, &nHic);
    readVectorInt(kmerPairFile, &kp_lm, &nkp);
    set_lm(kp_lm, nkp, k, &kp_l, &kp_m);
  }
  
  /* execute summation */  
  {
    P = calloc_errchk(nkp, sizeof(double *), "calloc P");
    for(b = 0; b < nkp; b++){
      P[b] = calloc_errchk(nkp, sizeof(double), "calloc P[]");
    }    
    q = calloc_errchk(nkp, sizeof(double), "calloc q");
    pairFreq = calloc_errchk(nkp, sizeof(int), "calloc pairFreq");

    for(n = 0; n < nHic; n++){
      for(lm = 0; lm < nkp; lm++){
	pairFreq[lm] = kmerFreq[h_i[n]][kp_l[lm]] * kmerFreq[h_j[n]][kp_m[lm]];
      }
      for(l = 0; l < nkp; l++){
	for(m = 0; m < nkp; m++){
	  P[l][m] += pairFreq[l] * pairFreq[m];
	}
      }
      for(lm = 0; lm < nkp; lm++){
	q[lm] -= h_mij[n] * pairFreq[lm];
      }
    }

  }

  /* write or show the results */
  {
    if(outDir == NULL){
      dumpResultsQP_P(stderr, (const double **)P, nkp);
      dumpResultsQP_q(stderr, (const double *)q, nkp);
    }else{
      sprintf(outFileP, "%s%s.P%ld.dat",
	      outDir, basename((char *)kmerPairFile), nkp);
      sprintf(outFileq, "%s%s.q%ld.dat",
	      outDir, basename((char *)kmerPairFile), nkp);

      if((fpP = fopen(outFileP, "w")) == NULL){
	fprintf(stderr, "error: fdopen %s\n%s\n",
		outFileP, strerror(errno));
	exit(EXIT_FAILURE);
      }
      
      fprintf(stderr, "Writing results to : %s\n", outFileP);
      
      dumpResultsQP_P(fpP, (const double **)P, nkp);
      
      fclose(fpP);

      if((fpq = fopen(outFileq, "w")) == NULL){
	fprintf(stderr, "error: fdopen %s\n%s\n",
		outFileq, strerror(errno));
	exit(EXIT_FAILURE);
      }
      
      fprintf(stderr, "Writing results to : %s\n", outFileq);
      
      dumpResultsQP_q(fpq, (const double *)q, nkp);
      
      fclose(fpq);

    }
  }

  /* free memory */
  {
    for(b = 0; b < nBin; b++){
      if(kmerFreq[b] != NULL){
	free(kmerFreq[b]);
      }
    }
    free(kmerFreq);
    free(h_i);
    free(h_j);
    free(h_mij);
    free(kp_lm);
    free(kp_l);
    free(kp_m);
    for(b = 0; b < nkp; b++){
      free(P[b]);
    }
    free(P);
    free(q);
    free(pairFreq);
  }

  return 0;
}

int show_usage(const char *progName){
  printf("usage:\n\n");
  printf("example: %s -f %s -H %s -p %s -k %d -r %d -o %s\n", 
	 progName,
	 "./data/hoge.dat",
	 "./data/hic.dat",
	 "./data/kmerpairfile",
	 5,
	 1000,
	 "./data/");
  return 0;
}

int check_params(const char *freqFile,
		 const char *hicFile,
		 const char *kmerPairFile,
		 const int k,
		 const int res,
		 const char *outDir,
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
    printf("  Hi-C data: %s\n", hicFile);
  }

  if(kmerPairFile == NULL){
    fprintf(stderr, "[ERROR] kmerPairFile is not specified\n");
    errflag++;
  }else if(errflag == 0){
    printf("  k-mer pair file: %s\n", kmerPairFile);
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

  if(outDir == NULL){
    fprintf(stderr, "[WARNING] output directory is not specified\n");
    fprintf(stderr, "          results will be written to stderr\n");
  }else if(errflag == 0){
    printf("  Output dir:    %s\n", outDir);
  }

  if(errflag > 0){
    show_usage(progName);
    exit(EXIT_FAILURE);
  }

  return 0;
}

int main(int argc, char **argv){
  char *freqFile = NULL, *hicFile = NULL, *kmerPairFile = NULL, *outDir = NULL;
  int k = 0, res = 0;

  int opt = 0, opt_idx = 0;
  struct option long_opts[] = {
    {"help",    no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    /* options */
    {"frequency", required_argument, NULL, 'f'},
    {"hic",       required_argument, NULL, 'H'},
    {"kmerpair",  required_argument, NULL, 'p'},
    {"k",         required_argument, NULL, 'k'},
    {"res",       required_argument, NULL, 'r'},
    {"out",       required_argument, NULL, 'o'},
    {0, 0, 0, 0}
  };

  while((opt = getopt_long(argc, argv, "hvf:H:p:k:r:o:",
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
      case 'p': /* kmerpair */
	kmerPairFile = optarg;
	break;
      case 'k': /* k */
	k = atoi(optarg);
	break;
      case 'r': /* res */
	res = atoi(optarg);
	break;
      case 'o': /* out */
	outDir = optarg;
	break;
    }
  }

  check_params(freqFile, hicFile, kmerPairFile, k, res, outDir, argv[0]);

#if 0
  if(outDir != NULL){
    outFile = calloc_errchk(F_NAME_LEN, sizeof(char), "calloc outFile");
    if(bitmode == 0){
      sprintf(outFile, "%s%s.k%d.t%f.T%ld.stamps",
	      outDir, basename(hicFile), k, threshold, T);
    }else{
      sprintf(outFile, "%s%s.k%d.t%f.T%ld.bit.stamps",
	      outDir, basename(hicFile), k, threshold, T);
    }
  }
#endif
  main_sub(freqFile, hicFile, kmerPairFile, k, res, outDir);
 
  return 0;
}
