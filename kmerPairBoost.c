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

#include "calloc_errchk.h"
#include "io.h"

#define F_NAME_LEN 128
#define BUF_SIZE 4096

inline void adaDataOnTheFly(const int h_i,
			 const int h_j,
			 const int **kmerFreq,
			 const int k,
			 short *x){
  const unsigned long nkmers = 1 << (2 * k);
  unsigned long l, m;
  for(l = 0; l < nkmers; l++){
    for(m = 0; m < nkmers; m++){
      x[l * nkmers + m] = (((kmerFreq[h_i][l] * kmerFreq[h_j][m]) > 0) ? 1 : 0);
    }
  }   
  return;
}

inline short adaDataOnTheFlyAxis(const int h_i,
			      const int h_j,
			      const int **kmerFreq,
			      const int k,
			      const int axis){
  const unsigned long nkmers = 1 << (2 * k);
  return (((kmerFreq[h_i][axis / nkmers] * kmerFreq[h_j][axis % nkmers]) > 0) ? 1 : 0);
}


int adaboostLearn(const int *y,
		  const int *h_i,
		  const int *h_j,
		  const int **kmerFreq,
		  const int k,
		  const unsigned long T,
		  const unsigned long N,
		  const unsigned long dim,
		  unsigned long **adaAxis,
		  int **adaSign,
		  double **adaBeta){
  short *x, xaxis;
  double *w, *p, *err, wsum, epsilon, min, max;
  unsigned long t, i, d, argmind, argmaxd;

  *adaAxis = calloc_errchk(T, sizeof(unsigned long), "calloc adaAxis");
  *adaSign = calloc_errchk(T, sizeof(int), "calloc adaSign");
  *adaBeta = calloc_errchk(T, sizeof(double), "calloc adaBeta");

  x = calloc_errchk(dim, sizeof(short), "calloc x");
  w = calloc_errchk(N, sizeof(double), "calloc w");
  p = calloc_errchk(N, sizeof(double), "calloc p");
  err = calloc_errchk(dim, sizeof(double), "calloc count");

  for(i = 0; i < N; i++){
    w[i] = 1.0 / N;
  }

  for(t = 0; t < T; t++){
    /* step 1 : compute normalized weights p[] */
    {
      wsum = 0;
      for(i = 0; i < N; i++){
	wsum += w[i];
      }
      for(i = 0; i < N; i++){
	p[i] = w[i] / wsum;
      }
    }

    /* step 2 : find the most appropriate axis (weak lerner) */
    {
      for(d = 0; d < dim; d++){
	err[d] = 0;
      }
      for(i = 0; i < N; i++){
	/* compute data point x */
	adaDataOnTheFly(h_i[i], h_j[i], kmerFreq, k, x);
	for(d = 0; d < dim; d++){
	  if(x[d]!= y[i]){
	    err[d] += p[i];
	  }
	}
      }

      {
	max = min = err[0];
	argmaxd = argmind = 0;
	for(d = 1; d < dim; d++){
	  if(err[d] < min){
	    min = err[d];
	    argmind = d;
	  }else if(err[d] > max){
	    max = err[d];
	    argmaxd = d;
	  }
	}
      }

      {
	if(max + min > 1.0){
	  /** 
	   * min > 1 - max 
	   *  argmaxd is the best axis
	   */
	  (*adaAxis)[t] = argmaxd;
	  (*adaSign)[t] = 1;
	  epsilon = 1 - max;
	}else{
	  /*  argmind is the best axis */
	  (*adaAxis)[t] = argmind;
	  (*adaSign)[t] = 0;
	  epsilon = min;       
	}      
      }
    }

    /* step 3: compute new weithgts */
    {
      (*adaBeta)[t] = epsilon / (1 - epsilon);
      for(i = 0; i < N; i++){
	/* compute data */
	xaxis = adaDataOnTheFlyAxis(h_i[i], h_j[i], kmerFreq, k, (*adaAxis)[t]);
	if(((*adaSign)[t] == 0 && xaxis == y[i]) ||
	   ((*adaSign)[t] == 1 && xaxis != y[i])){
	  w[i] *= (*adaBeta)[t];
	}
      }
    }
    
  }

  free(x);
  free(w);
  free(p);
  free(err);
  return 0;
}

int binarization(const double *source,
		 const double threshold,
		 const unsigned long length,
		 int **target){
  unsigned long i;
  *target = calloc_errchk(length, sizeof(int), "calloc target");
  for(i = 0; i < length; i++){
    (*target)[i] = ((source[i] > threshold) ? 1 : 0);
  }
  return 0;
}

int dump_results(FILE *fp,
		 const unsigned long *adaAxis,
		 const int *adaSign,
		 const double *adaBeta,
		 const unsigned long T,
		 const int k){
  const unsigned long nkmers = 1 << (2 * k);
  unsigned long t;
  char **kmerStrings;

  setKmerStrings(k, &kmerStrings);
  
  for(t = 0; t < T; t++){
    fprintf(fp, "%e\t%s\t%s\t%ld\t%d\n",
	    adaBeta[t],
	    kmerStrings[adaAxis[t] / nkmers],
	    kmerStrings[adaAxis[t] % nkmers],
	    adaAxis[t],
	    adaSign[t]);
  }

  for(t = 0; t < nkmers; t++){
    free(kmerStrings[t]);
  }
  free(kmerStrings);

  return 0;
}

int main_sub(const char *freqFile,
	     const char *hicFile,
	     const int k,
	     const int res,
	     const double threshold,
	     const unsigned long T,
	     const char *outDir){

  int **kmerFreq;
  unsigned long nBin;

  int *h_i, *h_j;
  double *h_mij;
  int *y;
  unsigned long nHic;

  unsigned long *adaAxis;
  int *adaSign;
  double *adaBeta;

  readTableInt(freqFile, "\t", 1 << (2 * k), &kmerFreq, &nBin);
  readHic(hicFile, res, &h_i, &h_j, &h_mij, &nHic);
  binarization(h_mij, threshold, nHic, &y);
  free(h_mij);

  adaboostLearn((const int *)y, 
		(const int *)h_i, (const int*)h_j, 
		(const int **)kmerFreq, 
		k, T, (const unsigned long)nHic, 
		1 << (4 * k),
		&adaAxis, &adaSign, &adaBeta);
		
  dump_results(stderr, adaAxis, adaSign, adaBeta, T, k);

  //char **kmerStrings;
  //setKmerStrings(k, &kmerStrings);

#if 0
  FILE *fpin, *fpout;
  int **feature;
  const unsigned long pairFeatureDim = 1 << (4 * k);
  const unsigned long featureDim = 1 << (2 * k);
  unsigned long featureLen, hicLine = 0, hicHigh = 0;
  //  unsigned long *freqBackGround, *freqHighContact;
  char outFile[F_NAME_LEN], buf[BUF_SIZE], mijbuf[128];
  char **kmerStrings;
  double mij;
  int i, j;   /* index for genomic bins */
  unsigned long l, m;  /* index for k-mers */
#endif


#if 0
  freqBackGround = calloc_errchk(sizeof(unsigned long), pairFeatureDim, "calloc freqBackGround");
  freqHighContact = calloc_errchk(sizeof(unsigned long), pairFeatureDim, "calloc freqHighContact");


  sprintf(outFile, "%s.k%d.t%f.kmerPairOdds", hicFile, k, threshold);

#endif


#if 0
  if((fpin = fopen(hicFile, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    hicFile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fpin) != NULL){
    sscanf(buf, "%d\t%d\t%s", &i, &j, (char *)(&mijbuf));
    mij = strtod(mijbuf, NULL);
    i /= res;
    j /= res;

    if(feature[i] != NULL && feature[j] != NULL){
      if(mij >= threshold){
	for(l = 0; l < featureDim; l++){
	  for(m = 0; m < featureDim; m++){
	    freqHighContact[l * featureDim + m] += feature[i][l] * feature[j][m];
	    freqBackGround[l * featureDim + m] += feature[i][l] * feature[j][m];
	  }
	}
	hicHigh++;
	hicLine++;
      }else{
	for(l = 0; l < featureDim; l++){
	  for(m = 0; m < featureDim; m++){
	    freqBackGround[l * featureDim + m] += feature[i][l] * feature[j][m];
	  }
	}
	hicLine++;
      }    
      if((hicLine % 1000000) == 0){
	fprintf(stderr, "proceeded %ld lines\n", hicLine);
      }
    }
  }

  printf("k-mer pair odds file is saved to %s\n", outFile);

  fclose(fpin);

  if((fpout = fopen(outFile, "w")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    outFile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  for(l = 0; l < featureDim; l++){
    for(m = 0; m < featureDim; m++){
      fprintf(fpout, "%f\t%s\t%s\t%ld\t%e\t%e\n", 
	      (1.0 * freqHighContact[l * featureDim + m] * hicLine) / (1.0 * freqBackGround[l * featureDim + m] * hicHigh),
	      kmerStrings[l],
	      kmerStrings[m],
	      l * featureDim + m,	   
	      1.0 * freqHighContact[l * featureDim + m] / (hicHigh * res * res),
	      1.0 * freqBackGround[l * featureDim + m] / (hicLine * res * res)
	      );
    }
  }

  fclose(fpout);

  free(freqBackGround);
  free(freqHighContact);
#endif

  return 0;
}

int show_usage(const char *progName){
  printf("usage:\n\n");
  printf("example: %s -f %s -H %s -k %d -r %d, -t %f -o %s\n", 
	 progName,
	 "../data/hoge.dat",
	 "../data/hic.dat",
	 5,
	 1000,
	 10.0,
	 "./");
  return 0;
}

int check_params(const char *freqFile,
		 const char *hicFile,
		 const int k,
		 const int res,
		 const double threshold,
		 const unsigned long T,
		 const char *outDir,
		 const char *progName){
  int errflag = 0;

  if(freqFile == NULL){
    fprintf(stderr, "input frequency file is not specified\n");
    errflag++;
  }else{
    printf("fasta file:    %s\n", freqFile);
  }

  if(hicFile == NULL){
    fprintf(stderr, "Hi-C data directory is not specified\n");
    errflag++;
  }else{
    printf("Hi-C data dir: %s\n", hicFile);
  }

  if(k <= 0){
    fprintf(stderr, "k is not specified\n");
    errflag++;
  }else{
    printf("k: %d\n", k);
  }

  if(res <= 0){
    fprintf(stderr, "resolution is not specified\n");
    errflag++;
  }else{
    printf("resolution: %d\n", res);
  }

  if(threshold <= 0){
    fprintf(stderr, "threshold is not specified\n");
    errflag++;
  }else{
    printf("threshold: %e\n", threshold);
  }

  if(T <= 0){
    fprintf(stderr, "number of weak lerners is not specified\n");
    errflag++;
  }else if(1 << (4 * k) < T){
    fprintf(stderr, "number of weak lerners(T = %ld) exceeds 16^k (%d)\n", T, 1 << (4 * k));
    errflag++;
  }else{
    printf("num. of weak lerners: %ld\n", T);
  }

  if(outDir == NULL){
    fprintf(stderr, "output directory is not specified\n");
    errflag++;
  }else{
    printf("Output dir:    %s\n", outDir);
  }

  if(errflag > 0){
    show_usage(progName);
    exit(EXIT_FAILURE);
  }

  return 0;
}

int main(int argc, char **argv){
  char *freqFile = NULL, *hicFile = NULL, *outDir = NULL;
  int k = 0, res = 0;
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
    {0, 0, 0, 0}
  };

  while((opt = getopt_long(argc, argv, "hvf:H:k:r:t:T:o:",
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
    }
  }

  check_params(freqFile, hicFile, k, res, threshold, T, outDir, argv[0]);

  main_sub(freqFile, hicFile, k, res, threshold, T, outDir);
 
  return 0;
}
