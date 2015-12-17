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

#define F_NAME_LEN 128
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

int readFeature(const char *freqFile,
		const int k,
		int ***feature,
		unsigned long *featureLen){
  FILE *fp;
  char buf[BUF_SIZE];
  char *tok;
  const char *dlim = "\t";
  long i = 0, l = 0;
  const unsigned long featureDim = 1 << (2 * k);

  *featureLen = wc(freqFile);
  
  *feature = calloc_errchk(sizeof(int *), *featureLen, "calloc feature");

  if((fp = fopen(freqFile, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    freqFile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fp)){
    if(buf[0] == '*'){
      (*feature)[l] = NULL;
    }else{
      /** memory allocation */
      (*feature)[l] = calloc_errchk(sizeof(int), featureDim, "calloc feature[l]");
      i = 0;
      tok = strtok(buf, dlim);
      while(tok != NULL && i < featureDim){
	(*feature)[l][i++] = atoi(tok);
	//fprintf(stderr, "%d ", (*feature)[l][i - 1]);
	tok = strtok(NULL, dlim);
      }
      //fprintf(stderr, "\n");
    }
    l++;
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


int main_sub(const char *freqFile,
	     const char *hicFile,
	     const int k,
	     const int res,
	     const double threshold,
	     const char *outDir){
  FILE *fpin, *fpout;
  int **feature;
  const unsigned long pairFeatureDim = 1 << (4 * k);
  const unsigned long featureDim = 1 << (2 * k);
  unsigned long featureLen, hicLine = 0, hicHigh = 0;
  unsigned long *freqBackGround, *freqHighContact;
  char outFile[F_NAME_LEN], buf[BUF_SIZE], mijbuf[128];
  char **kmerStrings;
  double mij;
  int i, j;   /* index for genomic bins */
  long l, m;  /* index for k-mers */

  setKmerStrings(k, &kmerStrings);

  freqBackGround = calloc_errchk(sizeof(unsigned long), pairFeatureDim, "calloc freqBackGround");
  freqHighContact = calloc_errchk(sizeof(unsigned long), pairFeatureDim, "calloc freqHighContact");
  sprintf(outFile, "%s.k%d.t%f.kmerPairOdds", hicFile, k, threshold);

  readFeature(freqFile, k, &feature, &featureLen);

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
		 const char *outDir,
		 const char *progName){

  if(freqFile == NULL){
    fprintf(stderr, "input frequency file is not specified\n");
    show_usage(progName);
    exit(EXIT_FAILURE);
  }else{
    printf("fasta file:    %s\n", freqFile);
  }

  if(hicFile == NULL){
    fprintf(stderr, "Hi-C data directory is not specified\n");
    show_usage(progName);
    exit(EXIT_FAILURE);
  }else{
    printf("Hi-C data dir: %s\n", hicFile);
  }

  if(k < 0){
    fprintf(stderr, "k is not specified\n");
    show_usage(progName);
    exit(EXIT_FAILURE);
  }else{
    printf("k: %d\n", k);
  }

  if(res < 0){
    fprintf(stderr, "resolution is not specified\n");
    show_usage(progName);
    exit(EXIT_FAILURE);
  }else{
    printf("resolution: %d\n", res);
  }

  if(threshold < 0){
    fprintf(stderr, "threshold is not specified\n");
    show_usage(progName);
    exit(EXIT_FAILURE);
  }else{
    printf("t: %f\n", threshold);
  }

  if(outDir == NULL){
    fprintf(stderr, "output directory is not specified\n");
    show_usage(progName);
    exit(EXIT_FAILURE);
  }else{
    printf("Output dir:    %s\n", outDir);
  }

  return 0;
}

int main(int argc, char **argv){
  char *freqFile = NULL, *hicFile = NULL, *outDir = NULL;
  int k = -1, res = -1;
  double threshold = -1;

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
    {"out",       required_argument, NULL, 'o'},
    {0, 0, 0, 0}
  };

  while((opt = getopt_long(argc, argv, "hvf:H:k:r:t:o:",
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
      case 'o': /* out */
	outDir = optarg;
	break;
    }
  }

  check_params(freqFile, hicFile, k, res, threshold, outDir, argv[0]);

  main_sub(freqFile, hicFile, k, res, threshold, outDir);
 
  return 0;
}
