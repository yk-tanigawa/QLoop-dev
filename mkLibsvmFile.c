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
		long *featureLen){
  FILE *fp;
  char buf[BUF_SIZE];
  char *tok;
  char *dlim = "\t";
  long i = 0, l = 0;
  long featureDim = 1 << (2 * k);

  *featureLen = wc(freqFile);
  
  if((*feature = calloc(sizeof(int *), *featureLen)) == NULL){
    perror("error(calloc) feature");
    exit(EXIT_FAILURE);
  }

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
      if(((*feature)[l] = calloc(sizeof(int), featureDim)) == NULL){
	fprintf(stderr, "error(calloc) fature[%ld]\n%s\n",
		l, strerror(errno));
	exit(EXIT_FAILURE);
      }

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

int main_sub(const char *freqFile,
	     const char *hicFile,
	     const int k,
	     const int res,
	     const double threshold,
	     const char *outDir){
  FILE *fpin, *fpout;
  int **feature;
  long featureLen, featureDim = 1 << (2 * k);
  char outFile[F_NAME_LEN], buf[BUF_SIZE], mijbuf[128];
  double mij;
  int i, j;   /* index for genomic bins */
  long l, m;  /* index for k-mers */
  long z;

  sprintf(outFile, "%s.libsvm", hicFile);

  readFeature(freqFile, k, &feature, &featureLen);
  //fprintf(stderr, "%ld\n", featureLen);

  if((fpin = fopen(hicFile, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    hicFile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  if((fpout = fopen(outFile, "w")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    outFile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, BUF_SIZE, fpin) != NULL){
    sscanf(buf, "%d\t%d\t%s", &i, &j, (char *)(&mijbuf));
    mij = strtod(mijbuf, NULL);
    i /= res;
    j /= res;
    if(feature[i] != NULL && feature[j] != NULL){
      fprintf(fpout, "%d", ((mij >= threshold)? 1 : 0));
      z = 1;
      for(l = 0; l < featureDim; l++){
	for(m = 0; m < featureDim; m++){
	  fprintf(fpout, " %ld:%d", z++, 
		  ((feature[i][l] * feature[j][m] > 0) ? 1 : 0));
	}
      }
      fprintf(fpout, "\n");
    }
  }

  printf("libSVM format file is saved to %s\n", outFile);

  fclose(fpout);
  fclose(fpin);

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
