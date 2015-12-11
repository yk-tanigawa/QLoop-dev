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

int readFeature(const char *freqFile,
		const int k,
		int **feature);


int main_sub(const char *freqFile,
	     const char *hicFile,
	     const int k,
	     const double threshold,
	     const char *outDir){
  return 0;
}

int show_usage(const char *progName){
  printf("usage:\n\n");
  printf("example: %s -f %s -H %s -k %d -t %f -o %s\n", 
	 progName,
	 "../data/hoge.dat",
	 "../data/hic.dat",
	 5,
	 10.0,
	 "./");
  return 0;
}

int check_params(const char *freqFile,
		 const char *hicFile,
		 const int k,
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
  //  char freqFile[F_NAME_LEN], hicFile[F_NAME_LEN], outDir[F_NAME_LEN];
  char *freqFile = NULL, *hicFile = NULL, *outDir = NULL;
  int k = -1;
  double threshold = -1;

  int opt = 0, opt_idx = 0;
  struct option long_opts[] = {
    {"help",    no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    /* options */
    {"frequency", required_argument, NULL, 'f'},
    {"hic",       required_argument, NULL, 'H'},
    {"k",         required_argument, NULL, 'k'},
    {"threshold", required_argument, NULL, 't'},
    {"out",       required_argument, NULL, 'o'},
    {0, 0, 0, 0}
  };

  while((opt = getopt_long(argc, argv, "hvf:H:k:t:o:",
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
      case 't': /* threshold */
	threshold = atof(optarg);
	break;
      case 'o': /* out */
	outDir = optarg;
	break;
    }
  }

  check_params(freqFile, hicFile, k, threshold, outDir, argv[0]);
  main_sub(freqFile, hicFile, k, threshold, outDir);
 
  return 0;
}
