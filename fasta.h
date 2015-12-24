#ifndef __FASTA_H__
#define __FASTA_H__

#include "constant.h"
#include "cmd_args.h"
#include "mywc.h"
#include "calloc_errchk.h"
#include "diffSec.h"

int fasta_read(const char *fasta_file, 
	       char **seq_head,
	       char **seq,
	       unsigned long *seq_len){
  unsigned long i = 0;

  *seq_len = mywc_b(fasta_file);
  {
    FILE *fp;
    char buf[BUF_SIZE];

    /* file open */
    if((fp = fopen(fasta_file, "r")) == NULL){
      fprintf(stderr, "error: fdopen %s\n%s\n",
	      fasta_file, strerror(errno));
      exit(EXIT_FAILURE);
    }

    /* get sequence header */
    {
      *seq_head = calloc_errchk(FASTA_HEADER_LEN, sizeof(char), "seq_head");
      fscanf(fp, "%s", *seq_head);

    }

    /* get sequence body */ 
    {
      *seq = calloc_errchk(*seq_len, sizeof(char), "seq");
      while(fscanf(fp, "%s", buf) != EOF && i < *seq_len){
	strncpy(*seq, buf, strlen(buf));
	*seq += strlen(buf);
	i += strlen(buf);
      }
    }
    fclose(fp);
  }
  *seq_len = i;

  return 0;
}

inline int c2i(const char c){
  switch(c){
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
   default:
     fprintf(stderr, "input genomic sequence contains unknown char : %c\n", c);
     exit(EXIT_FAILURE);
  }
}

int count_kmer_freq(const command_line_arguements *cmd_args){
  char *seq_head, *seq;
  unsigned long seq_len;

  struct timeval ts, tg;
  gettimeofday(&ts, NULL);

  fasta_read(cmd_args->fasta_file, 
	     &seq_head, &seq, &seq_len);

  gettimeofday(&tg, NULL);
  printf("%f\n", diffSec(ts, tg));

  return 0;
}

#if 0
int computeFeatureSub(const int k,
		      const int left,
		      const int right,
		      const char *sequence,
		      int **feature_i){
  int i, kmer = 0;
  const int bitMask = (1 << (2 * k)) - 1;

  /* check if the sequence contains N */
  for(i = left; i < right; i++){
    if(sequence[i] == 'N' || sequence[i] == 'n'){
      return -1;
    }
  }

  *feature_i = calloc_errchk(1 << (2 * k), sizeof(unsigned int), 
			     "error(calloc) feature[i]");
  
  for(i = left; i < left + k - 1; i++){
    kmer = (kmer << 2);
    kmer += c2i(sequence[i]);
  }

  /* count k-mer frequency */
  for(i = left + k - 1; i < right; i++){
    kmer = (kmer << 2);
    kmer += c2i(sequence[i]);
    kmer = kmer & bitMask;
    (*feature_i)[kmer] += 1;
  }

  return 0;
}
		       
int computeFeature(const int k,
		   const int binSize,
		   const char *sequence,
		   int ***feature,
		   int *binNum){
  int i;
  *binNum = strlen(sequence) / binSize;

  if((*feature = calloc(sizeof(int *), *binNum)) == NULL){
    perror("error(calloc) feature:");
    exit(EXIT_FAILURE);
  }
  
  for(i = 0; i < *binNum; i++){
    computeFeatureSub(k, i * binSize, (i + 1) * binSize + k, sequence, &((*feature)[i]));
  }
  return 0;
}

int sequencePrep(const int k,		 
		 const int binSize,
		 const char *fastaName, 
		 int ***feature,
		 int *binNum){
  /* parameters */
  char sequenceHead[80];
  char *sequence;
  long fileSize;

  /* allocate memory space for sequence */
  fileSize = getFileSize(fastaName);
  if((sequence = calloc(1, fileSize)) == NULL){
    perror("error(calloc) sequence:");
    exit(EXIT_FAILURE);
  }
  
  readFasta(fastaName, sequenceHead, sequence);	    
	    
  printf("loaded genome sequence %s (length = %ld)\n",
	 sequenceHead, (long)strlen(sequence));

  computeFeature(k, binSize, sequence, feature, binNum);

  free(sequence);

  return 0;
}

int featureSave(const char *outDir,
		const int k,
		const int res,
		const int chr,
		const int binNum,
		const int **feature){
  FILE *fp;
  long i, j, jMax;
  char fileName[100];

  jMax = (1 << (2 * k));

  sprintf(fileName, "%sk%d.res%d.chr%d.dat",
	  outDir, k, res, chr);


  if((fp = fopen(fileName, "w")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    fileName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  for(i = 0; i < binNum; i++){
    if(feature[i] == NULL){
      fprintf(fp, "*\n");
    }else{
      for(j = 0; j < jMax; j++){
	fprintf(fp, "%d\t", feature[i][j]);
      }
      fprintf(fp, "\n");
    }
  }
  
  fclose(fp);

  return 0;
}



int main_sub(const char *fastaName,
	     const char *hicDir,
	     const char *outDir,
	     const int k,
	     const int res,
	     const int chr,
	     const long minDist,
	     const long maxDist,
	     const char *normalizeMethod,
	     const char *expectedMethod){

  int **feature;
  char *hicFileRaw;
  double *normalize;
  double *expected;
  int binNum;
  int i;

  /* compute k-mer frequencies for each bin */
  sequencePrep(k, res, fastaName, &feature, &binNum);

  /* load Hi-C normalize vectors */
  hicNormPrep(hicDir, res, chr, maxDist, binNum,
	      normalizeMethod, expectedMethod,
	      &hicFileRaw, &normalize, &expected);

  /* normalize Hi-C data and write it into a file */
  HicPrep(k, res, chr, normalizeMethod, expectedMethod,
	  minDist / res, maxDist / res, 
	  (const int **)feature, hicFileRaw, 
	  normalize, expected, outDir);

  /* write feature vector */
  featureSave(outDir, k, res, chr, binNum, (const int **)feature);


  /* free the allocated memory and exit */
  free(normalize);
  free(expected);

  for(i = 0; i < binNum; i++){
    if(feature[i] != NULL){                
      free(feature[i]);
    }
  }
  
  free(feature);

  return 0;
}

int check_params(const char *fastaName,
		 const char *hicDir,
		 const char *outDir,
		 const int k,
		 const int res,
		 const int chr,
		 const long minDist,
		 const long maxDist,
		 const char *normalizeMethod,
		 const char *expectedMethod){
  if(fastaName == NULL){
    fprintf(stderr, "input fasta file is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("fasta file:    %s\n", fastaName);
  }
  if(hicDir == NULL){
    fprintf(stderr, "Hi-C data directory is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("Hi-C data dir: %s\n", hicDir);
  }
  if(outDir == NULL){
    fprintf(stderr, "output directory is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("Output dir:    %s\n", outDir);
  }
  if(k == -1){
    fprintf(stderr, "k is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("k: %d\n", k);
  }
  if(res == -1){
    fprintf(stderr, "resolution is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("resolution: %d\n", res);
  }
  if(chr == -1){
    fprintf(stderr, "chromosome number is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("chromosome: %d\n", chr);
  }
  if(minDist == -1){
    fprintf(stderr, "minimum distance is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("min distance: %ld\n", minDist);
  }
  if(maxDist == -1){
    fprintf(stderr, "maximum distance is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("Max distance: %ld\n", maxDist);
  }
  if(normalizeMethod == NULL){
    fprintf(stderr, "warning: normalize method is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("Normalization: %s\n", normalizeMethod);
  }
  if(expectedMethod == NULL){
    fprintf(stderr, "warning: expected value calculation method is not specified\n");
    exit(EXIT_FAILURE);
  }else{
    printf("Expectation:   %s\n", expectedMethod);
  }
  return 0;
}

int main(int argc, char **argv){
  char *fastaName = NULL; //"../data/GRCh37.ch21.fasta";
  char *hicDir = NULL; //"../data/GM12878_combined/";
  char *outDir = NULL; //"../out/";
  int k = -1; //1
  int res = -1; //1000;
  int chr = -1; //21;
  long minDist = -1; //10000;
  long maxDist = -1; //1000000;
  char *normalizeMethod = NULL; //"KR";
  char *expectedMethod = NULL; //"KR";

  struct option long_opts[] = {
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    /* options */
    {"fasta", required_argument, NULL, 'f'},
    {"hic", required_argument, NULL, 'H'},
    {"out", required_argument, NULL, 'o'},
    {"k", required_argument, NULL, 'k'},
    {"res", required_argument, NULL, 'r'},
    {"chr", required_argument, NULL, 'c'},
    {"min", required_argument, NULL, 'm'},
    {"max", required_argument, NULL, 'M'},
    {"norm", required_argument, NULL, 'n'},
    {"expected", required_argument, NULL, 'e'},
    {0, 0, 0, 0}
  };

  int opt = 0;
  int opt_idx = 0;

  while((opt = getopt_long(argc, argv, "hvf:H:o:k:r:c:m:M:n:e:", long_opts, &opt_idx)) != -1){
    switch (opt){
      case 'h': /* help */
	printf("usage:\n\n");
	printf("example: ./chrom -f ../data/GRCh37.ch21.fasta -H ../data/GM12878_combined/ -o ../out/ -k1 -r1000 -c21 -m10000 -M1000000 -n\"KR\" -e\"KR\"\n");
	exit(EXIT_SUCCESS);
      case 'v': /* version*/
	printf("version: 0.10\n");
	exit(EXIT_SUCCESS);
      case 'f': /* fasta */
	fastaName = optarg;
	break;
      case 'H': /* hic */
	hicDir = optarg;
	break;
      case 'o': /* out */
	outDir = optarg;
	break;
      case 'k': /* k */
	k = atoi(optarg);
	break;
      case 'r': /* res */
	res = atoi(optarg);
	break;
      case 'c': /* chr */
	chr = atoi(optarg);
	break;
      case 'm': /* min */
	minDist = atol(optarg);
	break;
      case 'M': /* max */
	maxDist = atol(optarg);
	break;
      case 'n': /* norm */
	normalizeMethod = optarg;
	break;
      case 'e': /* expected */
	expectedMethod = optarg;
	break;
    }
  }

  check_params(fastaName, hicDir, outDir, k, res, chr, minDist, maxDist, normalizeMethod, expectedMethod);
  main_sub(fastaName, hicDir, outDir, k, res, chr, minDist, maxDist, normalizeMethod, expectedMethod);

  return 0;
}
#endif

#endif
