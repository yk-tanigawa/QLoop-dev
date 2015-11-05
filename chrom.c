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

#define BENCH 1

#if BENCH

#include <sys/time.h>

int dumpResults(struct timeval *t, long *l){
  time_t diff_sec = (t[1].tv_sec - t[0].tv_sec);
  suseconds_t diff_usec = (t[1].tv_usec - t[0].tv_usec);
  double diff = (double)diff_sec + ((double)diff_usec / 1000000);
  printf("# of data points : \t%ld of %ld\n", l[1], l[0]);
  printf("computation time : \t%.6f [sec.]\n", diff);
  printf("time per data point: \t%e [sec. / data point]\n", diff / l[1]);
  return 0;
}

#endif


/* genomic sequence */

long getFileSize(const char *fname){
  int fd;
  struct stat stbuf;

  if((fd = open(fname, O_RDONLY)) == -1){
    fprintf(stderr, "error: open %s\n%s\n", 
	    fname, strerror(errno));      
    exit(EXIT_FAILURE);
  }

  if(fstat(fd, &stbuf) == -1){
    fprintf(stderr, "error: fstat %s\n%s\n",
	    fname, strerror(errno));
    exit(EXIT_FAILURE);
  }

  close(fd);

  return stbuf.st_size;
}

int readFasta(const char *fastaName, 
	      char *sequenceHead,
	      char *sequence){
  FILE *fp;
  char buf[256];

  /* file open */
  if((fp = fopen(fastaName, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    fastaName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  /* get sequence header */
  fscanf(fp, "%s", sequenceHead);
 
  /* get sequence body */
  while(fscanf(fp, "%s", buf) != EOF) {
    strncpy(sequence, buf, strlen(buf));
    sequence += strlen(buf);
  }

  fclose(fp);
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

  if((*feature_i = calloc(sizeof(int), 1 << (2 * k))) == NULL){
    perror("error(calloc) feature[i]");
    exit(EXIT_FAILURE);
  }
  
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

/* Hi-C contact frequency matrix  */

inline char *res2str(const int res){
  switch(res){
    case 1000:
      return "1kb";
    default:
      fprintf(stderr, "resolution size %d is not supported\n", res);
     exit(EXIT_FAILURE);
  }
}

int hicFileNames(const char *hicDir,
		 const int res,
		 const int chr,
		 const char *normalize,
		 const char *expected,
		 char **hicFileRaw,
		 char **hicFileNormalize,
		 char **hicFileExpected){
  char fileHead[100];

  if((*hicFileRaw = calloc(sizeof(char), 100)) == NULL){
    perror("error(calloc) hicFileRaw");
    exit(EXIT_FAILURE);
  }
  if((*hicFileNormalize = calloc(sizeof(char), 100)) == NULL){
    perror("error(calloc) hicFileNormalize");
    exit(EXIT_FAILURE);
  }
  if((*hicFileExpected = calloc(sizeof(char), 100)) == NULL){
    perror("error(calloc) hicFileExpected");
    exit(EXIT_FAILURE);
  }
  
  sprintf(fileHead, "%s%s_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%s.",
	  hicDir, res2str(res), chr, chr, res2str(res));
    
  sprintf(*hicFileRaw, "%s%s", fileHead, "RAWobserved");

  if(normalize != NULL){
    sprintf(*hicFileNormalize, "%s%s%s", fileHead, normalize, "norm");
  }

  if(expected != NULL){
    sprintf(*hicFileExpected, "%s%s%s", fileHead, expected, "expected");
  }

  return 0;
}

int readDouble(const char *fileName, double **array, const int lineNum){
  FILE *fp;
  int i = 0;
  char buf[50];

  if((*array = calloc(sizeof(double), lineNum)) == NULL){
    perror("error(calloc) readDouble");
    exit(EXIT_FAILURE);
  }

  if((fp = fopen(fileName, "r")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    fileName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fgets(buf, 50, fp) != NULL){
    if(i < lineNum){
      (*array)[i++] = strtod(buf, NULL);
    }else{
      break;
    }
  }

  fclose(fp);

  return 0;
}

int hicPrep(const char *hicDir, 
	    const int res, 
	    const int chr, 
	    const long maxDist,
	    const int binNum,
	    const char *normalizeMethod,
	    const char *expectedMethod,
	    char **hicFileRaw,
	    double **normalize, 
	    double **expected){	    

  char *hicFileNormalize;
  char *hicFileExpected;

  hicFileNames(hicDir, res, chr, normalizeMethod, expectedMethod,
	       hicFileRaw,
	       &hicFileNormalize,
	       &hicFileExpected);

  if(normalizeMethod != NULL){
    readDouble(hicFileNormalize, normalize, binNum);
    free(hicFileNormalize);
  }

  if(expectedMethod != NULL){
    readDouble(hicFileExpected, expected, maxDist / res);
    free(hicFileExpected);
  }

  return 0;
}

int setOutFileName(const char *outDir,
		   const int res,
		   const int chr,
		   const int k,		   
		   char **outFileP,
		   char **outFileq){
  *outFileP = calloc(sizeof(char), 100);
  *outFileq = calloc(sizeof(char), 100);
  sprintf(*outFileP, "%sres%d.chr%d.k%d.P.out", outDir, res, chr, k);
  sprintf(*outFileq, "%sres%d.chr%d.k%d.q.out", outDir, res, chr, k);
  return 0;
}

int save2file(const char *fileName,
	      const double *ary,
	      const long length){
  FILE *fp;
  long i = 0;

  if((fp = fopen(fileName, "w")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    fileName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  fprintf(fp, "%ld\n", length);

  for(i = 0; i < length; i++){
    fprintf(fp, "%e\n", ary[i]);
  }
  
  fclose(fp);

  return 0;
}
	  
/* computation of P and q */

int computeParamsNormExp(const int k, 
			 const int chr,
			 const int res,
			 const long minBinDist,
			 const long maxBinDist,
			 const int **feature, 
			 const char *hicFileRaw, 
			 const double *normalize,
			 const double *expected,
			 const char*outDir){
  FILE *fp;
  int i, j, l, m, n;
  double mij;
  double *P, *q, *dij;
  char *outFileP, *outFileq;
  char buf[100], mijbuf[20];
#if BENCH
  struct timeval tv[2];
  long benchLines[2];
#endif
  const long Psize = (1 << (8 * k - 1)) + (1 << (4 * k - 1));
  const long qsize = 1 << (4 * k);
  const long fsize = 1 << (2 * k);
  const double Pcoef = 1.0 / ((res + k - 1) * (res + k - 1));

  /* Hi-C file open */
  if((fp = fopen(hicFileRaw, "r")) == NULL){
    fprintf(stderr, "error: fopen %s\n%s\n",
	    hicFileRaw, strerror(errno));
    exit(EXIT_FAILURE);
  }

  /* memory allocation */
  if((P = calloc(sizeof(double), Psize)) == NULL){
    perror("error(calloc) P");
    exit(EXIT_FAILURE);
  }

  if((q = calloc(sizeof(double), qsize)) == NULL){
    perror("error(calloc) P");
    exit(EXIT_FAILURE);
  }

  if((dij = calloc(sizeof(double), qsize)) == NULL){
    perror("error(calloc) dij");
    exit(EXIT_FAILURE);
  }


  printf("Hi-C file %s is open\n", hicFileRaw);

#if BENCH
  benchLines[0] = benchLines[1] = 0;
  if(gettimeofday(&(tv[0]), NULL) == -1){
    perror("gettimeofday");
    exit(EXIT_FAILURE);
  }
#endif

  /* compute P and q */

  while(fgets(buf, 100, fp) != NULL){

#if BENCH
    benchLines[0]++;
#endif

    sscanf(buf, "%d\t%d\t%s", &i, &j, (char *)(&mijbuf));
    i /= res;
    j /= res;
    if(minBinDist <= abs(i - j) && 
       abs(i - j) <= maxBinDist &&
       feature[i] != NULL &&
       feature[j] != NULL){

      /* normalize & O/E conversion */
      mij = strtod(mijbuf, NULL) / (normalize[i] * normalize[j] * expected[abs(i - j)]);

      if(!isnan(mij) && !isinf(mij)){

#if BENCH
	benchLines[1]++;
#endif

	/* compute dij */
	for(m = 0; m < fsize; m++){
	  for(l = 0; l < fsize; l++){
	    dij[m * fsize + l] = feature[j][m] * feature[i][l];
	  }
	}

	/* add to q */
	for(m = 0; m < qsize; m++){
	  q[m] += mij * dij[m];
	}
	
	/* add to P */
	n = 0; /* n is an index for P */
	for(m = 0; m < qsize; m++){
	  for(l = m; l < qsize; l++){
	    P[n++] += Pcoef * dij[m] * dij[l];
	  }
	}
      }
    }
  }

#if BENCH
  if(gettimeofday(&(tv[1]), NULL) == -1){
    perror("gettimeofday");
    exit(EXIT_FAILURE);
  }
#endif

  fclose(fp);

  setOutFileName(outDir, res, chr, k, &outFileP, &outFileq);

  printf("Writing P to %s\n", outFileP);
  printf("Writing q to %s\n", outFileq);

  save2file(outFileP, P, Psize);
  save2file(outFileq, q, qsize);

  free(outFileP);
  free(outFileq);
  
  free(P);
  free(q);

#if BENCH
  dumpResults(tv, benchLines);
#endif

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

  sequencePrep(k, res, fastaName, &feature, &binNum);

  hicPrep(hicDir, res, chr, maxDist, binNum,
	  normalizeMethod, expectedMethod,
	  &hicFileRaw, &normalize, &expected);

  computeParamsNormExp(k, chr, res, minDist / res, maxDist / res, 
		       (const int **)feature, hicFileRaw, normalize, expected, outDir);

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
