#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
  struct stat stbuf;

  /* file open */
  if((fp = fopen(fastaName, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    fastaName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  /* get sequence header */
  fscanf(fp, "%s", sequenceHead);

 
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
  int i, kmer = 0, bitMask = (1 << (2 * k)) - 1;

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

#if DEBUG
  if((char *dump = calloc(sizeof(char), 1001)) == NULL){
    perror("error(calloc) debug func");
    exit(EXIT_FAILURE);
  }

#endif

  /* allocate memory space for sequence */
  fileSize = getFileSize(fastaName);
  if((sequence = calloc(1, fileSize)) == NULL){
    perror("error(calloc) sequence:");
    exit(EXIT_FAILURE);
  }
  
  readFasta(fastaName, sequenceHead, sequence);	    
	    
  printf("loaded genome sequence %s (length = %d)\n",
	 sequenceHead, strlen(sequence));

  computeFeature(k, binSize, sequence, feature, binNum);

#if DEBUG
  strncpy(dump, &(sequence[9412000]), 1000);
  printf("%s\n", dump);
#endif

  free(sequence);

  return 0;
}

/* Hi-C contact frequency matrix  */

char *res2str(const int res){
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
  
  char fileHead[100];
  sprintf(fileHead, "%s%s_resolution_intrachromosomal/chr%d/MAPQGE30/chr%d_%s.",
	  hicDir, res2str(res), chr, chr, res2str(res));
    
  sprintf(*hicFileRaw, "%s%s", fileHead, "RAWobserved");

  if(normalize != NULL){
    sprintf(*hicFileNormalize, "%s%s%s", fileHead, normalize, "norm");
  }

  if(expected != NULL){
    sprintf(*hicFileExpected, "%s%s%s", fileHead, expected, "expected");
  }
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
	    const int maxDist,
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
	  
/* main */

int computeParamsNormExp(const int k, 
			 const int res,
			 const int minBinDist,
			 const int maxBinDist,
			 const int **feature, 
			 const char *hicFileRaw, 
			 const double *normalize,
			 const double *expected){
  FILE *fp;
  char buf[100], mijbuf[20];
  int i, j, l, m, n;
  double mij;

  double *P, *q, *dij;
  long Psize = (1 << (8 * k - 1)) + (1 << (4 * k - 1));
  long qsize = 1 << (4 * k);
  long fsize = 1 << (2 * k);
  double Pcoef = 1.0 / ((res + k - 1) * (res + k - 1));

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

  while(fgets(buf, 100, fp) != NULL){
    sscanf(buf, "%d\t%d\t%s", &i, &j, &mijbuf);
    i /= res;
    j /= res;
    if(minBinDist <= abs(i - j) && 
       abs(i - j) <= maxBinDist &&
       feature[i] != NULL &&
       feature[j] != NULL){
      mij = strtod(mijbuf, NULL) / (normalize[i] * normalize[j] * expected[abs(i - j)]);
      if(!isnan(mij) && !isinf(mij)){
	printf("%d %d %.20f\n", i, j, mij);
	/* compute dij */
	for(m = 0; m < fsize; m++){
	  for(l = 0; l < fsize; l++){
	    dij[m * fsize + l] = feature[j][m] * feature[i][l];
	  }
	}

	#if 0
	for(m = 0; m < qsize; m++){
	  printf("%d ", dij[m]);
	}
	printf("\n");
	#endif

	/* add to q */
	for(m = 0; m < qsize; m++){
	  q[m] += mij * dij[m];
	}
	
	#if 0
	for(m = 0; m < qsize; m++){
	  printf("%f ", q[m]);
	}
	printf("\n");
	#endif

	n = 0; /* n is an index for P */
	for(m = 0; m < qsize; m++){
	  for(l = m; l < qsize; l++){
	    P[n++] += Pcoef * dij[m] * dij[l];
	  }
	}

      }
    }
  }

  printf("\n");


  for(m = 0; m < qsize; m++){
    printf("%f ", q[m]);
  }
  printf("\n");
  printf("\n");
  for(m = 0; m < Psize; m++){
    printf("%f ", P[m]);
  }
  printf("\n");


  fclose(fp);


  free(P);
  free(q);
}



int main(void){
  char *fastaName = "../data/GRCh37.ch21.fasta";
  char *hicDir = "../data/GM12878_combined/";
  int k = 1;
  int binSize = 1000;
  int res = 1000;
  int chr = 21;
  int minDist = 10000;
  int maxDist = 1000000;
  char *normalizeMethod = "KR";
  char *expectedMethod = "KR";

  int **feature;
  char *hicFileRaw;
  double *normalize;
  double *expected;
  int binNum;

  int i;

#if DEBUG
  int l;
#endif

  sequencePrep(k, binSize, fastaName, &feature, &binNum);

  hicPrep(hicDir, res, chr, maxDist, binNum,
	  normalizeMethod, expectedMethod,
	  &hicFileRaw, &normalize, &expected);

  computeParamsNormExp(k, res, minDist / res, maxDist / res, 
		       (const int **)feature, hicFileRaw, normalize, expected);

  printf("%s\n", hicFileRaw);


  /* free the allocated memory and exit */
  free(normalize);
  free(expected);

  for(i = 0; i < binNum; i++){
    if(feature[i] != NULL){
        
      #if DEBUG
      printf("%d :", i);
      for(l = 0; l < (1 << (2 * k)); l++){
	printf("%d ", feature[i][l]);
      }
      printf("\n");
      #endif
        
      free(feature[i]);
    }
  }
  
  free(feature);

  return 0;
}

#define DEBUG 0

