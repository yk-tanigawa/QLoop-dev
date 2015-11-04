#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

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
  char *dump = calloc(sizeof(char), 1001);
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
  *hicFileRaw =
    "../data/GM12878_combined/1kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_1kb.RAWobserved";
  *hicFileNormalize = 
    "../data/GM12878_combined/1kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_1kb.KRnorm";
  *hicFileExpected =
    "../data/GM12878_combined/1kb_resolution_intrachromosomal/chr21/MAPQGE30/chr21_1kb.KRexpected";
}

int readDouble(const char *fileName, double **array, const int lineNum){
  *array = calloc(sizeof(double), lineNum);
  FILE *fp;
  double buf;
  int i = 0;
  
  if((fp = fopen(fileName, "r")) == NULL){
    fprintf(stderr, "error: fdopen %s\n%s\n",
	    fileName, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while(fscanf(fp, "%lf", &buf) != EOF) {
    if(i < lineNum){
      (*array)[i++] = buf;
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

  if(hicFileNormalize != NULL){
    readDouble(hicFileNormalize, normalize, maxDist / res);
  }
  if(hicFileExpected != NULL){
    readDouble(hicFileExpected, expected, binNum);
  }

  return 0;

}
	  

int main(void){
  char *fastaName = "../data/GRCh37.ch21.fasta";
  char *hicDir = "../data/GM12878_combined/";
  int k = 5;
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

