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
	    
  printf("loaded genome sequence %s (length = %d)\n",
	 sequenceHead, strlen(sequence));

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
	  
/* computation of P and q */

int computeParamsNormExp(const int k, 
			 const int chr,
			 const int res,
			 const int minBinDist,
			 const int maxBinDist,
			 const int **feature, 
			 const char *hicFileRaw, 
			 const double *normalize,
			 const double *expected,
			 const char*outDir){
  FILE *fp;
  char buf[100], mijbuf[20];
  int i, j, l, m, n;
  double mij;
  double *P, *q, *dij;
  const long Psize = (1 << (8 * k - 1)) + (1 << (4 * k - 1));
  const long qsize = 1 << (4 * k);
  const long fsize = 1 << (2 * k);
  const double Pcoef = 1.0 / ((res + k - 1) * (res + k - 1));
  char *outFileP, *outFileq;

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

      /* normalize & O/E conversion */
      mij = strtod(mijbuf, NULL) / (normalize[i] * normalize[j] * expected[abs(i - j)]);

      if(!isnan(mij) && !isinf(mij)){
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

  fclose(fp);

  setOutFileName(outDir, res, chr, k, &outFileP, &outFileq);
  save2file(outFileP, P, Psize);
  save2file(outFileq, q, qsize);

  free(outFileP);
  free(outFileq);
  

  free(P);
  free(q);
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

  fprintf(fp, "%d\n", length);

  for(i = 0; i < length; i++){
    fprintf(fp, "%e\n", ary[i]);
  }
  
  fclose(fp);

  return 0;
}


int main(void){
  char *fastaName = "../data/GRCh37.ch21.fasta";
  char *hicDir = "../data/GM12878_combined/";
  char *outDir = "../out/";
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

  sequencePrep(k, binSize, fastaName, &feature, &binNum);

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

#define DEBUG 0

