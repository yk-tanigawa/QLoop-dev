#ifndef __FASTA_H__
#define __FASTA_H__

#include "constant.h"
#include "mywc.h"
#include "calloc_errchk.h"
#include "cmd_args.h"

/**
 * This header file contains some functions to perform the following tasks
 * - read FASTA format file and store seqeunce to memory 
 * - compute k-mer frequencies for bins
 */

int fasta_read(const char *, char **, char **, unsigned long *);
int c2i(const char);
int set_kmer_freq_odds(const cmd_args *, double ***);
		  


/* read fasta file */
int fasta_read(const char *fasta_file, 
	       char **seq_head,
	       char **seq,
	       unsigned long *seq_len){

  *seq_len = mywc_b(fasta_file);

  {
    FILE *fp;
    char buf[BUF_SIZE];

    /* file open */
    if((fp = fopen(fasta_file, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      fasta_file, strerror(errno));
      exit(EXIT_FAILURE);
    }

    /* get sequence header */
    {
      unsigned int i;
      *seq_head = calloc_errchk(FASTA_HEADER_LEN, sizeof(char), "seq_head");
      fscanf(fp, "%s", *seq_head);
      for(i = 0; i < strlen(*seq_head) - 1; i++){
	(*seq_head)[i] = (*seq_head)[i + 1];
      }
      (*seq_head)[i] = '\0';
    }

    /* get sequence body */ 
    {
      unsigned long i = 0;
      *seq = calloc_errchk(*seq_len, sizeof(char), "seq");
      while(fscanf(fp, "%s", buf) != EOF && i < *seq_len){	
	strncpy(&((*seq)[i]), buf, strlen(buf));
	i += strlen(buf);
      }
      *seq_len = i;
      (*seq)[i++] = '\0';
      if((*seq = realloc(*seq, i * sizeof(char))) == NULL){
	fprintf(stderr, "realloc: seq\n");
	exit(EXIT_FAILURE);
      }
    }    
    fclose(fp);

  }
  return 0;
}

/* functionn to convert nucleotide letter to binary coded number */
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

int set_features(const cmd_args *args,
		       double ***features){
  char *seq_head, *seq;
  unsigned long seq_len, bin_num;

  {
    /* read fasta file */
    fasta_read(args->fasta_file, 
	       &seq_head, &seq, &seq_len);
    
    bin_num = (seq_len / args->res);
  
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "sequence: %s (len = %ld) => %ld bins)\n", 
	    seq_head, seq_len, bin_num);
  }

  {
    /* allocate memory for k-mer frequency odds table */  
    *features = calloc_errchk(bin_num, sizeof(double *),
				    "features");			      
  }

  {
    const int k = args->k;
    const int res = args->res;
    const int margin = args->margin;
    const unsigned long bin_min = (long)((margin + res - 1) / res);      
    const unsigned long bin_max = (long)((seq_len - margin - k + 1) / res);
    const unsigned int bit_mask = (1 << (2 * k)) - 1;
    unsigned long bin, valid_bin_num = 0;
    unsigned int *kmer_freq_sum = calloc_errchk(bit_mask + 1,
						sizeof(unsigned int),
						"calloc kmer_freq_sum");

    /* count k-mer frequency */
    for(bin = bin_min; bin < bin_max; bin++){
      unsigned int contain_n = 0, i;
      for(i = bin * res - margin;
	  i < (bin + 1) * res + k - 1 + margin; i++){
	if(seq[i] == 'N' || seq[i] == 'n'){	
	  contain_n = 1;
	  break;
	}
      }
      if(contain_n != 0){
	(*features)[bin] = NULL;
      }else{       
	unsigned int kmer = 0;
	/* For bins not containing 'N', allocate memory */
	(*features)[bin] = calloc_errchk(bit_mask + 1,
					       sizeof(double),
					       "calloc (*features)[]");
	/* convert first (k-1)-mer to bit-encoded sequence */
	for(i = bin * res - margin;
	    i < bin * res - margin + k - 1; i++){
	  kmer <<= 2;
	  kmer += (c2i(seq[i]) & 3);
	}
	/* count k-mer frequency */
	for(;
	    i < (bin + 1) * res + k - 1 + margin; i++){
	  kmer <<= 2;
	  kmer += (c2i(seq[i]) & 3);
	  (*features)[bin][(kmer & bit_mask)] += 1.0;
	  kmer_freq_sum[(kmer & bit_mask)] += 1;
	}
	valid_bin_num++;
      }    
    }

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "# of valid bins : %ld\n", valid_bin_num);

    /* conver to k-mer frequency odds */
    for(bin = bin_min; bin < bin_max; bin++){
      if((*features)[bin] != NULL){
	unsigned int kmer;
	double sum = 0;
	for(kmer = 0; kmer < bit_mask + 1; kmer++){
#if 0
	  (*features)[bin][kmer] /= (kmer_freq_sum[kmer]);
#endif
	  sum += (*features)[bin][kmer] * (*features)[bin][kmer];
	}

	/* normalize L_1 norm */
	for(kmer = 0; kmer < bit_mask + 1; kmer++){
	  (*features)[bin][kmer] /= (bit_mask + 1);
	}

#if 0
	/* normalize so that ||(*feature)[bin]||^2 = 1 */
	for(kmer = 0; kmer < bit_mask + 1; kmer++){
	  (*features)[bin][kmer] /= sum;
	}
#endif
      }
    }
  }

  fprintf(stderr, "%s [INFO] ", args->prog_name);
  fprintf(stderr, "computation of k-mer frequency odds finished\n");

  return 0;
}


#endif
