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

int set_kmer_freq_odds(const cmd_args *args,
		       double ***kmer_freq_odds){
  unsigned int **kmer_freq;
  const unsigned int bit_mask = (1 << (2 * (args->k))) - 1;
  char *seq_head, *seq;
  unsigned long seq_len, bin_num, bin;
  unsigned int i, contain_n, kmer;

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
    /* allocate memory for k-mer frequency table */  
    kmer_freq = calloc_errchk(bin_num, sizeof(unsigned int *),
			      "kmer_freq");			      
    *kmer_freq_odds = calloc_errchk(bin_num, sizeof(double *),
				    "kmer_freq_odds");			      
  }

  {
    /* count k-mer frequency */
    for(bin = (int)((args->margin) / (args->res)) + 1;
	bin < bin_num - (int)((args->margin) / (args->res)) - 1; bin++){
      /* [1, bin_num - 1) is because we check adjascent bins to calc. k-mer freq */

      contain_n = 0;
      for(i = bin * (args->res) - (args->margin); 
	  i < (bin + 1) * (args->res) + (args->k) - 1 + (args->margin); i++){
	if(seq[i] == 'N' || seq[i] == 'n'){	
	  contain_n = 1;
	  break;
	}
      }
      if(contain_n != 0){
	kmer_freq[bin] = NULL;
      }else{       
	/* For bins not containing 'N', allocate memory */
	kmer_freq[bin] = calloc_errchk(bit_mask + 1,
				       sizeof(unsigned int),
				       "calloc kmer_freq[]");
	kmer = 0;
	
	/* convert first (k-1)-mer to bit-encoded sequence */
	for(i = bin * (args->res) - (args->margin);
	    i < bin * (args->res) - (args->margin) + (args->k) - 1; i++){
	  kmer <<= 2;
	  kmer += (c2i(seq[i]) & 3);
	}
	/* count k-mer frequency */
	for(;
	    i < (bin + 1) * (args->res) + (args->k) - 1 + (args->margin); i++){
	  kmer <<= 2;
	  kmer += (c2i(seq[i]) & 3);
	  kmer_freq[bin][(kmer & bit_mask)] += 1;
	}
      }    
    }
  }

  { /* conver k-mer frequency to k-mer frequency odds */

    unsigned int *kmer_freq_sum = calloc_errchk(bit_mask + 1,
						sizeof(unsigned int),
						"calloc kmer_freq_sum");
    int valid_bin_num = 0;
    unsigned int kmer = 0;

    /* count the sum of k-mer frequency */
    for(bin = (int)((args->margin) / (args->res)) + 1;
	bin < bin_num - (int)((args->margin) / (args->res)) - 1; bin++){
      if(kmer_freq[bin] != NULL){
	valid_bin_num ++;
	for(kmer = 0; kmer < bit_mask + 1; kmer++){
	  kmer_freq_sum[kmer] += kmer_freq[bin][kmer];
	}
      }
    }

    for(bin = (int)((args->margin) / (args->res)) + 1;
	bin < bin_num - (int)((args->margin) / (args->res)) - 1; bin++){
      if(kmer_freq[bin] != NULL){
	(*kmer_freq_odds)[bin] = calloc_errchk(bit_mask + 1,
					       sizeof(double),
					       "calloc kmer_freq_odds[]");
	for(kmer = 0; kmer < bit_mask + 1; kmer++){
	  (*kmer_freq_odds)[bin][kmer] = (1.0 * valid_bin_num * kmer_freq[bin][kmer]) / ((args->res) + 2 * (args->margin) * kmer_freq_sum[kmer]);	  
	}	
	
      }
      free(kmer_freq[bin]);
    }
    free(kmer_freq);
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "# of valid bins : %d\n", valid_bin_num);
  }
  fprintf(stderr, "%s [INFO] ", args->prog_name);
  fprintf(stderr, "computation of k-mer frequency odds finished\n");

  return 0;
}



#endif
