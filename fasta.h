#ifndef __FASTA_H__
#define __FASTA_H__

#include "constant.h"
#include "cmd_args.h"
#include "mywc.h"
#include "calloc_errchk.h"
#include "diffSec.h"
#include "io.h"

/**
 * This header file contains some functions to perform the following tasks
 * - read FASTA format file and store seqeunce to memory 
 * - compute k-mer frequencies for bins
 */


/* arguements for function */
typedef struct _kmer_freq_count_args{
  int thread_id;
  unsigned long begin;
  unsigned long end;
  unsigned int k;
  unsigned int res;
  char *seq;
  unsigned int ***kmer_freq;
} kmer_freq_count_args;


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

#if 1
int set_kmer_freq(const command_line_arguements *cmd_args,
		  unsigned int ***kmer_freq){
  const unsigned int bit_mask = (1 << (2 * (cmd_args->k))) - 1;
  char *seq_head, *seq;
  unsigned long seq_len, bin_num, bin;
  unsigned int i, contain_n, kmer;
  int *interval;
  unsigned long interval_num;
  unsigned long interval_i = 0;

  /* read fasta file */
  fasta_read(cmd_args->fasta_file, 
	     &seq_head, &seq, &seq_len);

  bin_num = (seq_len / cmd_args->res);

  fprintf(stderr, "%s: info: sequence: %s (%ld : %ld)\n", 
	  cmd_args->prog_name, seq_head, seq_len, bin_num);


  /* allocate memory for k-mer frequency table */  
  *kmer_freq = calloc_errchk(bin_num, sizeof(unsigned int *),
			     "kmer_freq");			      

  /* get genomic intervals of interest */
  read_int(cmd_args->interval_file,
	   &interval, &interval_num);

  /* count k-mer frequency */
  for(bin = 1; bin < bin_num - 1; bin++){
    /* [1, bin_num - 1) is because we check adjascent bins to calc. k-mer freq */

    if((interval_i >= interval_num) ||
       (((interval_i % 2) == 0) && (((int)bin) < interval[interval_i]))){
      /*
       * *kmer_freq_vector will be NULL if either of the following holds
       * 1] We already checked all the intervals of interest
       * 2] We are waiting for a next interval ( interval_i % 2 == 0) and
       *    this bin is still not in the next interval (bin < interval[interval_i])
       */
      (*kmer_freq)[bin] = NULL;

    }else{           
      if(((int)bin) == interval[interval_i]){
	/*
	 * If this bin is at the boundary of the current interval, 
	 * increment the index variable interval_i
	 */
	interval_i++;
      }

      contain_n = 0;
      for(i = (bin - 1)* cmd_args->res; 
	  i < (bin + 2) * cmd_args->res + cmd_args->k - 1; i++){
	if(seq[i] == 'N' || seq[i] == 'n'){	
	  contain_n = 1;
	  break;
	}
      }
      if(contain_n != 0){
	(*kmer_freq)[bin] = NULL;
      }else{
	
	/* For bins not containing 'N', allocate memory */
	(*kmer_freq)[bin] = calloc_errchk((1 << (2 * (cmd_args->k))),
					  sizeof(unsigned int),
					  "calloc kmer_freq[]");
	kmer = 0;
	/* convert first (k-1)-mer to bit-encoded sequence */
	for(i = (bin - 1) * cmd_args->res;
	    i < (bin - 1)* cmd_args->res + cmd_args->k - 1; i++){
	  kmer <<= 2;
	  kmer += (c2i(seq[i]) & 3);
	}
	/* count k-mer frequency */
	for(;
	    i < (bin + 2) * cmd_args->res + cmd_args->k - 1; i++){
	  kmer <<= 2;
	  kmer += (c2i(seq[i]) & 3);
	  (*kmer_freq)[bin][(kmer & bit_mask)] += 1;
	}
      }    
    }
  }

  return 0;
}
#endif

#if 0
void *kmer_freq_count(void *args){
  kmer_freq_count_args *params = (kmer_freq_count_args *)args;
  unsigned int bin = 0, i = 0, kmer = 0, contain_n = 0;
  const unsigned int bit_mask = (1 << (2 * (params->k))) - 1;

#if 1
  fprintf(stderr, "thread %d: start [%ld : %ld]\n",
	  params->thread_id, params->begin, params->end);  
#endif

  for(bin = params->begin; bin <= params->end; bin++){
    contain_n = 0;

    /* check if this genome bin contains 'N' and allocate memory */
    for(i = bin * params->res;
	i < (bin + 1) * params->res + params->k - 1; i++){
      if((params->seq)[i] == 'N' || (params->seq)[i] == 'n'){	
	contain_n = 1;
	break;
      }
    }
  
    if(contain_n != 0){

      /* if the bin contains 'N', then skip */
      (*(params->kmer_freq))[bin] = NULL;

    }else{

      /* For bins not containing 'N', allocate memory */
      (*(params->kmer_freq))[bin] = calloc_errchk((1 << (2 * (params->k))),
						  sizeof(unsigned int),
						  "calloc kmer_freq[]");
      kmer = 0;
      /* convert first (k-1)-mer to bit-encoded sequence */
      for(i = bin * params->res; 
	  i < bin * params->res + params->k - 1; i++){
	kmer <<= 2;
	kmer += c2i((params->seq)[i]);
      }
      /* count k-mer frequency */
      for(i = bin * params->res;
	  i < (bin + 1) * params->res + params->k - 1; i++){
	kmer <<= 2;
	kmer += c2i((params->seq)[i]);
	(*(params->kmer_freq))[bin][(kmer & bit_mask)] += 1;
      }
    }
  }
  return NULL;
}

int set_kmer_freq_pthread(const command_line_arguements *cmd_args,
		  unsigned int ***kmer_freq){
  char *seq_head, *seq;
  unsigned long seq_len, bin_num;
 
  /* read fasta file */
  fasta_read(cmd_args->fasta_file, 
	     &seq_head, &seq, &seq_len);
  fprintf(stderr, "%s: info: sequence: %s (%ld)\n", 
	  cmd_args->prog_name, seq_head, seq_len);
  
  /* count k-mer frequency */
  {
    kmer_freq_count_args *params;
    pthread_t *threads = NULL;
    int i;
    bin_num = (seq_len / cmd_args->res); 

    /* allocate memory for k-mer frequency table */  
    *kmer_freq = calloc_errchk(bin_num, sizeof(unsigned int *),
			       "kmer_freq");			      

    /* allocate memory for threads */
    params = calloc_errchk(cmd_args->exec_thread_num,			   
			   sizeof(kmer_freq_count_args),
			   "calloc: hic_prep_thread_args");
    threads = calloc_errchk(cmd_args->exec_thread_num,			   
			    sizeof(pthread_t),
			    "calloc: threads");
    /* set variables */
    for(i = 0; i < cmd_args->exec_thread_num; i++){
      params[i].thread_id = i;
      params[i].begin = ((i == 0) ? 0 : params[i - 1].end + 1);
      params[i].end = ((i == (cmd_args->exec_thread_num - 1)) ?
		       bin_num - 1 :
		       ((bin_num / cmd_args->exec_thread_num) * (i + 1) - 1));
      params[i].k = cmd_args-> k;
      params[i].res = cmd_args->res;
      params[i].seq = seq;
      params[i].kmer_freq = kmer_freq;
    }

    /* pthread create */
    for(i = 0; i < cmd_args->exec_thread_num; i++){
      pthread_create(&threads[i], NULL, 
		     kmer_freq_count,
		     (void*)&params[i]);
    }

    /* pthread join */
    for(i = 0; i < cmd_args->exec_thread_num; i++){
      pthread_join(threads[i], NULL);
    }
  }

  {
    unsigned int i, j;
    for(i = 0; i < bin_num; i++){
      fprintf(stderr, "%10d: ", i);
      if(kmer_freq[i] != NULL){
	for(j = 0; j < (1 << (2 * (cmd_args->k))); j++){
	  fprintf(stderr, "%d", kmer_freq[i][j] > 0 ? 1 : 0);
	}
	fprintf(stderr, "\n");
      }else{
	fprintf(stderr, "*\n");
      }
    }
  }


  free(seq);
  free(seq_head);

  show_info(stderr, cmd_args->prog_name, 
	    "Complete k-mer frequency table computation");

  return 0;
}
#endif
#endif
