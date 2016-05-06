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
#include <libgen.h>
#include <pthread.h>

#include "calloc_errchk.h"

#include "cmd_args.h"
#include "constant.h"
#include "filename.h"
#include "show_msg.h"
#include "hic.h"
#include "fasta.h"
#include "threshold.h"
#include "adaboost.h"
#include "qp.h"

int debug_dump_kmer_freq(const unsigned int **kmer_freq,
			 const unsigned int k,
			 const int len){
  int i, j;
  for(i = 0; i < len; i++){
    fprintf(stderr, "%d ", i);
    
    if(kmer_freq[i] == NULL){
      fprintf(stderr, "*\n");
    }else{
      for(j = 0; j < (1 << (2 * k)); j++){
	fprintf(stderr, "%d", kmer_freq[i][j] > 0 ? 1 : 0);
      }
      fprintf(stderr, "\n");
    }
  }
  return 0;
}



int main_sub(const command_line_arguements *args){

#if 0
  unsigned int **kmer_freq;

  set_kmer_freq(args, &kmer_freq);
  
  debug_dump_kmer_freq((const unsigned int **)kmer_freq, args->k, 100);
#endif

  unsigned int **kmer_freq;
  hic *hic;
  adaboost *model;
  canonical_kp *kp;
  thresholds *th;
  double **P, *q;
  filenames *fnames;

  set_filenames(args, &fnames);

  set_kmer_freq(args, &kmer_freq);
  hic_prep(args, &hic);
  hic_check_kmer(hic, (const unsigned int **)kmer_freq, args->prog_name);
  hic_pack(hic, args->prog_name);    


  set_canonical_kmer_pairs(args->k, &kp);
  set_thresholds(hic->mij, 1000, hic->nrow, &th);
  write_histo(args, th, fnames->histo);

  adaboost_learn(args,
		 (const unsigned int **)kmer_freq,
		 hic,
		 get_threshold(args, th, args->percentile),
		 kp,
		 &model,
		 fnames->adaboost);
  qp_prep(args,
	  (const unsigned int **)kmer_freq,
	  hic,
	  kp,
	  model,
	  &P, &q, 
	  get_threshold(args, th, 0.005),
	  get_threshold(args, th, 0.995),
	  fnames->qp_P,
	  fnames->qp_q);
  return 0;
}



int check_args(const command_line_arguements *args){
  int errflag = 0;

  if(args->k <= 0){	       
    show_error(stderr, args->prog_name, "k is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: k: %d\n", args->prog_name, args->k);
  }

  if(args->res <= 0){
    show_error(stderr, args->prog_name, "resolution is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: resolution: %d\n", args->prog_name, args->res);
  }

  if(args->min_size <= 0){
    show_error(stderr, args->prog_name, "minimum size is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: minimum size: %d\n", args->prog_name, args->min_size);
  }

  if(args->max_size <= 0){
    show_error(stderr, args->prog_name, "maximum size is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: maximum size: %d\n", args->prog_name, args->max_size);
  }

  if(args->iteration_num <= 0){
    show_error(stderr, args->prog_name, "iteration number is not specified");
    errflag++;
  }else if((unsigned long)(1 << (4 * args->k)) < args->iteration_num){
    show_error(stderr, args->prog_name, "iteration number is invalid");
    fprintf(stderr, "number of weak lerners(T = %ld) exceeds 16^k (%d)\n",
	    args->iteration_num, 1 << (4 * args->k));
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: number of iterations in AdaBoost: %ld\n",
	    args->prog_name, args->iteration_num);
  }

  if(args->percentile <= 0){
    show_error(stderr, args->prog_name, "percentile threshold is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: percentile threshold: %e\n",
	    args->prog_name, args->percentile);
  }

  if(args->fasta_file != NULL){
    fprintf(stderr, "%s: info: fasta file: %s\n",
	    args->prog_name, args->fasta_file);
  }

  if(args->hicRaw_dir != NULL){
    fprintf(stderr, "%s: info: Hi-C raw data directory: %s/\n", 
	    args->prog_name, args->hicRaw_dir);
  }

  if(args->kmerFreq_file != NULL){
    fprintf(stderr, "%s: info: k-mer frequency count file: %s\n", 
	    args->prog_name, args->kmerFreq_file);
  }

  if(args->hic_file != NULL){
    fprintf(stderr, "%s: info: pre-processed Hi-C file: %s\n",
	    args->prog_name, args->hic_file);
  }

  if(args->boost_oracle_file != NULL){
    fprintf(stderr, "%s: info: oracle file for AdaBoost: %s\n",
	    args->prog_name, args->boost_oracle_file);
  }

  if(args->interval_file != NULL){
    fprintf(stderr, "%s: info: interval file to specify the region of interest: %s\n",
	    args->prog_name, args->interval_file);
  }

  if(args->output_dir == NULL){
    show_warning(stderr, args->prog_name, "output directory is not specified");
    show_warning(stderr, args->prog_name, "results will be written to stdout");
  }else if(errflag == 0){
    fprintf(stderr, "%s: info: output dir: %s/\n", 
	    args->prog_name, args->output_dir);
  }

  if(args->exec_thread_num > 0){	       
    fprintf(stderr, "%s: info: thread num: %d\n", 
	    args->prog_name, args->exec_thread_num);
  }

  if(errflag > 0){
    show_usage(stderr, args->prog_name);
    exit(EXIT_FAILURE);
  }
  return 0;

}

int main(int argc, char **argv){  
  command_line_arguements *args;
  int opt = 0, opt_idx = 0;
  struct option long_opts[] = {
    {"help",          no_argument, NULL, 'h'},
    {"version",       no_argument, NULL, 'v'},
    /* parameters */
    {"chr",           required_argument, NULL, 'c'},
    {"k",             required_argument, NULL, 'k'},
    {"res",           required_argument, NULL, 'r'},
    {"min_size",      required_argument, NULL, 'm'},
    {"max_size",      required_argument, NULL, 'M'},
    {"iteration_num", required_argument, NULL, 'i'},
    {"percentile",    required_argument, NULL, 'p'},
    {"norm",          required_argument, NULL, 'n'},
    {"expected",      required_argument, NULL, 'e'},
    /* input */
    {"fasta",         required_argument, NULL, 'g'},
    {"hicRaw",        required_argument, NULL, 'R'},
    {"kmerFreq",      required_argument, NULL, 'f'},
    {"hic",           required_argument, NULL, 'H'},
    {"boostOracle",   required_argument, NULL, 'O'},
    {"interval",      required_argument, NULL, 'I'},
    /* output */
    {"out",           required_argument, NULL, 'o'},
    /* exec_mode */
    {"quite",         no_argument,       NULL, 'q'},
    {"skipPrep",      no_argument,       NULL, 's'},
    {"QPonly",        no_argument,       NULL, 'Q'},
    {"thread_num",    required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  args = calloc_errchk(1, sizeof(command_line_arguements), 
		       "calloc: command line args");

  while((opt = getopt_long(argc, argv, "hvc:r:k:m:M:i:p:n:e:g:R:f:H:O:o:qsQt:",
			   long_opts, &opt_idx)) != -1){
    switch (opt){
      case 'h': /* help */
	show_usage(stdout, argv[0]);
	exit(EXIT_SUCCESS);
      case 'v': /* version*/
	fprintf(stdout, "version: 0.41\n");
	exit(EXIT_SUCCESS);
      /* parameters */
      case 'c': /* chr */
	args->chr = atoi(optarg);
	break;
      case 'k': /* k */
	args->k = atoi(optarg);
	break;
      case 'r': /* res */
	args->res = atoi(optarg);
	break;
      case 'm': /* min_size */
	args->min_size = atoi(optarg);
	break;
      case 'M': /* max_size */
	args->max_size = atoi(optarg);
	break;
      case 'i': /* iteration_num */
	args->iteration_num = atol(optarg);
	break;
      case 'p': /* percentile */
	args->percentile = atof(optarg);
	break;
      case 'n': /* norm */
	args->norm = optarg;
	break;
      case 'e': /* expected */
	args->exp = optarg;
	break;
      /* input */
      case 'g': /* fasta */
	args->fasta_file = optarg;
	break;
      case 'R': /* hicRaw */
	args->hicRaw_dir = optarg;
	if((args->hicRaw_dir)[strlen(args->hicRaw_dir) - 1] == '/'){
	  (args->hicRaw_dir)[strlen(args->hicRaw_dir) - 1] = '\0';
	}
	break;
      case 'f': /* kmerFreq */
	args->kmerFreq_file = optarg;
	break;
      case 'H': /* hic */
	args->hic_file = optarg;
	break;
      case 'O': /* boostOracle */
	args->boost_oracle_file = optarg;
	break;
      case 'I': /* interval */
	args->interval_file = optarg;
	break;
      /* output */
      case 'o': /* out */
	args->output_dir = optarg;
	if((args->output_dir)[strlen(args->output_dir) - 1] == '/'){
	  (args->output_dir)[strlen(args->output_dir) - 1] = '\0';
	}
	break;
      /* exec_mode */
      case 'q': /* quite */
	args->exec_mode_quite = 1;
	break;
      case 's': /* skipPrep */
	args->exec_mode_skip_prep = 1;
	break;
      case 'Q': /* QPonly */
	args->exec_mode_QP_only = 1;
	break;
      case 't': /* thread_num */
	args->exec_thread_num = atoi(optarg);
	break;
    }
  }

  {
    args->prog_name = calloc_errchk(strlen(argv[0]), sizeof(char),
				    "args->prog_name");
    strncpy(args->prog_name, argv[0], strlen(argv[0]));
  }

  /* set exec_thread_num */
  if(args->exec_thread_num <= 0){
    args->exec_thread_num = (int)sysconf(_SC_NPROCESSORS_ONLN);
  }

  check_args(args);
  main_sub(args);

  free(args);
  return 0;
}


