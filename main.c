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
#include "hic.h"
#include "fasta.h"


int show_usage(FILE *fp, const char *prog_name){
  fprintf(fp, "%s:Info: Usage\n", prog_name);
  return 0;
}

int show_error(FILE *fp,
		const char *prog_name,
		const char *errmsg){
  fprintf(fp, "%s: error: %s\n", prog_name, errmsg);
  return 0;
}

int show_warning(FILE *fp,
		 const char *prog_name,
		 const char *warnmsg){
  fprintf(fp, "%s: warning: %s\n", prog_name, warnmsg);
  return 0;
}


int main_sub(const command_line_arguements *args,
	     const char *prog_name){

#if 0
  hic *hic;
  hic_prep(args, &hic);
#endif

  unsigned int **kmer_freq;
  set_kmer_freq(args, &kmer_freq);

  return 0;
}



int check_args(const command_line_arguements *args,
	       const char *prog_name){
  int errflag = 0;

  if(args->k <= 0){	       
    show_error(stderr, prog_name, "k is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: k: %d\n", prog_name, args->k);
  }

  if(args->res <= 0){
    show_error(stderr, prog_name, "resolution is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: resolution: %d\n", prog_name, args->res);
  }

  if(args->min_size <= 0){
    show_error(stderr, prog_name, "minimum size is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: minimum size: %d\n", prog_name, args->min_size);
  }

  if(args->max_size <= 0){
    show_error(stderr, prog_name, "maximum size is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: maximum size: %d\n", prog_name, args->max_size);
  }

  if(args->iteration_num <= 0){
    show_error(stderr, prog_name, "iteration number is not specified");
    errflag++;
  }else if((unsigned long)(1 << (4 * args->k)) < args->iteration_num){
    show_error(stderr, prog_name, "iteration number is invalid");
    fprintf(stderr, "number of weak lerners(T = %ld) exceeds 16^k (%d)\n",
	    args->iteration_num, 1 << (4 * args->k));
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: number of iterations in AdaBoost: %ld\n",
	    prog_name, args->iteration_num);
  }

  if(args->percentile <= 0){
    show_error(stderr, prog_name, "percentile threshold is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: percentile threshold: %e\n",
	    prog_name, args->percentile);
  }

  if(args->fasta_file != NULL){
    fprintf(stderr, "%s: Info: fasta file: %s\n",
	    prog_name, args->fasta_file);
  }

  if(args->hicRaw_dir != NULL){
    fprintf(stderr, "%s: Info: Hi-C raw data directory: %s/\n", 
	    prog_name, args->hicRaw_dir);
  }

  if(args->kmerFreq_file != NULL){
    fprintf(stderr, "%s: Info: k-mer frequency count file: %s\n", 
	    prog_name, args->kmerFreq_file);
  }

  if(args->hic_file != NULL){
    fprintf(stderr, "%s: Info: pre-processed Hi-C file: %s\n",
	    prog_name, args->hic_file);
  }

  if(args->boost_oracle_file != NULL){
    fprintf(stderr, "%s: Info: oracle file for AdaBoost: %s\n",
	    prog_name, args->boost_oracle_file);
  }

  if(args->output_dir == NULL){
    show_warning(stderr, prog_name, "output directory is not specified");
    show_warning(stderr, prog_name, "results will be written to stdout");
  }else if(errflag == 0){
    fprintf(stderr, "%s: Info: output dir: %s/\n", 
	    prog_name, args->output_dir);
  }

  if(args->exec_thread_num > 0){	       
    fprintf(stderr, "%s: Info: thread num: %d\n", 
	    prog_name, args->exec_thread_num);
  }

  if(errflag > 0){
    show_usage(stderr, prog_name);
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
	fprintf(stdout, "version: 0.20\n");
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

  check_args(args, argv[0]);

  /* set exec_thread_num */
  if(args->exec_thread_num <= 0){
    args->exec_thread_num = (int)sysconf(_SC_NPROCESSORS_ONLN);
  }

  main_sub(args, argv[0]);

  free(args);
  return 0;
}


