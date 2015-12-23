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

#include "calloc_errchk.h"

typedef struct _command_line_arguements {
  /* parameters */
  int chr;
  unsigned int k;
  unsigned int res;
  unsigned int min_size;
  unsigned int max_size;
  unsigned long iteration_num;
  double percentile;
  char *norm;
  char *exp;
  /* input */
  char *fasta_file;
  char *hicRaw_file;
  char *kmerFreq_file;
  char *hic_file;
  char *boost_oracle_file;
  /* output */
  char *output_dir;
  /* exec_mode */
  int exec_mode_quite;
  int exec_mode_skip_prep;
} command_line_arguements;


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
  return 0;
}

int check_args(const command_line_arguements *args,
	       const char *prog_name){
  int errflag = 0;

  if(args->k <= 0){	       
    show_error(stderr, prog_name, "k is not specified");
    errflag++;
  }else if(errflag == 0){
    printf("  k: %d\n", args->k);
  }

  if(args->res <= 0){
    show_error(stderr, prog_name, "resolution is not specified");
    errflag++;
  }else if(errflag == 0){
    printf("  resolution: %d\n", args->res);
  }

  if(args->min_size <= 0){
    show_error(stderr, prog_name, "minimum size is not specified");
    errflag++;
  }else if(errflag == 0){
    printf("  minimum size: %d\n", args->min_size);
  }

  if(args->max_size <= 0){
    show_error(stderr, prog_name, "maximum size is not specified");
    errflag++;
  }else if(errflag == 0){
    printf("  maximum size: %d\n", args->max_size);
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
    printf("  num. of weak lerners: %ld\n", args->iteration_num);
  }

  if(args->percentile <= 0){
    show_error(stderr, prog_name, "percentile threshold is not specified");
    errflag++;
  }else if(errflag == 0){
    printf("  percentile threshold: %e\n", args->percentile);
  }

  if(args->output_dir == NULL){
    show_warning(stderr, prog_name, "output directory is not specified");
    show_warning(stderr, prog_name, "results will be written to stdout");
  }else if(errflag == 0){
    printf("  Output dir:    %s\n", args->output_dir);
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
    {"quite",         no_argument, NULL, 'q'},
    {"skipPrep",      no_argument, NULL, 's'},
    {0, 0, 0, 0}
  };

  args = calloc_errchk(1, sizeof(command_line_arguements), 
		       "calloc: command line args");

  while((opt = getopt_long(argc, argv, "hvc:r:k:m:M:i:p:n:e:g:R:f:H:O:o:qs",
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
	args->hicRaw_file = optarg;
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
	break;
      /* exec_mode */
      case 'q': /* quite */
	args->exec_mode_quite = 1;
	break;
      case 's': /* skipPrep */
	args->exec_mode_skip_prep = 1;
	break;
    }
  }

  check_args(args, argv[0]);

  main_sub(args, argv[0]);

  free(args);
  return 0;
}

