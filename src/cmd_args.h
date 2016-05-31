#ifndef __CMD_ARGS__
#define __CMD_ARGS__

#include <string.h>
#include <getopt.h>
#include <unistd.h>

#include "calloc_errchk.h"

	      
typedef struct _cmd_args {
  /* parameters */
  int k;
  int res;
  int margin;
  int iter1;
  int iter2;
  double acc;
  /* input */
  char *fasta_file;
  char *hic_file;
  char *kmer_pair;
  /* output */
  char *out_file;
  /* saved results */
  char *pri_file;
  char *sec_file;
  /* exec_mode */
  int verbose_level;
  int thread_num;
  char *prog_name;
} cmd_args;


int show_usage(FILE *, const char *);
int cmd_args_chk(const cmd_args *);
int cmd_args_chk_pred(const cmd_args *);
int cmd_args_parse(const int, char **, cmd_args **);
		   

int show_usage(FILE *fp, 
	       const char *prog_name){
  fprintf(fp, "%s [INFO] ", prog_name);
  fprintf(fp, "usage:\n");
  fprintf(fp, 
	  "%s -k k --res r [--margin M] --iter1 n --iter2 m --acc a --fasta f --hic H --kmer c --out o [--pri p] [--sec s] [--verbose V] --thread_num t \n",
	  prog_name);
  return 0;
}

int show_usage_pred(FILE *fp, 
		    const char *prog_name){
  fprintf(fp, "%s [INFO] ", prog_name);
  fprintf(fp, "usage:\n");
  fprintf(fp, 
	  "%s -k k --res r [--margin M] --fasta f --hic H --kmer c --out o --pri p [--verbose V] --thread_num t \n",
	  prog_name);
  return 0;
}

int cmd_args_chk(const cmd_args *args){
  int errflag = 0;

  /* parameters */

  if(args->k <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "k is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "k", args->k);
  }

  if(args->res <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "res is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "res", args->res);
  }

  if(args->margin < 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "margin is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "margin", args->margin);
  }

  if(args->iter1 <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "iter1 is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "iter1", args->iter1);
  }

  if(args->iter2 <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "iter2 is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "iter2", args->iter2);
  }

  if(args->acc <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "acc is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %f\n", "acc", args->acc);
  }

  /* input */

  if(args->fasta_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "fasta file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "fasta_file", args->fasta_file);
  }

  if(args->hic_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "hic file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "hic_file", args->hic_file);
  }

  if(args->kmer_pair == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "k-mer pair file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "kmer_pair", args->kmer_pair);
  }

  /* output */

  if(args->out_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "output file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s.{pri, sec}\n", "out_file", args->out_file);
  }

  /* saved results */

  if(args->pri_file != NULL && errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "pri_file", args->pri_file);
  }

  if(args->sec_file != NULL && errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "sec_file", args->sec_file);
  }

  /* exec mode */
  
  if(args->thread_num <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "thread_num is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "thread_num", args->thread_num);
  }

  if(errflag > 0){
    show_usage(stderr, args->prog_name);
    exit(EXIT_FAILURE);
  }


  return 0;
}

int cmd_args_chk_pred(const cmd_args *args){
  int errflag = 0;

  /* parameters */

  if(args->k <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "k is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "k", args->k);
  }

  if(args->res <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "res is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "res", args->res);
  }

  if(args->margin < 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "margin is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "margin", args->margin);
  }

  /* input */

  if(args->fasta_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "fasta file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "fasta_file", args->fasta_file);
  }

  if(args->hic_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s", "hic file is not specified");
    fprintf(stderr, "%s\n", " (we require Hi-C file to obtain the target coordinates)");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "hic_file", args->hic_file);
  }

  if(args->kmer_pair == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "k-mer pair file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "kmer_pair", args->kmer_pair);
  }

  /* output */

  if(args->out_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "output file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s.cmp\n", "out_file", args->out_file);
  }

  /* saved results */

  if(args->pri_file == NULL){
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "pri file is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %s\n", "pri_file", args->pri_file);
  }

  /* exec mode */
  
  if(args->thread_num <= 0){	       
    fprintf(stderr, "%s [ERROR] ", args->prog_name);
    fprintf(stderr, "%s\n", "thread_num is not specified");
    errflag++;
  }else if(errflag == 0){
    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "%s : %d\n", "thread_num", args->thread_num);
  }

  if(errflag > 0){
    show_usage_pred(stderr, args->prog_name);
    exit(EXIT_FAILURE);
  }


  return 0;
}


int cmd_args_parse(const int argc, char **argv,	       
		   cmd_args **args){
  
  int opt = 0, opt_idx = 0;
  struct option long_opts[] = {
    {"help",    no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    /* parameters */
    {"k",         required_argument, NULL, 'k'},
    {"res",       required_argument, NULL, 'r'},
    {"margin",    required_argument, NULL, 'M'},
    {"iter1",     required_argument, NULL, 'n'},
    {"iter2",     required_argument, NULL, 'm'},
    {"acc",       required_argument, NULL, 'a'},
    /* input */
    {"fasta",     required_argument, NULL, 'f'},
    {"hic",       required_argument, NULL, 'H'},
    {"kmer",      required_argument, NULL, 'c'},
    /* output */
    {"out",       required_argument, NULL, 'o'},
    /* saved results*/
    {"pri",       required_argument, NULL, 'p'},
    {"sec",       required_argument, NULL, 's'},
    /* exec_mode */
    {"verbose",   required_argument, NULL, 'V'},
    {"thread",    required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  *args = calloc_errchk(1, sizeof(cmd_args), 
			"calloc: command line args");

  while((opt = getopt_long(argc, argv, "hvk:r:M:n:m:a:f:H:c:o:p:s:V:t:",
			   long_opts, &opt_idx)) != -1){
    switch (opt){
      case 'h': /* help */
	show_usage(stdout, argv[0]);
	exit(EXIT_SUCCESS);
      case 'v': /* version*/
	fprintf(stdout, "version: 0.56\n");
	exit(EXIT_SUCCESS);

      /* parameters */
      case 'k': /* k */
	(*args)->k = atoi(optarg);
	break;
      case 'r': /* res */
	(*args)->res = atoi(optarg);
	break;
      case 'M': /* margin */
	(*args)->margin = atoi(optarg);
	break;
      case 'n': /* iter1 */
	(*args)->iter1 = atoi(optarg);
	break;
      case 'm': /* iter2 */
	(*args)->iter2 = atoi(optarg);
	break;
      case 'a': /* acceleration */
	(*args)->acc = atof(optarg);
	break;

      /* input */
      case 'f': /* fasta */
	(*args)->fasta_file = optarg;
	break;
      case 'H': /* hic */
	(*args)->hic_file = optarg;
	break;
      case 'c': /* kmer */
	(*args)->kmer_pair = optarg;
	break;

      /* output */
      case 'o': /* out */
	(*args)->out_file = optarg;
	break;

      /* saved results */
      case 'p': /* pri */
	(*args)->pri_file = optarg;
	break;
      case 's': /* sec */
	(*args)->sec_file = optarg;
	break;

      /* exec_mode */
      case 'V': /* verbose */
	(*args)->verbose_level = atoi(optarg);
	break;
      case 't': /* thread_num */
	(*args)->thread_num = atoi(optarg);
	break;
    }
  }

  {
    (*args)->prog_name = calloc_errchk(strlen(argv[0]), sizeof(char),
				    "(*args)->prog_name");
    strncpy((*args)->prog_name, argv[0], strlen(argv[0]));
  }

  /* set thread_num */
  if((*args)->thread_num <= 0){
    (*args)->thread_num = (int)sysconf(_SC_NPROCESSORS_ONLN);
  }

  /* set margin */
  if((*args)->margin <= 0){
    (*args)->margin = 0;
  }

  return 0;

}


#endif
