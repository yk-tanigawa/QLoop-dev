#ifndef __CMD_ARGS__
#define __CMD_ARGS__


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
  char *hicRaw_dir;
  char *kmerFreq_file;
  char *hic_file;
  char *boost_oracle_file;
  /* output */
  char *output_dir;
  /* exec_mode */
  int exec_mode_quite;
  int exec_mode_skip_prep;
  int exec_mode_QP_only;
  int exec_thread_num;
  char *prog_name;
} command_line_arguements;


#endif
