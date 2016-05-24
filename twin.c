#include <stdio.h>

#include "src/constant.h"
#include "src/cmd_args.h"
#include "src/fasta.h"

int main(int argc, char **argv){  
  cmd_args *args;

  cmd_args_parse(argc, argv, &args);
  cmd_args_chk(args);

  unsigned int **kmer_freq;
  set_kmer_freq((const cmd_args *)args, &kmer_freq);


  return 0;
}
