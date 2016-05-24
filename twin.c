#include <stdio.h>

#include "src/constant.h"
#include "src/cmd_args.h"
#include "src/fasta.h"
#include "src/hic.h"

int main(int argc, char **argv){  
  cmd_args *args;
  {
    cmd_args_parse(argc, argv, &args);
    cmd_args_chk(args);
  }

  double **kmer_freq_odds;
  set_kmer_freq_odds((const cmd_args *)args, &kmer_freq_odds);

  hic *data;
  hic_read((const cmd_args *)args, &data);

  return 0;
}
