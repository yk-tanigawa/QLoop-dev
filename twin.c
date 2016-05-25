#include <stdio.h>

#include "src/constant.h"
#include "src/cmd_args.h"
#include "src/fasta.h"
#include "src/hic.h"
#include "src/kmer.h"

int main(int argc, char **argv){  
  cmd_args *args;
  {
    cmd_args_parse(argc, argv, &args);
    cmd_args_chk(args);
  }

  double **kmer_freq_odds;
  hic *data;
  canonical_kp *ckps;
  {
    set_kmer_freq_odds((const cmd_args *)args, &kmer_freq_odds);
    hic_read((const cmd_args *)args, &data);
    canonical_kp_read((const cmd_args *)args, &ckps);
  }
  return 0;
}
