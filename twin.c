#include <stdio.h>

#include "src/constant.h"
#include "src/cmd_args.h"
#include "src/fasta.h"
#include "src/hic.h"
#include "src/kmer.h"
#include "src/l2boost.h"

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

  l2boost *model;
  {
    FILE *fp_out;
    if((fp_out = fopen(args->out_file, "w")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      args->out_file, strerror(errno));
      exit(EXIT_FAILURE);
    }
    l2boost_step_dump_head(fp_out);
    l2boost_step_dump((const l2boost *)model, 0, 0, 0, 0, 0, fp_out);

    if(args->pri_file != NULL){
      l2boost_load((const cmd_args *)args,
		   (const canonical_kp *)ckps,	       
		   (const unsigned int)args->iter1,
		   (const char *)args->pri_file,
		   &model,
		   fp_out);      
    }else{
      l2boost_init((const canonical_kp *)ckps,	       
		   (const unsigned int)args->iter1,
		   &model);
    }

    l2boost_train((const cmd_args *)args,		
		  (const double **)kmer_freq_odds,
		  (const hic *)data,
		  (const canonical_kp *)ckps,
		  (const double)args->acc,
		  &model,
		  fp_out);

    fclose(fp_out);

  }
  return 0;
}
