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

  double **features;
  hic *data;
  canonical_kp *ckps;
  {
    set_features((const cmd_args *)args, &features);
    hic_read((const cmd_args *)args, &data);
    canonical_kp_read((const cmd_args *)args, &ckps);
  }

  boost *model;
  {
    FILE *fp_out;
    if((fp_out = fopen(args->out_file, "w")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      args->out_file, strerror(errno));
      exit(EXIT_FAILURE);
    }


    boost_init((const cmd_args *)args,
		 (const canonical_kp *)ckps,	       
		 (const unsigned int)args->iter1,
		 (const char *)args->pri_file,
		 &model,
		 fp_out);      

    l2_train((const cmd_args *)args,		
		  (const double **)features,
		  (const hic *)data,
		  (const canonical_kp *)ckps,
		  (const double)args->acc,
		  &model,
		  fp_out);

    fclose(fp_out);

  }
  return 0;
}
