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
  kmer *kmers;
  {
    set_features((const cmd_args *)args, &features);
    hic_read((const cmd_args *)args, &data);
    kmer_read((const cmd_args *)args, &kmers);
  }

#if 0
  l2boost *model;
  {
    FILE *fp_out;
    if((fp_out = fopen(args->out_file, "w")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      args->out_file, strerror(errno));
      exit(EXIT_FAILURE);
    }
    l2boost_step_dump_head(fp_out);
    fprintf(fp_out, "%d\t%ld\t%e\t%e\t%f\t%f\n",
	    0, (long int)0, 0.0, 1.0, 0.0, 0.0);

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
		  (const double **)features,
		  (const hic *)data,
		  (const canonical_kp *)ckps,
		  (const double)args->acc,
		  &model,
		  fp_out);

    fclose(fp_out);

  }
#endif
  return 0;
}
