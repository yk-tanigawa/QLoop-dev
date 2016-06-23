#include <stdio.h>

#include "src/constant.h"
#include "src/mywc.h"
#include "src/cmd_args.h"
#include "src/fasta.h"
#include "src/hic.h"
#include "src/kmer.h"
#include "src/l2boost.h"
#include "src/pred.h"

int main(int argc, char **argv){  
  cmd_args *args;
  {
    cmd_args_parse(argc, argv, &args);
    cmd_args_chk_pred(args);
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
  double *pred;
  {

    boost_init((const cmd_args *)args,
	       (const canonical_kp *)ckps,	   
	       NULL,
	       (const unsigned int)mywc(args->pri_file),
	       (const char *)args->pri_file,
	       &model,
	       stderr);      

    predict((const cmd_args *)args,		
	    (const double **)features,
	    (const hic *)data,
	    (const canonical_kp *)ckps,
	    (const boost *)model,	 
	    &pred,
	    stderr);

    pred_cmp_file((const cmd_args *)args,
		  (const hic *)data,
		  (const double *)pred,
		  stderr);


  }
  return 0;
}
