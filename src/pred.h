#ifndef __PRED_H__
#define __PRED_H__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "constant.h"
#include "cmd_args.h"
#include "mywc.h"
#include "calloc_errchk.h"

/* Hi-C data */
typedef struct _hic {
  unsigned long nrow;
  unsigned int *i;
  unsigned int *j;
  double *mij;
} hic;

int hic_read(const cmd_args *, hic **);
	     
/**
 * read Hi-C data from a file 
 */

int hic_read(const cmd_args *args,
	     hic **data){

  {
    /* allocate memory */
    *data         = calloc_errchk(1, sizeof(hic), "calloc hic");   
    (*data)->nrow = mywc(args->hic_file);  
    (*data)->i    = calloc_errchk((*data)->nrow, sizeof(unsigned int),
				  "calloc hic (*data)->i");
    (*data)->j    = calloc_errchk((*data)->nrow, sizeof(unsigned int),
				  "calloc hic (*data)->j");
    (*data)->mij  = calloc_errchk((*data)->nrow, sizeof(double), 
				  "calloc hic (*data)->mij");
  }

  /* read from a file */
  {
    FILE *fp;
    char buf[BUF_SIZE], tmp_mij_str[BUF_SIZE];
    unsigned int tmp_i = 0, tmp_j = 0;
    unsigned long row = 0;

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "start reading Hi-C file from %s\n",
	    args->hic_file);

    if((fp = fopen(args->hic_file, "r")) == NULL){
      fprintf(stderr, "error: fopen %s\n%s\n",
	      args->hic_file, strerror(errno));
      exit(EXIT_FAILURE);
    }
    
    while(fgets(buf, BUF_SIZE, fp) && row < (*data)->nrow){
      sscanf(buf, "%d\t%d\t%s", &tmp_i, &tmp_j, (char *)(&tmp_mij_str));	 
      ((*data)->mij)[row] = strtod(tmp_mij_str, NULL);	  
      if(tmp_i <= tmp_j){
	((*data)->i)[row] = tmp_i / args->res;
	((*data)->j)[row] = tmp_j / args->res;
      }else{
	((*data)->i)[row] = tmp_j / args->res;
	((*data)->j)[row] = tmp_i / args->res;
      }
      row++;
    }
  
    fclose(fp);

    fprintf(stderr, "%s [INFO] ", args->prog_name);
    fprintf(stderr, "# of Hi-C data points = %ld\n",
	    (*data)->nrow);

  }

  return 0;
}

#endif
