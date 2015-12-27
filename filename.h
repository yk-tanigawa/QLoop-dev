#ifndef __FILENAME_H__
#define __FILENAME_H__

#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include "constant.h"
#include "calloc_errchk.h"

typedef struct _filenames {
  char *kmer_freq;
  char *hic;
  char *histo;
  char *adaboost;
  char *qp_P;
  char *qp_q;
} filenames;

int show_filenames(FILE *fp,
		   filenames *fnames){
  fprintf(fp, "%s\n", fnames->kmer_freq);
  fprintf(fp, "%s\n", fnames->hic);
  fprintf(fp, "%s\n", fnames->histo);
  fprintf(fp, "%s\n", fnames->adaboost);
  fprintf(fp, "%s\n", fnames->qp_P);
  fprintf(fp, "%s\n", fnames->qp_q);
  return 0;
}


int set_filenames(const command_line_arguements *args,
		  filenames **fnames){
  char *header;
  *fnames = calloc_errchk(1, sizeof(filenames), "calloc: filenames");

  /**
./main: info: fasta file: /data/yt/GRCh37.ch21.fasta
./main: info: Hi-C raw data directory: /data/yt/GM12878_combined/
  */

  { /* kmer_freq */
    char buf[F_NAME_LEN];
    char *pos;

    (*fnames)->kmer_freq = calloc_errchk(F_NAME_LEN, sizeof(char),
					 "fnames->kmer_freq");

    {
      strncpy(buf, basename(args->fasta_file), F_NAME_LEN);
      if((pos = strstr(buf, ".fasta")) != NULL){
	*pos = '\0';
      }else if((pos = strstr(buf, ".fa")) != NULL){
	*pos = '\0';
      }
    }

    sprintf((*fnames)->kmer_freq, "%s/%s.k%d.res%dk.freq", 
	    args->output_dir, buf, args->k, (args->res) / 1000);	    	   
  }

  { /* common header */
    header = calloc_errchk(F_NAME_LEN, sizeof(char), "calloc: fnames: header");
    sprintf(header, 
	    "%s/chr%d.m%dk.M%dk.%s.%s",
	    args->output_dir,
	    args->chr,
	    (args->min_size) / 1000,
	    (args->max_size) / 1000,
	    args->norm,
	    args->exp);
  }

  { /* Hi-C */
    (*fnames)->hic = calloc_errchk(F_NAME_LEN, sizeof(char),
					"fnames->hic");
    sprintf((*fnames)->hic, "%s.hic", header);
  }

  { /* histo */
    (*fnames)->histo = calloc_errchk(F_NAME_LEN, sizeof(char),
					"fnames->histo");
    sprintf((*fnames)->histo, "%s.histo", header);
  }

  { /* AdaBoost */
    (*fnames)->adaboost = calloc_errchk(F_NAME_LEN, sizeof(char),
					"fnames->adaboost");
    sprintf((*fnames)->adaboost, "%s.k%d.res%dk.p%d.T%ld.stamps", header, 
	    args->k, (args->res) / 1000, (int)(100 * (args->percentile)),
	    args->iteration_num);
  }
  { /* QP */
    (*fnames)->qp_P = calloc_errchk(F_NAME_LEN, sizeof(char),
				    "fnames->qp_P");
    (*fnames)->qp_q = calloc_errchk(F_NAME_LEN, sizeof(char),
				    "fnames->qp_q");
    sprintf((*fnames)->qp_P, "%s.k%d.res%dk.p%d.T%ld.P", header,
	    args->k, (args->res) / 1000, (int)(100 * (args->percentile)),
	    args->iteration_num);
    sprintf((*fnames)->qp_q, "%s.k%d.res%dk.p%d.T%ld.q", header,
	    args->k, (args->res) / 1000, (int)(100 * (args->percentile)),
	    args->iteration_num);
  }

  return 0;
}


#endif
