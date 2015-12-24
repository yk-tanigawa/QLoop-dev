#ifndef __SHOW_MSG__
#define __SHOW_MSG__

#include <stdio.h>

int show_usage(FILE *fp, 
	       const char *prog_name){
  fprintf(fp, "%s:Info: Usage\n", prog_name);
  return 0;
}

int show_error(FILE *fp,
	       const char *prog_name,
	       const char *msg){
  fprintf(fp, "%s: error: %s\n", prog_name, msg);
  return 0;
}

int show_warning(FILE *fp,
		 const char *prog_name,
		 const char *msg){
  fprintf(fp, "%s: warning: %s\n", prog_name, msg);
  return 0;
}

int show_info(FILE *fp,
		 const char *prog_name,
		 const char *msg){
  fprintf(fp, "%s: info: %s\n", prog_name, msg);
  return 0;
}


#endif
