#include <stdio.h>

#include "src/cmd_args.h"

int main(int argc, char **argv){  
  cmd_args *args;

  cmd_args_parse(argc, argv, &args);

  cmd_args_chk(args);



  return 0;
}
