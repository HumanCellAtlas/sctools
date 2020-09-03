/** 
 *  @file   main.cpp 
 *  @brief  The main function for the fastqprocessing step 
 *  @author Kishori Konwar 
 *  @date   2020-08-27 
 ***********************************************/

#include "fastqprocess.h"
#include "utilities.h"

/* Flag set by ‘--verbose’. */
int main (int argc, char **argv)
{

  INPUT_OPTIONS options;

  read_options(argc, argv, options);

  WHITE_LIST_DATA *white_list_data = read_white_list(options.white_list_file);

  process_inputs(options, white_list_data);

  return 0;

}

