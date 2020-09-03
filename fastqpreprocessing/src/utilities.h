#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__
#include "fastqprocess.h"
#include <math.h>
#include <cstdint>


typedef struct _input_options {
  _input_options() {
     barcode_length = -1;
     umi_length = -1;
     sample_id = "";
     bam_size = 1.0;
  }

  vector<string> I1s, R1s, R2s;
  string white_list_file;
  int barcode_length, umi_length; 
  double bam_size;
  string sample_id;
} INPUT_OPTIONS;


int process_inputs(const INPUT_OPTIONS &,  const WHITE_LIST_DATA *) ;
int32_t get_num_blocks(const INPUT_OPTIONS &options);

WHITE_LIST_DATA *read_white_list(const string &white_list_file);

void read_options(int, char **, INPUT_OPTIONS &); 
int32_t filesize(const char *filename);
int32_t getFileSize(const std::string &fileName);

void error (char *msg);

#endif
