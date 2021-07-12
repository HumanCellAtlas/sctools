/**
 *  @file   inputoptions.h
 *  @brief  Utility functions for input options processing
 *  @author Kishori Konwar
 *  @date   2020-08-26
 ***********************************************/
#ifndef __DATA_TYPES__
#define __DATA_TYPES__

#include <getopt.h>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

typedef std::tuple<std::string, std::string, std::string> TAGTUPLE;
typedef std::tuple<std::string, int, int>  QUEUETUPLE;

typedef std::pair<std::string, bool>  STRING_BOOL_PAIR;

typedef std::vector<std::string>  STRING_VECTOR;

typedef std::unordered_map <std::string, int64_t> STRING_INT32_MAP;


typedef struct _tags {
   char **tags;
} TAGS;

typedef struct _tags_holder {
    int num_tags;
    TAGS *tags;
    char *memorypool;

    char *allocated_memory(int size) {
       return 0;
    }

    char *double_memory() {
       return 0;
    } 
} TAGS_HOLDER;





// structure for correcting the barcodes
typedef struct _white_list_data {
    /// an unordered map from whitelist barcodes and 1-mutations
    /// to the index of the correct barcode
    STRING_INT32_MAP mutations;
    /// vector of whitelist barcodes
    STRING_VECTOR barcodes;
} WHITE_LIST_DATA;


/// Structure to hold input options for fastqprocess
typedef struct _input_options_fastqprocess {
  /// Initialize some of the values
  _input_options_fastqprocess() {
     barcode_length = -1;
     umi_length = -1;
     sample_id = "";
     bam_size = 1.0;
  }
  /// I1, R1 and R2 files name
  std::vector<std::string> I1s, R1s, R2s;
  /// Barcode white list file
  std::string white_list_file;
  //// chemistry dependent (V2/V3) barcode and UMI length
  int barcode_length, umi_length;
  /// Bam file size to split by (in GB)
  double bam_size;
  /// sample name
  std::string sample_id;
} INPUT_OPTIONS_FASTQPROCESS;


/**
 * @brief Reads the options to the fastqprocess program
 *
 * @param argc  no of arguments to the main function
 * @param argv arguments array to the main function
 * @param options the structure for holding the options for getopt
*/
void read_options_fastqprocess(int, char **, INPUT_OPTIONS_FASTQPROCESS &);


// Structure to hold input options for tagsort
typedef struct _input_options_tagsort {
  /// Initialize some of the values
  _input_options_tagsort() {
     temp_folder =  std::string("/tmp/");
     inmemory_chunk_size = 1.0;
  }
  /// name of the bam file 
  std::string bam_input;
  // temp folder for disk sorting
  std::string temp_folder;
  // output file
  std::string output_file;
  /// size of data to load for inmemory sorting (in GB)
  double inmemory_chunk_size;
  // tags to sort by
  std::vector<std::string> tags;

  // bam lib "htslib" or "libgenstat"
  std::string bamlib;
} INPUT_OPTIONS_TAGSORT;


/**
 * @brief Reads the options to the tagsort program
 *
 * @param argc  no of arguments to the main function
 * @param argv arguments array to the main function
 * @param options the structure for holding the options for getopt
*/
void read_options_tagsort(int, char **, INPUT_OPTIONS_TAGSORT &);

#endif
