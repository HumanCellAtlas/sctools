#ifndef __DATA_TYPES__
#define __DATA_TYPES__

/**
 *  @file   input_options.h
 *  @brief  Utility functions for input options processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <getopt.h>
#include <vector>
#include <string>
#include <unordered_map>
#include "globals.h"

using namespace std;

typedef std::tuple<std::string *, std::string *, std::string *> TRIPLET;

typedef struct TagCounter {
    TagCounter() {
      prev_tag = "";
    }
    std::unordered_map<std::string, int> data;
    std::string prev_tag;
    int count = 0;

    void clear() { data.clear(); }

    void update(const std::string &tag) {
       if (tag.compare(prev_tag)!=0)  {
          count++;
          prev_tag = tag;
       }
    }
} TAG_COUNTER;

typedef std::tuple<TRIPLET * /*  tuple<std::string *, std::string *, std::string *>*/,  
                   std::string /* reference */, 
                   std::string /* biotype */,  
                   int /* pos */, 
                   int /*rev strand   1 for yes, 0 otherwise*/, 
                   float /*avg barcode qual score */, 
                   float /* frac of barcode qual score >30 */, 
                   float /*avg qual seq */ , 
                   float /*fract of >30 score qual seq*/, 
                   int /*NH*/, 
                   int /*perfect molecule barcode, 1 is yes, 0 otherwise*/,
                   int /*spliced reads 1 yes, 0 otherwise*/,
                   int /*is duplicate */,
                   int /*perfect cell barcode 1 is yes, 0 otherwise*/,
                   float /* fraction of umi qual score > 30 */
                 > TAGTUPLE;

typedef std::tuple<std::string, int, int>  QUEUETUPLE;

typedef std::pair<std::string, bool> STRING_BOOL_PAIR;

typedef std::vector<std::string>  STRING_VECTOR;

typedef std::unordered_map <std::string, int64_t> STRING_INT64_MAP;


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
    STRING_INT64_MAP mutations;
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
     verbose_flag = 0;
  }
  //verbose flag
  unsigned int verbose_flag;
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

enum METRIC_TYPE {CELL, GENE};

// Structure to hold input options for tagsort
typedef struct _input_options_tagsort {
  /// Initialize some of the values
  _input_options_tagsort() {
     bam_input = "";
     gtf_file = "";
     temp_folder =  std::string("/tmp/");
     alignments_per_thread = NUM_ALNS_PER_THREAD;
     nthreads = 1;
     compute_metric = 0;
     output_sorted_info = 0;
     metric_type = "";
  }
  // metric type
  std::string metric_type;

  // output sorted info
  unsigned int output_sorted_info;
  // compute metric
  unsigned int compute_metric;
  /// name of the bam file 
  std::string bam_input;
  /// name of the gtf file 
  std::string gtf_file;
  // temp folder for disk sorting
  std::string temp_folder;
  // metric_output file
  std::string metric_output_file;
  // sorted tsv output file
  std::string sorted_output_file;
  // number of alignment per thread 
  unsigned int alignments_per_thread;
  // number of threads
  unsigned int nthreads;
  // barcode tag
  std::string barcode_tag;
  // umi tag
  std::string umi_tag;
  // gene tag
  std::string gene_tag;
  // order of the tags to sort by
  std::unordered_map<std::string, unsigned int> tag_order;

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
