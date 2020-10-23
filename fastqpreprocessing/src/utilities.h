/**
 *  @file   utilities.h
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-26
 ***********************************************/
#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__
#include <math.h>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <getopt.h>

typedef std::pair<std::string, bool>  STRING_BOOL_PAIR;

typedef std::vector<std::string>  STRING_VECTOR;

typedef std::unordered_map <std::string, int32_t> STRING_INT32_MAP;

// structure for correcting the barcodes
typedef struct _white_list_data {
    /// an unordered map from whitelist barcodes and 1-mutations
    /// to the index of the correct barcode
    STRING_INT32_MAP mutations;
    /// vector of whitelist barcodes
    STRING_VECTOR barcodes;
} WHITE_LIST_DATA;

/// Structure to hold input options
typedef struct _input_options {
  /// Initialize some of the values
  _input_options() {
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
} INPUT_OPTIONS;


/**
 * @brief Compute the number of bam files
 *
 * @details
 *  Adds up the size of the input file and divides by the size of
 * user specified bam file size

 * @param options Input options structure that contains file name
*/
int32_t get_num_blocks(const INPUT_OPTIONS &options);

/**
 * @brief Build barcode correction map white list barcodes & mutations
 *
 * @details
 * A barcode is computed by checking if it is either in the white
 * list or 1-mutation away from any white listed barcode. To check
 * whether a barcode is correct or to correct it, if 1-mutation away from
 * a barcode in the white list, we build a
 * a map is created with the barcodes and the 1-mutation. The keys are
 * barcodes or mutation and the values are index of the crrect barcode
 *
 * @param whilte_list_file  white list file from 10x genomics' cellranger
 * @return a stricture containing the barcode/1-mutation barcode to index
 *         of the correct barcode
*/
WHITE_LIST_DATA *read_white_list(const std::string &white_list_file);

/**
 * @brief Reads the options to the program
 *
 * @param argc  no of arguments to the main function
 * @param argv arguments array to the main function
 * @param options the structure for holding the options for getopt
*/
void read_options(int, char **, INPUT_OPTIONS &);

/**
 *  @brief Computes the size of a file in bytes
 *
 *  @param filename file name whose size is computed
 *  @return size of the file in bytes
*/
int32_t filesize(const char *filename);


/**
 *  @brief Computes the size of a file in bytes
 *
 *  @param filename file name whose size is computed
 *  @return size of the file in bytes
*/
int32_t getFileSize(const std::string &fileName);

/**
 * @brief Print system error and exit
 *
 * @param msg  error string to print
*/
void error(char *msg);

/**
 * @brief examines the existence and size of the input files
 *
*/
void _print_file_info(const std::vector<std::string> &fastqs, \
    const std::string &type);

#endif
