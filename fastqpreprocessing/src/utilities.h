#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__

/**
 *  @file   utilities.h
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <math.h>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <experimental/filesystem>

#include "datatypes.h"
/**
 * @brief Compute the number of bam files
 *
 * @details
 *  Adds up the size of the input file and divides by the size of
 * user specified bam file size

 * @param options Input options structure that contains file name
*/
int64_t get_num_blocks(const INPUT_OPTIONS_FASTQPROCESS &options);

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
 *  @brief Computes the size of a file in bytes
 *
 *  @param filename file name whose size is computed
 *  @return size of the file in bytes
*/
int64_t filesize(const char *filename);


/**
 *  @brief Computes the size of a file in bytes
 *
 *  @param filename file name whose size is computed
 *  @return size of the file in bytes
*/
int64_t getFileSize(const std::string &fileName);

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


/**
 * @brief this function generates a random string  of a specified length
 * consisting of alphanumeric characters
 *
 * @param length: length of the string
 * @return a random alphanumeric string of specified length
*/
std::string random_string(size_t length);

/**
 * @brief this function reads the lines in a  text file into a vector 
 * of strings
 *
 * @param file_name: file name
 * @return a vector of strings
*/

std::vector<std::string> read_lines(const std::string &file_name);

template<typename T>
inline void freeStlContainer(T& p_container)
{
   return;
/*
   T empty;
   using std::swap;
   swap(p_container, empty);
*/
}

#endif
