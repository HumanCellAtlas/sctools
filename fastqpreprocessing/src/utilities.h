#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__

/**
 *  @file   utilities.h
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

// structure for correcting the barcodes
struct WHITE_LIST_DATA
{
  // an unordered map from whitelist barcodes and 1-mutations
  // to the index of the correct barcode
  std::unordered_map<std::string, int64_t> mutations;
  // vector of whitelist barcodes
  std::vector<std::string> barcodes;
};

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
std::unique_ptr<WHITE_LIST_DATA> readWhiteList(std::string const& white_list_file);

/**
 * @brief Print system error and exit
 *
 * @param msg  error string to print
*/
void crashWithPerror(std::string msg);

void crash(std::string msg);

#endif
