#ifndef __HTSLIB_TAG_SORT__
#define __HTSLIB_TAG_SORT__

/**
 *  @file   htslib_tagsort.h
 *  @brief  Utility functions for input options processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <tuple>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unordered_map>
#include <memory>
#include <semaphore.h>
#include <set>

#include <htslib/sam.h>
#include "input_options.h"
#include "datatypes.h"
#include "sort_write.h"
#include "utilities.h"
#include "globals.h"


/**
 * @brief From the input bam create a list of txt files with the records (lines)
 * sorted according to the * tags
 *
 * @details
 * The input bam file is read chunk by chunk, sorted by the tags and the written
 * out as a text file in the sorted manner.
 *
 * @param options: INPUT_OPTIONS_TAGSORT the inputs to the program
 * @return a vector containing the file paths of the partial files
*/
void create_sorted_file_splits_htslib(INPUT_OPTIONS_TAGSORT& options);

#endif
