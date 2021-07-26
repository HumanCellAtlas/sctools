#ifndef __HTSLIB_TAG_SORT__
#define __HTSLIB_TAG_SORT__

//#include "utilities.h"
#include <tuple>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unordered_map>
#include <memory>
#include <semaphore.h>

#include <htslib/sam.h>
#include "inputoptions.h"
#include "datatypes.h"
#include "sort_write.h"
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

namespace htslib {

std::vector<string> create_sorted_file_splits_htslib(INPUT_OPTIONS_TAGSORT &options);

}
 

#endif
