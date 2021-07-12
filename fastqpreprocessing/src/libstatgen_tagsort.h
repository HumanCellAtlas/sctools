#ifndef __LIBSTATGEN_TAG_SORT__
#define __LIBSTATGEN_TAG_SORT__

#include "utilities.h"
#include "inputoptions.h"

#include <tuple>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "utilities.h"
#include "inputoptions.h"
#include "BaseAsciiMap.h"
#include "SamFile.h"
#include "SamValidation.h"


#ifdef LIBSTATGEN
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
using namespace std;
namespace libstatgen {
std::vector<std::string> create_sorted_file_splits_libstatgen(INPUT_OPTIONS_TAGSORT &options);
}
#endif

#endif
