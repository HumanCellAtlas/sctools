#ifndef __SORT_WRITE_H__
#define __SORT_WRITE_H__

#include <string>
#include <vector>
#include <semaphore.h>

#include "gzstream.h"
#include "datatypes.h"
#include "utilities.h"
#include "globals.h"

/**
 * @brief This function takes a vector of tuples of the tags, sorts them
 * in the dictionary order of the tags and then writes these in the same 
 * order to a txt file
 *
 * @details
 * The function take the vector of tags tuples, writes the sorted tuples into 
 * a file. The filename is generated randomly (enought to avoid collision with other files)
 * in the temp folder specified.
 *
 * @param tuple_records: vector<TAGTUPLE> &, reference to a vector of TAGTUPLES
 * @return a string for the random file name
*/
void  write_out_partial_txt_file(const vector<TAGTUPLE> &tuple_records,  \
          std::string const & tmp_folder,  std::vector<string> &partial_files);

#endif
