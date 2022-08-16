#ifndef __SORT_WRITE_H__
#define __SORT_WRITE_H__
/**
 *  @file   sort_write.h
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/
#include <string>
#include <vector>
#include <semaphore.h>

#include "utilities.h"
#include "globals.h"

using TRIPLET = std::tuple<std::string*, std::string*, std::string*>;

using TAGTUPLE = std::tuple<
    TRIPLET* /*  tuple<std::string *, std::string *, std::string *>*/,
    std::string /* reference */,
    std::string /* biotype */,
    int /* pos */,
    int /*rev strand   1 for yes, 0 otherwise*/,
    float /*avg barcode qual score */,
    float /* frac of barcode qual score >30 */,
    float /*avg qual seq */,
    float /*fract of >30 score qual seq*/,
    int /*NH*/,
    int /*perfect molecule barcode, 1 is yes, 0 otherwise*/,
    int /*spliced reads 1 yes, 0 otherwise*/,
    int /*is duplicate */,
    int /*perfect cell barcode 1 is yes, 0 otherwise*/,
    float /* fraction of umi qual score > 30 */
    >;

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
void write_out_partial_txt_file(std::vector<TAGTUPLE> const& tuple_records,
                                std::string const& tmp_folder);

#endif
