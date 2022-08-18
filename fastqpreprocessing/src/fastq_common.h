#include "FastQFile.h"
#include "FastQStatus.h"
#include "SamFile.h"
#include "SamValidation.h"

#include <functional>
#include <string>
#include <vector>

void fillSamRecordCommon(SamRecord* samRecord, FastQFile* fastQFileI1,
                         FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                         bool has_I1_file_list,
                         std::string const& barcode_seq, std::string const& barcode_quality,
                         std::string const& umi_seq, std::string const& umi_quality);

void mainCommon(
    std::string white_list_file, int num_writer_threads, std::string output_format,
    std::vector<std::string> I1s, std::vector<std::string> R1s, std::vector<std::string> R2s,
    std::function <void(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> sam_record_filler,
    std::function <std::string(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> barcode_getter);
