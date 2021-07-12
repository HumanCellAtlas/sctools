/**
 *  @file   fastqprocess.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/
//#include "libstatgen_tagsort.h"

// this file is a part of the archive file so these header do not
// need to be visible in the .h file above and creates name conflicts

#include "libstatgen_tagsort.h"


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
namespace libstatgen {
std::vector<string> create_sorted_file_splits_libstatgen(INPUT_OPTIONS_TAGSORT &options) {
    string input_bam = options.bam_input;
    string tmp_folder = options.temp_folder;
    
    // size of individual chunks to sort in memory as an approx 20 mil alignments makes 1 GB bam
    int num_align_per_file =  static_cast<int>(options.inmemory_chunk_size * 1000000);
    std::vector<string> partial_files;

    SamFile samIn;
    samIn.OpenForRead(input_bam.c_str());

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    SamRecord samRecord;

    long int i =0;
    const String *tagstr = NULL;
    vector<TAGTUPLE>  tuple_records;

    vector<SamRecord> samrecords;
    std::cout << "Running libstatgen" << std::endl;
   // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord)) {

       tagstr = samRecord.getStringTag(options.tags[0].c_str());
       std::string a = tagstr!=NULL ? std::string(tagstr->c_str()) : std::string("");

       tagstr = samRecord.getStringTag(options.tags[1].c_str());
       std::string b = tagstr!=NULL ? std::string(tagstr->c_str()) : std::string("");

       tagstr = samRecord.getStringTag(options.tags[2].c_str());
       std::string c = tagstr!=NULL ? std::string(tagstr->c_str()) : std::string("");

       if (i!=0 && i%num_align_per_file==0) {
           std::cout << "Batch number : " << partial_files.size() << std::endl;
           std::string split_file_path = "hello"; //write_out_partial_txt_file(tuple_records, tmp_folder);
           tuple_records.clear();
           partial_files.push_back(split_file_path);
          // if (partial_files.size() ==7) break;
       }

       tuple_records.push_back(std::make_tuple(a, b, c));
       i = i + 1;
    }

    if (tuple_records.size()>0) {
        std::string split_file_path = "hey"; // write_out_partial_txt_file(tuple_records, tmp_folder);
        tuple_records.clear();
        partial_files.push_back(split_file_path);
     }

    std::cout << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;

    return partial_files;
}

}
#endif
