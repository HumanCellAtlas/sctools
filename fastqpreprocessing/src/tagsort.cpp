/**
 *  @file   fastqprocess.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/

#include "fastqprocess.h"
#include "utilities.h"
#include "inputoptions.h"
#include "tagsort.h"
#include <tuple>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>

#define STRING_LEN  40
using namespace std;

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
/*
std::string write_out_partial_txt_file(const vector<TAGTUPLE> &tuple_records, std::string const & tmp_folder) {
    std::string tempfile = tmp_folder + string("/") + random_string(STRING_LEN) + std::string(".txt");

    ofstream output_fp;
    output_fp.open (tempfile); 
    std::vector<std::pair<std::string, int>> index_pairs;
    int k  = 0;
    for( auto it=tuple_records.begin(); it!=tuple_records.end(); it++, k++) { 
       // 'UB:CB:GX"   or "GX:CB:UB"
       index_pairs.push_back(std::make_pair(get<0>(*it) + std::string("\t") + get<1>(*it) + std::string("\t") + get<2>(*it) , k));
    }

    // sort using the lambda function
    std::sort(index_pairs.begin(), index_pairs.end(), sortbyfirst);

    for (auto it=index_pairs.begin(); it!=index_pairs.end(); it++) {
       output_fp << get<0>(tuple_records[it->second]) << "\t"
                 << get<1>(tuple_records[it->second]) << "\t"
                 << get<2>(tuple_records[it->second]) << std::endl; 
    }
    output_fp.close();

    return tempfile;
}
*/


std::string write_out_partial_txt_file(const vector<TAGTUPLE> &tuple_records, std::string const & tmp_folder) {
    std::string tempfile = tmp_folder + string("/") + random_string(STRING_LEN) + std::string(".txt");
/*
    ofstream output_fp;
    output_fp.open (tempfile); 
    std::vector<std::pair<std::string, int>> index_pairs;
    int k  = 0;
    for( auto it=tuple_records.begin(); it!=tuple_records.end(); it++, k++) { 
       // 'UB:CB:GX"   or "GX:CB:UB"
       index_pairs.push_back(std::make_pair(get<0>(*it) + std::string("\t") + get<1>(*it) + std::string("\t") + get<2>(*it) , k));
    }

    // sort using the lambda function
    std::sort(index_pairs.begin(), index_pairs.end(), sortbyfirst);

    stringstream str(stringstream::out|stringstream::binary);

    for (auto it=index_pairs.begin(); it!=index_pairs.end(); it++) {
       str << get<0>(tuple_records[it->second]) << "\t"
                 << get<1>(tuple_records[it->second]) << "\t"
                 << get<2>(tuple_records[it->second]) << std::endl; 
    }
    //output_fp.write(str.str().c_str(), str.str().length());


    output_fp.close();
*/

    return tempfile;
}

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
std::vector<string> create_sorted_file_splits(INPUT_OPTIONS_TAGSORT &options) {
    string input_bam = options.bam_input;
    string tmp_folder = options.temp_folder;
    
    // size of individual chunks to sort in memory as an approx 20 mil alignments makes 1 GB bam
    int num_align_per_file =  static_cast<int>(options.inmemory_chunk_size * 2000000);
    std::vector<string> partial_files;

    SamFile samIn;
    samIn.OpenForRead(input_bam.c_str());

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    SamRecord samRecord;
   // Keep reading records until ReadRecord returns false.
    TAGTUPLE tagstuple = std::make_tuple("", "", "");

    long int i =0;
    const String *tagstr = NULL;
    vector<TAGTUPLE>  tuple_records;

    vector<SamRecord> samrecords;
    while(samIn.ReadRecord(samHeader, samRecord)) {
       tagstr = samRecord.getStringTag(options.tags[0].c_str());
       std::string a = tagstr!=NULL ? std::string(tagstr->c_str()) : std::string("");


       tagstr = samRecord.getStringTag(options.tags[1].c_str());
       std::string b = tagstr!=NULL ? std::string(tagstr->c_str()) : std::string("");

       tagstr = samRecord.getStringTag(options.tags[2].c_str());
       std::string c = tagstr!=NULL ? std::string(tagstr->c_str()) : std::string("");


       if (i!=0 && i%num_align_per_file==0) {
           std::cout << "Batch number : " << partial_files.size() << std::endl;
           std::string split_file_path = write_out_partial_txt_file(tuple_records, tmp_folder);
           tuple_records.clear();
           partial_files.push_back(split_file_path);
           //if (partial_files.size() ==7) break;
       }

       tuple_records.push_back(std::make_tuple(a, b, c));
       i = i + 1;
    }

    if (tuple_records.size()>0) {
        std::string split_file_path = write_out_partial_txt_file(tuple_records, tmp_folder);
        tuple_records.clear();
        partial_files.push_back(split_file_path);
     }

    std::cout << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;

    return partial_files;
}

void fill_buffer(Context &contx) {
    contx.data[contx.i].clear();
    int k = 0;

    for (std::string line; std::getline(*(contx.file_handles[contx.i]), line); k++) {
       contx.data[contx.i].push_back(line);
        if (k == contx.BUF_SIZE-1)
           break;
    }

    contx.data_size[contx.i] = contx.data[contx.i].size();
    if (contx.data_size[contx.i] != 0) {
       contx.ptrs[contx.i] = 0;
       contx.isempty[contx.i] = false;
    } else {
       contx.ptrs[contx.i] = contx.BUF_SIZE;
       contx.isempty[contx.i] = true;
    }
}


void merge_partial_files(const std::vector<std::string> &partial_files, 
                         const std::string &output_file, 
                         const std::vector<std::string> &tags) {

    Context contx; 
    contx.NUM_PARTS = partial_files.size();
    contx.BUF_SIZE = 100000;

    for (auto i=0; i < partial_files.size(); i++) {
       ifstream *input_fp = new ifstream;
       input_fp->open(partial_files[i]); 
       contx.file_handles.push_back(input_fp);
    }

    // set the isempty for each file to false
    for (auto i=0; i < partial_files.size(); i++) {
       contx.isempty.push_back(false);
    }

    // set a vector of vectors of data for each file
    contx.data.reserve(contx.NUM_PARTS);
    for (auto i=0; i < partial_files.size(); i++) {
       contx.data.push_back(std::vector<std::string>());
    }

    // set the data_size of the buffer for each file to 0
    for (auto i=0; i < partial_files.size(); i++) {
       contx.data_size.push_back(0);
    }

    // set the pointer to f each buffer to contxt.BUF_SIZE
    for (auto i=0; i < partial_files.size(); i++) {
       contx.ptrs.push_back(contx.BUF_SIZE);
    }

    //fill the buffers
    for (auto i=0; i < contx.NUM_PARTS; i++) {
        contx.i = i;
        fill_buffer(contx);
    }

    // create the heap from the first batch loaded data
    contx.num_active_files= 0;
    for (auto i=0; i< contx.NUM_PARTS; i++){
        contx.i = i;
        if (contx.ptrs[i] != contx.BUF_SIZE) {
           contx.num_active_files += 1;
           string comp_tag = contx.data[i][contx.ptrs[i]];
           contx.heap.push(QUEUETUPLE(comp_tag, i, contx.ptrs[i]));
        }
    }

    //  now merge by pop an push     
    ofstream fout;
    fout.open(output_file); 

    // pop and push from the heap
    int k = 0;
    int num_alignments = 0;
    while (contx.num_active_files > 0) {
        // pop the smallest
       QUEUETUPLE qtuple = contx.heap.top();
       auto val = get<0>(qtuple);
       auto i = get<1>(qtuple);
       auto j = get<2>(qtuple);
       contx.heap.pop();

       string field  = contx.data[i][j];
       num_alignments += 1;
       //fout << field << std::endl;
         
       /*   add a new element to the heapq
           if there is no data then fill it unless the file is empty
            BUF_SIZE = 5
     
            | | | |  data_size = 4
                   *
     
            | | | | |   
                 *
       */
        if (contx.isempty[i] == false && contx.ptrs[i] == contx.data_size[i]) {
            contx.i = i;
            fill_buffer(contx);
        } 
 
         // make sure it is not empty
        if (contx.data_size[i] > 0) {
            std::string comp_tag = contx.data[i][contx.ptrs[i]];
            contx.heap.push(QUEUETUPLE(comp_tag, i, contx.ptrs[i]));
            contx.ptrs[i] += 1;
        } else { // one more file is fully read 
            contx.num_active_files -= 1;
        }
    }
 

    // close the input files as they are empty
    for (auto i=0; i < contx.file_handles.size(); i++) {
       contx.file_handles[i]->close();
    }

    // close output files as there is no more to write
    fout.close();

    // we no longer need the partial files
    for (auto i=0; i < partial_files.size(); i++) {
        if(remove(partial_files[i].c_str()) != 0)
           std::cerr << string("Error deleting file") <<  partial_files[i] << std::endl;
    }

}

/* Flag set by ‘--verbose’. */
int main (int argc, char **argv)
{

  INPUT_OPTIONS_TAGSORT options;

  read_options_tagsort(argc, argv, options);

  std::cout << "bam input " << options.bam_input << std::endl;
  std::cout << "temp folder " << options.temp_folder << std::endl;
  std::cout << "output file " <<  options.output_file << std::endl;
  std::cout << "temp folder " << options.inmemory_chunk_size << std::endl;
  std::cout << "tags:" << std::endl;

  for(auto it = options.tags.begin(); it != options.tags.end(); it++) { 
      std::cout << "\t" << *it << std::endl;
  }

  /* first create a list of sorted, and simplified sorted files */
  std::vector<string> partial_files = create_sorted_file_splits(options);

  /* now merge the sorted files to create one giant sorted file by using 
    a head to compare the values based on the tags used  */

  merge_partial_files(partial_files, options.output_file, options.tags);


  return 0;
}

