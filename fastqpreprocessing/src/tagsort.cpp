/**
 *  @file   fastqprocess.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/

#ifdef LIBSTATGEN
#include "libstatgen_tagsort.h"
#endif

#ifdef HTSLIB
#include "htslib_tagsort.h"
#endif

#include "tagsort.h"
#include <regex>

std::string special("TTTGTCATCTTGAGGT\tACTGTTTAAG\tENSMUSG00000071528.4,ENSMUSG00000071528.4\tchr19\tCODING\t47086155\t1\t72.125\t1\t38.5\t1\t1\t0\t0\t1");
int filling_counter = 0;
void fill_buffer(Context &contx) {
    contx.data[contx.i].clear();
    int k = 0;

    for (std::string line; k< contx.BUF_SIZE, std::getline(*(contx.file_handles[contx.i]), line); k++) {
       contx.data[contx.i].push_back(line);
       if (special==line) {
          std::cout << "filing special " << line << std::endl;
       }
       filling_counter += 1;
   //    if (k == contx.BUF_SIZE-1) break;
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
                         const std::string &output_file ) {

    Context contx; 
    contx.NUM_PARTS = partial_files.size();
    contx.BUF_SIZE = 100000;

    for (auto i=0; i < contx.NUM_PARTS; i++) {
       ifstream *input_fp = new ifstream;
       input_fp->open(partial_files[i]); 
       contx.file_handles.push_back(input_fp);
    }

    // set the isempty for each file to false
    for (auto i=0; i < contx.NUM_PARTS; i++) {
       contx.isempty.push_back(false);
    }

    // set a vector of vectors of data for each file
    for (auto i=0; i < contx.NUM_PARTS; i++) {
       contx.data.push_back(std::vector<std::string>());
    }

    // set the data_size of the buffer for each file to 0
    for (auto i=0; i < contx.NUM_PARTS; i++) {
       contx.data_size.push_back(0);
    }

    // set the pointer to f each buffer to contxt.BUF_SIZE
    for (auto i=0; i < contx.NUM_PARTS; i++) {
       contx.ptrs.push_back(contx.BUF_SIZE);
    }

    //fill the buffers
    for (auto i=0; i < contx.NUM_PARTS; i++) {
        contx.i = i;
        fill_buffer(contx);
    }

   std::regex rgx("\t");
   std::sregex_token_iterator end;

    // create the heap from the first batch loaded data
    contx.num_active_files= 0;
    for (auto i=0; i< contx.NUM_PARTS; i++){
        contx.i = i;
        if (contx.ptrs[i] != contx.BUF_SIZE) {
           contx.num_active_files += 1;
           std::sregex_token_iterator iter(contx.data[i][contx.ptrs[i]].begin(), contx.data[i][contx.ptrs[i]].end(), rgx, -1);
           std::stringstream comp_tag;
           for (auto k = 0; k < 3 && iter!=end; ++iter, k++) {
              if (k>0) comp_tag << "\t";
              comp_tag << *iter;
           }
           contx.heap.push(QUEUETUPLE(comp_tag.str(), i, contx.ptrs[i]++));
        }
    }
    //  now merge by pop an push     
    ofstream fout;
    fout.open(output_file); 

    // pop and push from the heap
    int k = 0;
    int num_alignments = 0;
    int i, j, i_1, j_1;

    stringstream str(stringstream::out|stringstream::binary);
    //while (contx.num_active_files > 0) {

    while (!contx.heap.empty()) {

       // read the top
       QUEUETUPLE qtuple = contx.heap.top();
       auto val = get<0>(qtuple);
       i = get<1>(qtuple);
       j = get<2>(qtuple);

        // pop it now
       contx.heap.pop();


       if (num_alignments%100000==0 && num_alignments > 0) {
           fout.write(str.str().c_str(), str.str().length());
           str.clear();
           str.str("");
       } 
       string field  = contx.data[i][j];
       //std::cout << "\twriting special " << field << std::endl;
       str << field << std::endl;
       num_alignments += 1;

       /*   add a new element to the heapq
            if there is no data then fill it unless the file is empty
       */

        // if ismpty is true means the file has been fully read 
        if (contx.isempty[i] == false && contx.ptrs[i] == contx.data_size[i]) {
            contx.i = i;
            fill_buffer(contx);
        } 

         // make sure it is not empty
        if (contx.data_size[i] > 0) {
           // std::string comp_tag = contx.data[i][contx.ptrs[i]];

           std::sregex_token_iterator iter(contx.data[i][contx.ptrs[i]].begin(), contx.data[i][contx.ptrs[i]].end(), rgx, -1);
           std::stringstream comp_tag;
           for (auto k = 0; k < 3 && iter!=end; ++iter, k++) {
              if (k>0) comp_tag << "\t";
              comp_tag << *iter;
           }

           contx.heap.push(QUEUETUPLE(comp_tag.str(), i, contx.ptrs[i]++));
        } else { // one more file is fully read 
            contx.num_active_files -= 1;
        }
    }
 
    // write out the remaining data
    fout.write(str.str().c_str(), str.str().length());

    // close the input files as they are empty
    for (auto i=0; i < contx.file_handles.size(); i++) {
       contx.file_handles[i]->close();
    }

    // close output files as there is no more to write
    fout.close();


    
    // we no longer need the partial files
    for (auto i=0; i < partial_files.size(); i++) {
   //     if(remove(partial_files[i].c_str()) != 0)
    //      std::cerr << string("Error deleting file") <<  partial_files[i] << std::endl;
    }
    std::cout << "Written "<< num_alignments << " alignments in total" << std::endl;

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

  for(auto it = options.tag_order.begin(); it != options.tag_order.end(); it++) { 
      std::cout << "\t" << it->first << "\t" << it->second << std::endl;
  }

  /* first create a list of sorted, and simplified sorted files */
  std::vector<string> partial_files;

#ifdef LIBSTATGEN
//  if (options.bamlib==std::string("LIBSTATGEN"))
      partial_files = libstatgen::create_sorted_file_splits_libstatgen(options);
#endif


#ifdef HTSLIB
 // if (options.bamlib==std::string("HTSLIB"))
      partial_files = htslib::create_sorted_file_splits_htslib(options);
#endif

  /* now merge the sorted files to create one giant sorted file by using 
    a head to compare the values based on the tags used  */

  merge_partial_files(partial_files, options.output_file);
  std::cout << "Aligments " <<  filling_counter << " loaded to buffer " << std::endl;


  return 0;
}

