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

int filling_counter = 0;
void fill_buffer(Context &contx) {
    contx.data[contx.i].clear();
    int k = 0;

    for (std::string line; k< contx.BUF_SIZE, std::getline(*(contx.file_handles[contx.i]), line); k++) {
       contx.data[contx.i].push_back(line);
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
    int i, j;
    stringstream str(stringstream::out|stringstream::binary);
    //while (contx.num_active_files > 0) {
    while (!contx.heap.empty()) {
        // pop the smallest
       QUEUETUPLE qtuple = contx.heap.top();
       auto val = get<0>(qtuple);
       i = get<1>(qtuple);
       j = get<2>(qtuple);
       contx.heap.pop();
       num_alignments += 1;

       string field  = contx.data[i][j];
       if (num_alignments%1000000==0) {
//           std::cout << "writing "<< num_alignments << "  of size " << str.str().length() << std::endl;
           fout.write(str.str().c_str(), str.str().length());
           str.clear();
           str.str("");
       } else {
           str << field << std::endl;
       }
       /*   add a new element to the heapq
           if there is no data then fill it unless the file is empty
            BUF_SIZE = 5
     
            | | | |  data_size = 4
                   *
     
            | | | | |   
                 *
       */
        // if ismpty is true means the file has been fully read 
        if (contx.isempty[i] == false && contx.ptrs[i] == contx.data_size[i]-1) {
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
        if(remove(partial_files[i].c_str()) != 0)
           std::cerr << string("Error deleting file") <<  partial_files[i] << std::endl;
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

  for(auto it = options.tags.begin(); it != options.tags.end(); it++) { 
      std::cout << "\t" << *it << std::endl;
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

  merge_partial_files(partial_files, options.output_file, options.tags);
  std::cout << "Aligments " <<  filling_counter << " loaded to buffer " << std::endl;


  return 0;
}

