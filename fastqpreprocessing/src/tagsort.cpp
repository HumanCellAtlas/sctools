/**
 *  @file   tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori M. Konwar
 *  @date   2020-08-27
 ***********************************************/

#include "htslib_tagsort.h"
#include "tagsort.h"
#include <regex>


extern std::vector<string> partial_files;
int filling_counter = 0;

/*
 * @brief fills the buffer for the files
 *
 * @param contx is the context of the file
*/
void fill_buffer(Context &contx) {
    contx.data[contx.i].clear();
    int k = 0;

    // the order of the loop condition is iportant first make sure if you can accomodate then try to read, 
    // otherwise it might create a read but never processed
    for (std::string line; k < contx.BUF_SIZE &&  std::getline(*(contx.file_handles[contx.i]), line); k++) {
       contx.data[contx.i].push_back(line);
       filling_counter += 1;
    }

//    std::cout << " Filling buffer " << contx.i << " with " << contx.data[contx.i].size() << " lines " << std::endl;
    contx.data_size[contx.i] = contx.data[contx.i].size();

    if (contx.data_size[contx.i] != 0) {
       contx.ptrs[contx.i] = 0;
       contx.isempty[contx.i] = false;
    } else {
       contx.ptrs[contx.i] = contx.BUF_SIZE;
       contx.isempty[contx.i] = true;
    }

#ifdef DEBUG
    std::cout << "-->" << std::endl;
    for (int m = 0; m < contx.NUM_PARTS; m++) {
       std::cout << "\t" << m << " : " << contx.data_size[m] << " : " << contx.ptrs[m] << std::endl;
    }
#endif

}

/*
 * @brief Merges the files that are already sorted
 *
 * @param partial_files
 * @param output_file
*/

bool merge_partial_files(const std::string &output_file ) {

    // input the buffer size and partial files
    Context contx(partial_files.size(), DATA_BUFFER_SIZE); 
    auto cmp = [](const QUEUETUPLE &a, const  QUEUETUPLE &b) {
        return get<0>(a) > get<0>(b);
     };

    std::priority_queue<QUEUETUPLE, std::vector<QUEUETUPLE>,  decltype(cmp) > heap(cmp);
 
    // open files 
    for (auto i=0; i < contx.NUM_PARTS; i++) {
#ifdef GZSTREAM
       igzstream *input_fp = new igzstream;
#else
       ifstream *input_fp = new ifstream;
#endif

       // set exceptionss to be thrown on failure
       //input_fp->exceptions(std::ifstream::failbit | std::ifstream::badbit);
       input_fp->open(partial_files[i].c_str()); 
       if(!input_fp->is_open()) {
            std::cerr << "ERROR failed to open the file " << partial_files[i] << std::endl;
            std::cerr << "      consider increasing the number of alignments per thread." << std::endl;
            std::cerr << "      That might require you to run of on a machine with more RAM." << std::endl;
            return false;
       }

       contx.file_handles.push_back(input_fp);
    }

    //fill the buffers
    for (auto i=0; i < contx.NUM_PARTS; i++) {
        contx.i = i;
        fill_buffer(contx);
    }
    
    std::regex rgx("\t");
    std::sregex_token_iterator end;

    // create the heap from the first batch loaded data
    contx.num_active_files = 0;
    for (auto i=0; i< contx.NUM_PARTS; i++){
        contx.i = i;
        if (contx.ptrs[i] != contx.BUF_SIZE) {
           std::sregex_token_iterator iter(contx.data[i][contx.ptrs[i]].begin(), contx.data[i][contx.ptrs[i]].end(), rgx, -1);
           std::stringstream comp_tag;
           for (auto k = 0; k < 3 && iter!=end; ++iter, k++) {
              if (k>0) comp_tag << "\t";
              comp_tag << *iter;
           }
           heap.push(QUEUETUPLE(comp_tag.str(), i, contx.ptrs[i]));
           contx.ptrs[i]++;
           contx.num_active_files += 1;
        }
    }

    //  now merge by pop an push     
#ifdef GZSTREAM
    ogzstream  fout;
#else
    ofstream fout;
#endif
    fout.open(output_file.c_str()); 

    // pop and push from the heap
    int num_alignments = 0;
    int i, j;

    stringstream str(stringstream::out|stringstream::binary);
    //while (contx.num_active_files > 0) {
    std::string  prev_comp_tag = "";
    while (!heap.empty()) {
       // read the top
       QUEUETUPLE qtuple = heap.top();
       std::string curr_comp_tag = get<0>(qtuple);
       assert(prev_comp_tag.compare(curr_comp_tag) <= 0);

#ifdef DEBUG
       contx.print_status();
       if (prev_comp_tag.compare(curr_comp_tag) <= 0 )  {
          std::cout << "Expecte " << prev_comp_tag << "\n\t\t" << curr_comp_tag << std::endl;
       } else {
          std::cout << "Anomaly " << prev_comp_tag << "\n\t\t" << curr_comp_tag << std::endl;
          exit(0);
       }
#endif

       i = get<1>(qtuple);  //buffer no
       j = get<2>(qtuple);  //the pointer into the ith buffer array

       // pop it now
       heap.pop();

       // start writing in chunks from the stream buffer
       if (num_alignments%contx.BUF_SIZE==0) {
           fout.write(str.str().c_str(), str.str().length());
           str.clear();
           str.str("");
       } 

       // load into stream buffer
       string field  = contx.data[i][j];
       str << field << std::endl;
       num_alignments += 1;

        // if ismpty is true means the file has been fully read 
       if (contx.isempty[i] == false && contx.ptrs[i] == contx.data_size[i]) {
           contx.i = i;
           fill_buffer(contx);
       } 

       //if (num_alignments%10000 == 0 )
       //std::cout << " heap size " << heap.size() << std::endl;
       /*  add a new element to the heapq
           if there is no data then fill it unless the file is empty
       */
       // make sure it is not empty
       if (contx.data_size[i] > 0) {
          // std::string comp_tag = contx.data[i][contx.ptrs[i]];
           std::sregex_token_iterator iter(contx.data[i][contx.ptrs[i]].begin(), contx.data[i][contx.ptrs[i]].end(), rgx, -1);
           std::stringstream comp_tag;
           for (auto k = 0; k < 3 && iter!=end; ++iter, k++) {
              if (k>0) comp_tag << "\t";
              comp_tag << *iter;
           }
           heap.push(QUEUETUPLE(comp_tag.str(), i, contx.ptrs[i]));
           contx.ptrs[i]++;
       } else { // one more file is fully read 
           contx.num_active_files -= 1;
       }
       prev_comp_tag = curr_comp_tag;
    }
 
    // write out the remaining data
    fout.write(str.str().c_str(), str.str().length());
    str.str("");
    str.clear();

    // close the input files as they are empty
    for (unsigned int i=0; i < contx.file_handles.size(); i++) {
       contx.file_handles[i]->close();
    }

    // close output files as there is no more to write
    fout.close();

    std::cout << "Written "<< num_alignments << " alignments in total" << std::endl;
    contx.clear();
    return true;
}

/* Flag set by ‘--verbose’. */
int main (int argc, char **argv)
{

  INPUT_OPTIONS_TAGSORT options;

  read_options_tagsort(argc, argv, options);

  std::cout << "bam input " << options.bam_input << std::endl;
  std::cout << "temp folder " << options.temp_folder << std::endl;
  std::cout << "output file " <<  options.output_file << std::endl;
  std::cout << "temp folder " << options.alignments_per_thread << std::endl;
  std::cout << "tags:" << std::endl;

  for(auto it = options.tag_order.begin(); it != options.tag_order.end(); it++) { 
      std::cout << "\t" << it->first << "\t" << it->second << std::endl;
  }

  /* first create a list of sorted, and simplified sorted files */
  htslib::create_sorted_file_splits_htslib(options);

  /* now merge the sorted files to create one giant sorted file by using 
    a head to compare the values based on the tags used  */
  std::cout << "Merging " <<  partial_files.size() << " sorted files!"<< std::endl;

  if(!merge_partial_files(options.output_file)) {
     std::cout << "Failed to complete the merging as the number of concurrently open files increased the max limit" << std::endl;
  }

  // we no longer need the partial files
  for (unsigned int i=0; i < partial_files.size(); i++) {
      if(remove(partial_files[i].c_str()) != 0)
        std::cerr << string("Error deleting file") <<  partial_files[i] << std::endl;
  }
  partial_files.clear();
  std::cout << "Aligments " <<  filling_counter << " loaded to buffer " << std::endl;

  return 0;
}

