/**
 *  @file   utilities.cpp
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include "utilities.h"

using namespace std;
namespace fs = std::experimental::filesystem;

/** @copydoc filesize */
int64_t filesize(const char *filename) {
    FILE *f = fopen(filename, "rb");  /* open the file in read only */

    int64_t size = 0;
    if (fseek(f, 0, SEEK_END) ==0 ) /* seek was successful */
       size = ftell(f);
    fclose(f);
    return size;
}

/** @copydoc read_white_list */
WHITE_LIST_DATA *read_white_list(const string &white_list_file) {
    char ATCG[] = {'A', 'C', 'G', 'T', 'N'};

    fstream newfile;
    WHITE_LIST_DATA *white_list_data = new  WHITE_LIST_DATA;

    // open a file to perform read operation using file object
    newfile.open(white_list_file, ios::in);
    int k = 0;
    if (newfile.is_open()) {   // checking whether the file is open
       string tp;

       // read data from file object and put it into string.
       while (getline(newfile, tp)) {
          //insert the barcode into the list
          white_list_data->barcodes.push_back(tp);

          for (unsigned int i=0; i < tp.size(); i++) {
            for (int j=0; j < 5; j++) {
              char c = tp[i];
              tp[i] = ATCG[j];
              /* if the  mutation in any of the positions
                 is already present update the new correct index,
                 else insert the barcode or update the index
                 This is done to have the same values for corrected barcodes
                 as in the python implementation
              */
              if (white_list_data->mutations.find(tp) ==
                  white_list_data->mutations.end()) {
                  white_list_data->mutations.insert({tp, k});
              } else {
                  white_list_data->mutations[tp] = k;
              }
              tp[i] = c;
            }
          }

          /* -1 suggests it is already a whitelisted barcode
             This is used, instead of the actual index, because when
             the barcode is seen with -1 then no correction is necessary.
             This is likely to avoid lots of map look up as most barcodes
             not erroneous.
          */
          white_list_data->mutations.at(tp) = -1;
          k++;
       }
       // close the file object.
       newfile.close();
    }
    return white_list_data;
}

/** @copydoc  _print_file_info */
void _print_file_info(const std::vector<std::string> &fastqs,
     const std::string &type) {
    if (fastqs.size()) {
        std::cout << "INFO " << type << " files:" << std::endl;
            for (unsigned int i= 0; i < fastqs.size(); i++) {
               if (fs::exists(fastqs[i].c_str())) {
                   std::cout << "\t " << fastqs[i]  <<  " exists, file size "
                        <<  filesize(fastqs[i].c_str())  <<  std::endl;
               } else {
                   std::cout << "ERROR " << fastqs[i] << " is missing!\n";
                   std::cerr << "ERROR " << fastqs[i] << " is missing!\n";
                   exit(1);
               }
           }
       }
}


/** @copydoc  get_num_blocks */
int64_t get_num_blocks(const INPUT_OPTIONS_FASTQPROCESS &options) {
    double tot_size = 0;
    for (unsigned int i= 0; i < options.R1s.size(); i++) {
        if (options.I1s.size()) {
            tot_size +=  filesize(options.I1s[i].c_str());
        }
        tot_size +=  filesize(options.R1s[i].c_str());
        tot_size +=  filesize(options.R2s[i].c_str());
    }

    return ceil((tot_size/(1024*1024*1024))
           /static_cast<double>(options.bam_size));
}

/** @copydoc error */
void error(char *msg) {
    perror(msg);
    exit(1);
}

/** @copydoc random_string **/
std::string random_string(size_t length)
{
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length, 0);
    std::generate_n(str.begin(), length, randchar);
    return str;
}


/** @copydoc read_lines **/
std::vector<std::string> read_lines(const std::string &file_name) {
  string line;
  std::vector<std::string> lines;
  ifstream myfile(file_name.c_str());
  if (myfile.is_open()) {
    while ( getline (myfile,line) ) {
      lines.push_back(line);
    }
    myfile.close();
  } else {
    throw std::runtime_error("Could not open file");
  }
  return lines;
}

/** @copydoc error_message **/
void error_message(const char *msg) {
      std::cerr << msg;
}


