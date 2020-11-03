/**
 *  @file   utilities.cpp
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/

#include "utilities.h"
#include <string>
#include <cstdint>
#include <experimental/filesystem>

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

          for (int i=0; i < tp.size(); i++) {
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

/** @copydoc read_options */
void read_options(int argc, char **argv, INPUT_OPTIONS &options) {
  int c;
  int i;

  int verbose_flag = 0;

  static struct option long_options[] = {
          /* These options set a flag. */
          {"verbose",           no_argument,       0, 'v'},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"barcode-length",    required_argument, 0, 'b'},
          {"umi-length",        required_argument, 0, 'u'},
          {"bam-size",          required_argument, 0, 'B'},
          {"sample-id",         required_argument, 0, 's'},
          {"I1",                required_argument, 0, 'I'},
          {"R1",                required_argument, 0, 'R'},
          {"R2",                required_argument, 0, 'r'},
          {"white-list",        required_argument, 0, 'w'},
          {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char *help_messages[] = {
           "verbose messages  ",
           "barcode length [required]",
           "UMI length [required]",
           "output BAM file in GB [optional: default 1 GB]",
           "sample id [required]",
           "I1 [optional]",
           "R1 [required]",
           "R2 [required]",
           "whitelist (from cellranger) of barcodes [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "b:u:B:s:I:R:r:w:v",
                          long_options,
                          &option_index)) !=- 1
                         )
  {
      // process the option or arguments
      switch (c) {
        case 'v':
            verbose_flag = 1;
            break;
        case 0:
          /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;
        case 'b':
            options.barcode_length = atoi(optarg);
            break;
        case 'u':
            options.umi_length = atoi(optarg);
            break;
        case 'B':
            options.bam_size = atof(optarg);
            break;
        case 's':
            options.sample_id = string(optarg);
            break;
        case 'I':
            options.I1s.push_back(string(optarg));
            break;
        case 'R':
            options.R1s.push_back(string(optarg));
            break;
        case 'r':
            options.R2s.push_back(string(optarg));
            break;
        case 'w':
            options.white_list_file = string(optarg);
            break;
        case '?':
        case 'h':
          i = 0;
          printf("Usage: %s [options] \n", argv[0]);
          while (long_options[i].name != 0) {
            printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name,
                   long_options[i].has_arg == no_argument?
                  "no argument" : "required_argument",
                  help_messages[i]);
            i = i + 1;
          }
          /* getopt_long already printed an error message. */
          return;
          break;
        default:
          abort();
        }
    }

  // Check the options
  // number of R1 and R2 files should be equal
  bool exit_with_error = false;
  if ((options.R1s.size() != options.R2s.size())) {
     std::cout << "ERROR: Unequal number of R1 and R2 fastq files in input: "
         << "R1 : " << options.R1s.size()
         << "R2 : " << options.R2s.size()
         << std::endl;

     std::cerr << "ERROR: Unequal number of R1 and R2 fastq files in input: "
         << "R1 : " << options.R1s.size()
         << "R2 : " << options.R2s.size()
         << std::endl;

     exit_with_error = true;
  }

  if (options.R1s.size() == 0) {
     std::cout << "ERROR: No R1 file provided\n";
     std::cerr << "ERROR: No R1 file provided\n";

     exit_with_error = true;
  }


  if ((options.I1s.size() != options.R1s.size()) && (options.I1s.size() != 0)) {
     std::cout << "ERROR: Either the number of I1 input files are equal\n"
                  "       to the number of R1 input files, or no I1 input files\n"
                  "       should not be provided at all.\n";
     std::cerr << "ERROR: Either the number of I1 input files are equal\n"
                  "       to the number of R1 input files, or no I1 input files\n"
                  "       should not be provided at all.\n";

     exit_with_error = true;
  }
  // Bam file size must be positive
  if (options.bam_size <= 0) {
     std::cout << "ERROR: Size of a bam file (in GB) cannot be negative\n";
     std::cerr << "ERROR: Size of a bam file (in GB) cannot be negative\n";
     exit_with_error = true;
  }

  // must have a sample id
  if (options.sample_id.size() == 0) {
     std::cout << "ERROR: Must provide a sample id or name\n";
     std::cerr << "ERROR: Must provide a sample id or name\n";
     exit_with_error = true;
  }

  // barcode length must be positive
  if (options.barcode_length <= 0) {
     std::cout << "ERROR: Barcode length must be a positive integer\n";
     std::cerr << "ERROR: Barcode length must be a positive integer\n";
     exit_with_error = true;
  }

  // UMI length must be positive
  if (options.umi_length <= 0) {
     std::cout << "ERROR: UMI length must be a positive integer\n";
     std::cerr << "ERROR: UMI length must be a positive integer\n";
     exit_with_error = true;
  }

  // just prints out the files
  if (verbose_flag) {
      if (options.I1s.size()) {
          _print_file_info(options.I1s, std::string("I1"));
      }

      if (options.R1s.size()) {
          _print_file_info(options.R1s, std::string("R1"));
      }

      if (options.R2s.size()) {
          _print_file_info(options.R2s, std::string("R2"));
      }
  }

  if (exit_with_error) {
     exit(1);
  }


}


/** @copydoc  _print_file_info */
void _print_file_info(const std::vector<std::string> &fastqs,
     const std::string &type) {
    if (fastqs.size()) {
        std::cout << "INFO " << type << " files:" << std::endl;
            for (int i= 0; i < fastqs.size(); i++) {
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
int64_t get_num_blocks(const INPUT_OPTIONS &options) {
    double tot_size = 0;
    for (int i= 0; i < options.R1s.size(); i++) {
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

