/**
 *  @file   inputoptions.cpp
 *  @brief  functions for optons and input checking
 *  @author Kishori Konwar
 *  @date   2021-07-07
 ***********************************************/
#include "inputoptions.h"
#include <experimental/filesystem>

using namespace std;
namespace fs = std::experimental::filesystem;

/** @copydoc read_options_tagsort */
void read_options_tagsort(int argc, char **argv, INPUT_OPTIONS_TAGSORT &options) {
  int c;
  int i;

  int verbose_flag = 0;

  static struct option long_options[] = {
          /* These options set a flag. */
          {"verbose",                   no_argument,       0, 'v'},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"bam-input",                  required_argument, 0, 'b'},
          {"temp-folder",                required_argument, 0, 't'},
          {"output",                     required_argument, 0, 'o'},
          {"inmemory-chunk-size",        required_argument, 0, 'p'},
          {"tags",                       required_argument, 0, 'T'},
          {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char *help_messages[] = {
           "verbose messages  ",
           "input bam file [required]",
           "temp folder for disk sorting [options: default /tmp]",
           "output file [required]",
           "size of chunks for in-memory sorting [optional: default 1 GB]",
           "tags to sort by [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "b:t:o:p:T:v",
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
            options.bam_input = string(optarg);
            break;
        case 't':
            options.temp_folder = string(optarg);
            break;
        case 'o':
            options.output_file = string(optarg);
            break;
        case 'p':
            options.inmemory_chunk_size = atof(optarg);
            break;
        case 'T':
            options.tags.push_back(string(optarg));
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
  // bam input file must be there
  bool exit_with_error = false;
  if (options.bam_input.size() == 0) {
     std::cout << "ERROR: Must specify a input file name " << std::endl;
     std::cerr << "ERROR: Must specify a input file name " << std::endl;
     exit_with_error = true;
     exit(1);
  }

  // output file must exist
  if (options.output_file.size() == 0) {
     std::cout << "ERROR: Must specify a output file name " << std::endl;
     std::cerr << "ERROR: Must specify a output file name " << std::endl;
     exit_with_error = true;
     exit(1);
  }

  // check for input file
  if (not fs::exists(options.bam_input.c_str())) {
      std::cout << "ERROR " << "bam_input" << options.bam_input << " is missing!\n";
      std::cerr << "ERROR " << "bam_input" << options.bam_input << " is missing!\n";
      exit_with_error = true;
      exit(1);
  }

  // check for the temp folder
  if (not fs::exists(options.temp_folder.c_str())) {
      std::cout << "ERROR " << "temp folder " << options.temp_folder <<  " is missing!\n";
      std::cerr << "ERROR " << "temp folder " << options.temp_folder <<  " is missing!\n";
      exit_with_error = true;
      exit(1);
  }

  // The size of a set of aligments for in-memory sorting must be positive
  if (options.inmemory_chunk_size <= 0) {
     std::cout << "ERROR: The size (in GB) of chunks for in-memory sorting must be positive\n";
     std::cerr << "ERROR: The size (in GB) of chunks for in-memory sorting must be positive\n";
     exit_with_error = true;
      exit(1);
  }

  if (exit_with_error) {
     exit(1);
  }
}


/** @copydoc read_options_fastqprocess */
void read_options_fastqprocess(int argc, char **argv, INPUT_OPTIONS_FASTQPROCESS &options) {
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

