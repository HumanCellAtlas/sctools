#ifndef __FASTQ_PROCESS_H__
#define __FASTQ_PROCESS_H__

#include <FastQFile.h>
#include "FastQStatus.h"
#include "BaseAsciiMap.h"
#include "SamFile.h"
#include "SamValidation.h"

#include <semaphore.h>
#include <thread>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <vector>
#include <functional>
#include <mutex>         


using namespace std;
typedef std::pair<string, bool>  STRING_BOOL_PAIR;

typedef vector<string>  STRING_VECTOR;

typedef std::unordered_map <std::string, int32_t> STRING_BOOL_MAP;

typedef struct _white_list_data {
    STRING_BOOL_MAP mutations;
    STRING_VECTOR barcodes;
} WHITE_LIST_DATA;
 
typedef struct SamRecordBins {
    SamRecord **samrecords;
    int32_t *num_records;
    vector<int32_t> **file_index;
    string sample_id;
    int32_t block_size;
    int16_t num_files;
    bool stop;
    int16_t num_threads;
    int32_t active_thread_no;
} SAM_RECORD_BINS;

void process_file(int32_t tindex, String filename, String filename1, String filename2, \
                   unsigned int barcode_length, unsigned int umi_length, \
                   const WHITE_LIST_DATA *, SAM_RECORD_BINS *) ;

void bam_writers(int32_t , SAM_RECORD_BINS *) ;
#endif
