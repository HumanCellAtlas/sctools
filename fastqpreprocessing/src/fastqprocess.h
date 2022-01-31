#ifndef __FASTQ_PROCESS_H__
#define __FASTQ_PROCESS_H__
/**
 *  @file   fastqprocess.h
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/
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

#include "utilities.h"
#include "input_options.h"


/// Samrecord bins to be accessed by all threads
typedef struct SamRecordBins {
    /// array or array of samrecords
    /// one array for each reader samrecords[r].
    /// Note that every samrecord[r][i] for 0 <= i < num_records[r]
    /// is destined to be written to one specfic output bam file.
    /// These index of the bam file is one of the vectors "file_index[r][b]",
    /// where b is one of the bam files index
    SamRecord **samrecords;

    /// number of records in individual reader threads that
    /// can be written, i.e., num_records[r] stores the number of
    /// such records in samrecrods[r]
    int32_t *num_records;

    /// An array of arrays of vector, one array "file_index[r]" for reader thread r.
    /// The value vector in file_index[r][i] stores the indices of the samrecords[r]
    /// where the record in samrecords[r][i], for 0 <= i < num_records[r]
    /// should be written. This information is used by the writer threads.
    vector<int32_t> **file_index;

    /// sample name
    std::string sample_id;
    int32_t block_size;

    /// number of output bam files, and one writer thread per bam file
    int16_t num_files;
    /// flag to stop the writer
    bool stop;
    /// the thread (reader) that is currently wanting to write
    int32_t active_thread_num;
} SAM_RECORD_BINS;


/**
 * @brief Processes the input fastq files
 *
 * @detail
 *  This function creates a set of readers (as many as there are files),
 * a set of writers to write the individual bam files, a set of
 * semaphores for readers to signal to writers when buffer of records
 * are ready, and another set of semaphores for writers to signal
 * readers then the buffer has been emptied out and, therefore, reader
 * can go ahead and fill with more records.
 *
 * @params options user input options
 * @params white_list_data data-structure to store barcode correction
 *         map and vector of correct barcodes
*/
void process_inputs(const INPUT_OPTIONS_FASTQPROCESS & options, \
                   const WHITE_LIST_DATA * white_list_data);

/**
 * @brief Process one triplet of file R1/R2 and I1 in a thread
 *
 * @detail
 *   This function will be run by a thread for each set of R1/R2 and I1
 * files.
 *
 * @param tindex reader thread index
 * @param filename name of I1 file
 * @param filename1 name of R1 file
 * @param filename2 name of R2 file
 * @param barcode_length length of a barcode
 * @param umi_length length of UMI
 * @param white_list_data  data-structure barcode-correction based on
 *                         white list
 * @param samrecord_bins  bins for samrecords from the reader threads
*/

void process_file(int32_t tindex, std::string filename, String filename1, \
                  String filename2,  unsigned int barcode_length, \
                  unsigned int umi_length, \
                  const WHITE_LIST_DATA *white_list_data, \
                  SAM_RECORD_BINS * samrecord_bins);

/**
 * @brief Function for the writer thread
 *
 * @detail
 *  Dependeing on the number of output bam files there are as many
 * writer thread as there are output bam files. Each writer thread
 * writers into only one bam file
 *
 * @param  windex  index of the writer thread
 * @param samrecord_bins  bins for samrecords from the reader threads
*/
void bam_writers(int32_t windex, SAM_RECORD_BINS *samrecord_bins);

/**
 * @brief Function for the writer thread
 *
 * @detail
 *  Dependeing on the number of output bam files there are as many
 * writer thread as there are output bam files. Each writer thread
 * writers into only one bam file
 *
 * @param  windex  index of the writer thread
 * @param samrecord_bins  bins for samrecords from the reader threads
*/
void fastq_writers(int32_t windex, SAM_RECORD_BINS *samrecord_bins);
#endif
