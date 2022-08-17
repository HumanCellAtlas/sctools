/**
 *  @file   fastqprocess.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/

#include "input_options.h"
#include "utilities.h"

#include <gzstream.h>
#include <iostream>
#include <fstream>

#include <cstdint>

#include "FastQFile.h"
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

// TODO DEDUP
/// Samrecord bins to be accessed by all threads
struct SAM_RECORD_BINS
{
  /// array or array of samrecords
  /// one array for each reader samrecords[r].
  /// Note that every samrecord[r][i] for 0 <= i < num_records[r]
  /// is destined to be written to one specfic output bam file.
  /// These index of the bam file is one of the vectors "file_index[r][b]",
  /// where b is one of the bam files index
  SamRecord** samrecords;

  /// number of records in individual reader threads that
  /// can be written, i.e., num_records[r] stores the number of
  /// such records in samrecrods[r]
  int32_t* num_records;

  /// An array of arrays of vector, one array "file_index[r]" for reader thread r.
  /// The value vector in file_index[r][i] stores the indices of the samrecords[r]
  /// where the record in samrecords[r][i], for 0 <= i < num_records[r]
  /// should be written. This information is used by the writer threads.
  std::vector<int32_t>** file_index;

  /// sample name
  std::string sample_id;
  int32_t block_size;

  /// number of output bam files, and one writer thread per bam file
  int16_t num_files;
  /// flag to stop the writer
  bool stop;
  /// the thread (reader) that is currently wanting to write
  int32_t active_thread_num;
};

/// number of samrecords per buffer in each reader
constexpr int kSamRecordBufferSize = 100000;

/// array of semaphores for readers
sem_t* g_semaphores = 0;

/// array of semaphores for writers
sem_t* g_semaphores_workers = 0;

/** @copydoc create_record_holders */
SAM_RECORD_BINS* create_samrecord_holders(int16_t nthreads, const std::string sample_id,
                                          int16_t num_files)
{
  // samrecord data to hold buffer for the reader
  SAM_RECORD_BINS* samrecord_data = new SAM_RECORD_BINS;
  if ((samrecord_data->samrecords = new SamRecord *[nthreads]) == 0)
    crash("Failed to allocate memory for the samRecords pointer arrays");

  // one samrecord per reader thread to re-use repeatedly while writing
  for (int i = 0; i < nthreads; i++)
    if ((samrecord_data->samrecords[i] = new SamRecord[kSamRecordBufferSize]) == 0)
      crash("Failed to allocate memory for the samRecords");

  // for each reader thread keep the number of records to write out
  if ((samrecord_data->num_records = new int[nthreads]) == 0)
    crash("Failed to allocate memory for the num records array");

  // for each thread  we allocate an array of indices (to final output files)
  if ((samrecord_data->file_index = new std::vector<int>* [nthreads]) == 0)
    crash("Failed to allocate memory for the pointer for array of vectors");

  // for each read thread allocate the index vector
  for (int i = 0; i < nthreads; i++)
    if ((samrecord_data->file_index[i] = new std::vector<int>[num_files]) == 0)
      crash("Failed to allocate memory for the vectors for index of file");

  // set the remaining data
  samrecord_data->sample_id = sample_id;
  samrecord_data->num_files = num_files;
  samrecord_data->stop = false;
  return samrecord_data;
}

void process_file(int tindex, std::string filenameI1, String filenameR1,
                  String filenameR2,  unsigned int barcode_length,
                  unsigned int umi_length,
                  const WHITE_LIST_DATA* white_list_data,
                  SAM_RECORD_BINS* samrecord_data);
void fastq_writers(int windex, SAM_RECORD_BINS* samrecord_data);
void bam_writers(int windex, SAM_RECORD_BINS* samrecord_data);

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
void process_inputs(InputOptionsFastqProcess const& options,
                    const WHITE_LIST_DATA* white_list_data)
{
  // number of files based on the input size
  int num_files = getNumBlocks(options);
  // create the data for the threads
  SAM_RECORD_BINS* samrecord_data =
    create_samrecord_holders(options.R1s.size(), options.sample_id, num_files);

  g_semaphores_workers = new sem_t[num_files];
  for (int i = 0; i < num_files; i++)
    sem_init((g_semaphores_workers + i), 0, 0);

  // create the bam file writers semaphores
  g_semaphores = new sem_t[num_files];
  for (int i = 0; i < num_files; i++)
    sem_init((g_semaphores + i), 0, 0);

  // execute the bam file writers threads
  std::vector<std::thread> writers;
  if (options.output_format=="BAM")
    for (int i = 0; i < num_files; i++)
      writers.emplace_back(bam_writers, i, samrecord_data);
  else if (options.output_format=="FASTQ")
    for (int i = 0; i < num_files; i++)
      writers.emplace_back(fastq_writers, i, samrecord_data);
  else
    crash("ERROR: Output-format must be either FASTQ or BAM");

  // execute the fastq readers threads
  std::vector<std::thread> readers;
  for (unsigned int i = 0; i < options.R1s.size(); i++)
  {
    assert(options.I1s.empty() || options.I1s.size() == options.R1s.size());
    // if there is no I1 file then send an empty file name
    std::string I1 = options.I1s.empty() ? "" : options.I1s[i];

    readers.emplace_back(process_file, i, I1, options.R1s[i].c_str(),
                         options.R2s[i].c_str(), options.barcode_length,
                         options.umi_length, white_list_data, samrecord_data);
  }

  // every reader thread joins.
  for (auto& reader : readers)
    reader.join();

  // set the stop flag for the writers
  samrecord_data->stop = true;

  // ask the writers to make one more loop in the while loop
  for (int j = 0; j < samrecord_data->num_files; j++)
    if (sem_post(&g_semaphores[j]) == -1)
      crashWithPerror("sem_post: g_semaphores");

  // wait for the writers to stop after they have seen the stop flag
  for (auto& writer : writers)
    writer.join();

  // destroy the g_semaphores
  for (int i = 0; i < samrecord_data->num_files; i++)
    sem_destroy(&g_semaphores[i]);

  // destroy the g_semaphores for g_semaphores_workers
  for (int i = 0; i < samrecord_data->num_files; i++)
    sem_destroy(&g_semaphores_workers[i]);

  // delete the records
  delete [] samrecord_data->num_records;
}

void fastq_writers(int windex, SAM_RECORD_BINS* samrecord_data)
{
  std::string r1_output_fname = "fastq_R1_" + std::to_string(windex) + ".fastq.gz";
  ogzstream r1_out(r1_output_fname.c_str());
  if (!r1_out)
    crash("ERROR: Failed to open R1 fastq file " + r1_output_fname + " for writing");

  std::string r2_output_fname = "fastq_R2_" + std::to_string(windex) + ".fastq.gz";
  ogzstream r2_out(r2_output_fname.c_str());
  if (!r2_out)
    crash("ERROR: Failed to open R2 fastq file " + r2_output_fname + " for writing");

  // keep writing forever, until there is a flag to stop
  while (true)
  {
    // wait until some data is ready from a reader thread
    if (sem_wait(&g_semaphores[windex]) == -1)
      crashWithPerror("sem_wait:g_semaphores");

    // write out the record buffers for the reader thread "active_thread_num"
    // that signalled that buffer is ready to be written
    SamRecord* samRecord = samrecord_data->samrecords[samrecord_data->active_thread_num];

    // go through the index of the samrecords that are stored for the current
    // writer, i.e., "windex" or the corresponding BAM file
    for (auto index : samrecord_data->file_index[samrecord_data->active_thread_num][windex])
    {
      //       samOut.WriteRecord(samHeader, samRecord[index]);
      r1_out << "@" << samRecord[index].getReadName() << std::endl
             << samRecord[index].getString("CR").c_str() << samRecord[index].getString("UR") << std::endl
             << "+" << std::endl
             << samRecord[index].getString("CY") << samRecord[index].getString("UY") << std::endl;
    }

    for (auto index : samrecord_data->file_index[samrecord_data->active_thread_num][windex])
    {
      //       samOut.WriteRecord(samHeader, samRecord[index]);
      r2_out << "@" << samRecord[index].getReadName() << std::endl
             << samRecord[index].getSequence() << std::endl
             << "+" << std::endl
             << samRecord[index].getQuality() << std::endl;
    }

    // lets the reads thread know that I am done writing the
    // buffer that are destined to be my file
    if (sem_post(&g_semaphores_workers[windex]) == -1)
      crashWithPerror("sem_post: g_semaphores_workers");

    // time to stop variable is valid
    if (samrecord_data->stop)
      break;
  }

  // close the fastq files
  r1_out.close();
  r2_out.close();
}

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
void bam_writers(int windex, SAM_RECORD_BINS* samrecord_data)
{
  std::string bam_out_fname = "subfile_" + std::to_string(windex) + ".bam";
  SamFile samOut;
  samOut.OpenForWrite(bam_out_fname.c_str());

  // Write the sam header.
  SamFileHeader samHeader;

  // add the HD tags for the header
  samHeader.setHDTag("VN", "1.6");
  samHeader.setHDTag("SO", "unsorted");

  // add the RG group tags
  SamHeaderRG* headerRG = new SamHeaderRG;
  headerRG->setTag("ID", "A");
  headerRG->setTag("SM", samrecord_data->sample_id.c_str());
  samHeader.addRG(headerRG);

  // add the header to the output bam
  samOut.WriteHeader(samHeader);

  // keep writing forever, until there is a flag to stop
  while (true)
  {
    // wait until some data is ready from a reader thread
    if (sem_wait(&g_semaphores[windex]) == -1)
      crashWithPerror("sem_wait:g_semaphores");

    // write out the record buffers for the reader thread "active_thread_num"
    // that signalled that buffer is ready to be written
    SamRecord* samRecord  =
      samrecord_data->samrecords[samrecord_data->active_thread_num];
    // go through the index of the samrecords that are stored for the current
    // writer, i.e., "windex" or the corresponding BAM file
    for (auto index : samrecord_data->file_index[samrecord_data->active_thread_num][windex])
      samOut.WriteRecord(samHeader, samRecord[index]);

    // lets the reads thread know that I am done writing the
    // buffer that are destined to be my file
    if (sem_post(&g_semaphores_workers[windex]) == -1)
      crashWithPerror("sem_post: g_semaphores_workers");

    // time to stop variable is valid
    if (samrecord_data->stop)
      break;
  }

  // close the bamfile
  samOut.Close();
}

/**
 * @brief fillSamRecord fill a SamRecord with the sequence and TAGs data
 *
 * @param samRecord  the SamRecord to fill with the data
 * @param fastQFileI1  the I1 fastq file
 * @param fastQFileR1  the R1 fastq file
 * @param fastQFileR2  the R2 fastq file
 * @param samRecord  the SamRecord to fill with the data
 * @param barcode_length length of UMI
 * @param has_I1_file_list a boolean indicating if I1 files are avaiable
*/
void fillSamRecord(SamRecord* samRecord, FastQFile& fastQFileI1,
                   FastQFile& fastQFileR1, FastQFile& fastQFileR2,
                   unsigned int barcode_length, unsigned int umi_length,
                   bool has_I1_file_list)
{
  // check the sequence names matching
  std::string a = std::string(fastQFileR1.myRawSequence.c_str());
  std::string b = std::string(fastQFileR1.myQualityString.c_str());

  // extract the raw barcode and UMI
  std::string barcode = a.substr(0, barcode_length);
  std::string UMI  = a.substr(barcode_length, umi_length);

  // extract raw barcode and UMI quality string
  std::string barcodeQString = b.substr(0, barcode_length);
  std::string UMIQString  = b.substr(barcode_length, umi_length);

  // reset the samrecord
  samRecord->resetRecord();

  // add read group and the sam flag
  samRecord->addTag("RG", 'Z', "A");
  samRecord->setFlag(4);

  // add idenfifier, sequence and quality score of the alignments
  samRecord->setReadName(fastQFileR2.mySequenceIdentifier.c_str());
  samRecord->setSequence(fastQFileR2.myRawSequence.c_str());
  samRecord->setQuality(fastQFileR2.myQualityString.c_str());

  // add barcode and quality
  samRecord->addTag("CR", 'Z', barcode.c_str());
  samRecord->addTag("CY", 'Z', barcodeQString.c_str());

  // add UMI
  samRecord->addTag("UR", 'Z', UMI.c_str());
  samRecord->addTag("UY", 'Z', UMIQString.c_str());

  // add raw sequence and quality sequence for the index
  if (has_I1_file_list)
  {
    std::string indexseq = std::string(fastQFileI1.myRawSequence.c_str());
    std::string indexSeqQual = std::string(fastQFileI1.myQualityString.c_str());
    samRecord->addTag("SR", 'Z', indexseq.c_str());
    samRecord->addTag("SY", 'Z', indexSeqQual.c_str());
  }
}
/**
   @brief getBukcetIndex computes the index for the bucket (of bam file)
 *    for a barcode and also add the correct barcode to the SamRecord
 *
 * @param barcode the barcode for which a bucket is computed
 * @param samRecord  the partially filled samrecord to add the corrected barcode
 * @param white_list_data the white list data for barcode correction
 * @param samrecord_data shared data for the various readers
 * @param n_barcode_corrected a variable keeping track of the number of barcodes corrected
 * @param n_barcode_correct a the number of barcodes so far are already correct
 * @param n_barcode_errrosv keeping track of the number of barcodes that are incorrectible
 *
 * @return the bucket number where the current SamRecord should go to
*/
int32_t getBucketIndex(const std::string& barcode, SamRecord* samRecord,
                       const WHITE_LIST_DATA* white_list_data, SAM_RECORD_BINS* samrecord_data,
                       int* n_barcode_corrected, int* n_barcode_correct, int* n_barcode_errors)
{

  std::string correct_barcode;
  // bucket barcode is used to pick the target bam file
  // This is done because in the case of incorrigible barcodes
  // we need a mechanism to uniformly distribute the alignments
  // so that no bam is oversized to putting all such barcode less
  // sequences into one particular. Incorregible barcodes are simply
  // added withouth the CB tag
  std::string bucket_barcode;
  if (white_list_data->mutations.find(barcode) != white_list_data->mutations.end())
  {
    if (white_list_data->mutations.at(barcode) == -1) // -1 means raw barcode is correct
    {
      correct_barcode = barcode;
      *n_barcode_correct += 1;
    }
    else
    {
      // it is a 1-mutation of some whitelist barcode so get the
      // barcode by indexing into the vector of whitelist barcodes
      correct_barcode =
        white_list_data->barcodes.at(white_list_data->mutations.at(barcode));
      *n_barcode_corrected += 1;
    }
    // is used for computing the file index
    bucket_barcode = correct_barcode;

    // corrected barcode should be added to the samrecord
    samRecord->addTag("CB", 'Z', correct_barcode.c_str());
  }
  else     // not possible to correct the raw barcode
  {
    *n_barcode_errors += 1;
    bucket_barcode = barcode;
  }
  // destination bam file index computed based on the bucket_barcode
  int32_t bucket = std::hash<std::string> {}(bucket_barcode.c_str()) %
                   samrecord_data->num_files;

  return bucket;
}

std::mutex g_block_mutex;

/**
 * @brief submit_block_tobe_written this function is for a reader to send
 *    signal to the writer to empty the current list of SamRecords that are
 *    ready to be written out
 *
 * @param samrecord_data  the samrecord data
 * @param tindex the index of the thread
*/
void submit_block_tobe_written(SAM_RECORD_BINS* samrecord_data, int tindex)
{
  g_block_mutex.lock();

  // it sets itself as the active thread who wants the
  // readers to clear the data
  samrecord_data->active_thread_num = tindex;

  // send a signal to every writer thread, i.e., write out data
  // data to any file where the samheader should be written to
  for (int32_t j = 0; j < samrecord_data->num_files; j++)
    if (sem_post(&g_semaphores[j]) == -1)
      crashWithPerror("sem_post: g_semaphores");

  // there is where I wait while the writers are writing
  for (int32_t j = 0; j < samrecord_data->num_files; j++)
    if (sem_wait(&g_semaphores_workers[j]) == -1)
      crashWithPerror("sem_wait: g_semaphores_workers");

  // they are done writing
  for (int j = 0; j < samrecord_data->num_files; j++)
    samrecord_data->file_index[tindex][j].clear();

  // not records to write in the current
  samrecord_data->num_records[tindex] = 0;

  // release the mutex, so that other readers might want to write
  g_block_mutex.unlock();
}

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
void process_file(int tindex, std::string filenameI1, String filenameR1,
                  String filenameR2,  unsigned int barcode_length,
                  unsigned int umi_length,
                  const WHITE_LIST_DATA* white_list_data,
                  SAM_RECORD_BINS* samrecord_data)
{
  /// setting the shortest sequence allowed to be read
  FastQFile fastQFileI1(4, 4);
  FastQFile fastQFileR1(4, 4);
  FastQFile fastQFileR2(4, 4);

  // open the I1 file
  bool has_I1_file_list = true;
  if (!filenameI1.empty())
  {
    if (fastQFileI1.openFile(String(filenameI1.c_str()), BaseAsciiMap::UNKNOWN) !=
        FastQStatus::FASTQ_SUCCESS)
    {
      std::cerr << "Failed to open file: " <<  filenameI1.c_str();
      abort();
    }
  }
  else
    has_I1_file_list = false;

  // open the R1 file
  if (fastQFileR1.openFile(filenameR1, BaseAsciiMap::UNKNOWN) !=
      FastQStatus::FASTQ_SUCCESS)
  {
    std::cerr << "Failed to open file: " <<  filenameR1.c_str();
    abort();
  }

  // open the R2 file
  if (fastQFileR2.openFile(filenameR2, BaseAsciiMap::UNKNOWN) !=
      FastQStatus::FASTQ_SUCCESS)
  {
    std::cerr << "Failed to open file: " <<  filenameR2.c_str();
    abort();
  }

  // point to the array of records already allocated for this reader
  SamRecord* samRecord  = samrecord_data->samrecords[tindex];

  // Keep reading the file until there are no more fastq sequences to process.
  int i = 0;
  int n_barcode_errors = 0;
  int n_barcode_corrected = 0;
  int n_barcode_correct = 0;
  int r = 0;
  printf("Opening the thread in %d\n", tindex);


  while (fastQFileR1.keepReadingFile())
  {
    if ((!has_I1_file_list ||
         (
           has_I1_file_list &&
           fastQFileI1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
         )
        ) &&
        fastQFileR1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS &&
        fastQFileR2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS)
    {
      i = i + 1;

      // prepare the rth samrecord with the sequence, sequence quality
      SamRecord* samrec = samRecord + r;
      // barcode and UMI and their quality sequences
      fillSamRecord(samrec, fastQFileI1, fastQFileR1, fastQFileR2,
                    barcode_length, umi_length, has_I1_file_list);

      // extract the raw barcode and UMI
      std::string a = std::string(fastQFileR1.myRawSequence.c_str());
      std::string barcode = a.substr(0, barcode_length);

      // bucket barcode is used to pick the target bam file
      // This is done because in the case of incorrigible barcodes
      // we need a mechanism to uniformly distribute the alignments
      // so that no bam is oversized to putting all such barcode less
      // sequences into one particular. Incorregible barcodes are simply
      // added withouth the CB tag
      int32_t bucket =  getBucketIndex(barcode, samrec, white_list_data,
                                       samrecord_data, &n_barcode_corrected,
                                       &n_barcode_correct, &n_barcode_errors);

      samrecord_data->num_records[tindex]++;
      // store the index of the record in the right vector that is done
      // to serve as a bucket B to hold the indices of samrecords that are
      // going to bamfile B
      samrecord_data->file_index[tindex][bucket].push_back(r);

      samrecord_data->num_records[tindex]++;
      // write a block of samrecords
      r = r + 1;

      // Once kSamRecordBufferSize samrecords are read, or there
      // is no more sequence to be read from the file then it is time to
      // signal the writer threads to clear the buffer to the bam files
      // only one reader should be successfull in doing so.

      // However, if a block of data is not read or end of the FASTQ file is
      // seen then continue reading the FASTQ files and keep creating the
      // sam records in its buffer. This is the same behavior across all
      // reader threads
      if (r == kSamRecordBufferSize || !fastQFileR1.keepReadingFile())
      {
        submit_block_tobe_written(samrecord_data, tindex);

        // start reading a new block of FASTQ sequences
        r = 0;
      }
      if (i % 10000000 == 0)
      {
        printf("%d\n", i);
        std::string a = std::string(fastQFileR1.myRawSequence.c_str());
        printf("%s\n", fastQFileR1.mySequenceIdLine.c_str());
        printf("%s\n", fastQFileR2.mySequenceIdLine.c_str());
      }
    }  //  if successful read of a sequence
  }

  // Finished processing all of the sequences in the file.
  // Close the input files.
  if (has_I1_file_list)
    fastQFileI1.closeFile();
  fastQFileR1.closeFile();
  fastQFileR2.closeFile();
  printf("Total barcodes:%d\n correct:%d\ncorrected:%d\nuncorrectible"
         ":%d\nuncorrected:%lf\n",
         i, n_barcode_correct, n_barcode_corrected, n_barcode_errors,
         n_barcode_errors/static_cast<double>(i) *100);
}


/* Flag set by ‘--verbose’. */
int main(int argc, char** argv)
{
  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);

  std::cout << "reading whitelist file " << options.white_list_file << "...";
  std::unique_ptr<WHITE_LIST_DATA> white_list_data = readWhiteList(options.white_list_file);
  std::cout << "done" << std::endl;

  process_inputs(options, white_list_data.get());
  return 0;
}

