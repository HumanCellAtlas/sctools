/** 
 *  @file   fastqprocess.cpp 
 *  @brief  functions for file processing 
 *  @author Kishori Konwar 
 *  @date   2020-08-27 
 ***********************************************/

#include "fastqprocess.h"
#include "utilities.h"
#include <cstdint>

/// maximum file length
#define MAX_FILE_LENGTH 500

/// number of samrecords per buffer in each reader
#define SAMRECORD_BUFFER_SIZE 100000

/// mutex 
std::mutex mtx;

/// array of semaphores for readers
sem_t *semaphores = 0;

/// array of semaphores for rwriters
sem_t *semaphores_workers = 0;

/** @copydoc create_record_holders */
SAM_RECORD_BINS * create_samrecord_holders(int16_t nthreads, 
                                           int32_t block_size, 
                                           const std::string sample_id, 
                                           int16_t num_files) {

    // samrecord data to hold buffer for the reader
    SAM_RECORD_BINS *samrecord_data = new SAM_RECORD_BINS;
    if ((samrecord_data->samrecords = new SamRecord *[nthreads]) == 0) {
        std::cerr << "Failed to allocate memory for the " 
                     "samRecords pointer arrays" << std::endl;
        return 0;
    }
    
    // one samrecord per reader thread to re-use repeatedly while writing
    for (int i = 0; i < nthreads; i++) {
       if ((samrecord_data->samrecords[i] = new SamRecord[block_size]) == 0) {
          std::cerr << "Failed to allocate memory for the " 
                       "samRecords" << std::endl;
          return 0;
         }
    }

    // for each reader thread keep the number of records to write out
    if ((samrecord_data->num_records = new int[nthreads]) == 0) {
       std::cerr << "Failed to allocate memory for the num " 
                    "records array" << std::endl;
       return 0;
    }

    // for each thread  we allocate an array of indices (to final output files(
    if ((samrecord_data->file_index = new vector<int> *[nthreads]) == 0) {
       std::cerr << "Failed to allocate memory for the pointer for " 
                    "array of vectors" << std::endl;
       return 0;
    }

    // for each read thread allocate the index vector
    for (int i = 0; i < nthreads; i++) {
      if ((samrecord_data->file_index[i] = new vector<int>[num_files]) == 0) {
         std::cerr << "Failed to allocate memory for the vectors for " 
                      "index of file" << std::endl;
         return 0;
      }
    }

    // set the remaining data
    samrecord_data->block_size = block_size;
    samrecord_data->sample_id = sample_id;
    samrecord_data->num_files = num_files;
    samrecord_data->stop = false;
    return samrecord_data;
}

/** @copydoc process_inputs */
void process_inputs(const INPUT_OPTIONS &options, 
                   const WHITE_LIST_DATA *white_list_data) {
     
     int block_size = SAMRECORD_BUFFER_SIZE;

     // number of files based on the input size  
     int num_files = get_num_blocks(options);
     // create the data for the threads
     SAM_RECORD_BINS *samrecord_data = 
               create_samrecord_holders(options.R1s.size(), block_size, 
                                        options.sample_id, num_files);
     semaphores_workers = new sem_t[num_files];
     for (int i = 0; i < num_files; i++) {
        sem_init((semaphores_workers + i), 0, 0);
     }

     // create the bam file writers semaphores
     semaphores = new sem_t[num_files];
     for (int i = 0; i < num_files; i++) {
        sem_init((semaphores + i), 0, 0);
     }

     // execute the bam file writers threads
     std::thread *writers = new std::thread[num_files];
     for (int i = 0; i < num_files; i++) {
        writers[i] = std::thread(bam_writers, i, samrecord_data);
     }

     // execute the fastq readers threads
     std::thread  *readers = new std::thread[options.R1s.size()];
     for (int i = 0; i < options.R1s.size(); i++) {
        std::string  I1;
        // if there is no I1 file then send an empty file name
        if (options.I1s.size() > 0) {
           I1 = std::string(options.I1s[i].c_str());
        } else {
           I1 = std::string("");
        }
        readers[i] = std::thread(process_file, i, I1, 
                       options.R1s[i].c_str(), options.R2s[i].c_str(), 
                       options.barcode_length, options.umi_length,
                       white_list_data, samrecord_data);
     }

     // every reader thread joins.
     for (int i = 0; i < options.R1s.size(); i++) {
        readers[i].join();
     }
     // set the stop flag for the writers
     samrecord_data->stop = true;

     // ask the writers to make one more loop in the whilte loop
     for (int j = 0; j < samrecord_data->num_files; j++) {
        if (sem_post(&semaphores[j]) == -1)
           error("sem_post: semaphores");
     }

     // wait for the writers to stop after they have seen the stop flag
     for (int i = 0; i < samrecord_data->num_files; i++) {
        writers[i].join();
     }

     // destroy the semaphores
     for (int i = 0; i < samrecord_data->num_files; i++) {
         sem_destroy(&semaphores[i]);
     }

     // destroy the semaphores for semaphores_workers
     for (int i = 0; i < samrecord_data->num_files; i++) {
         sem_destroy(&semaphores_workers[i]);
     }
     
     // delete reader and writer threads
     delete [] readers;
     delete [] writers;
}

/** @copydoc bam_writers */
void bam_writers(int windex, SAM_RECORD_BINS *samrecord_data) {
    SamFile samOut;
    std::string outputfile;
    
    // name of the output file
    char buf[MAX_FILE_LENGTH];
    sprintf(buf, "subfile_%d.bam", windex);
    outputfile = buf;

    // open to write the outputfile
    samOut.OpenForWrite(outputfile.c_str());

    // Write the sam header.
    SamFileHeader samHeader;

    // add the HD tags for the header
    samHeader.setHDTag("VN", "1.6");
    samHeader.setHDTag("SO", "unsorted");

    // add the RG group tags
    SamHeaderRG *headerRG = new SamHeaderRG;
    headerRG->setTag("ID", "A");
    headerRG->setTag("SM", samrecord_data->sample_id.c_str());
    samHeader.addRG(headerRG);
   
    // add the header to the output bam
    samOut.WriteHeader(samHeader);

    // keep writing foever, until there is a flag to stop
    while (true) {
      // wait untile some data is ready from a reader thread
      if (sem_wait(&semaphores[windex]) == -1)
         error("sem_wait:semaphores");

      // write out the record buffers for the reader thread "active_thread_no"
      // that signalled that buffer is ready to be written
      SamRecord *samRecord  = 
         samrecord_data->samrecords[samrecord_data->active_thread_no];
      // go through the index of the samrecords that are stored for the current 
      // writer, i.e., "windex" or the corresponding BAM file
      for (auto index : samrecord_data->file_index[samrecord_data->active_thread_no][windex]) {
           samOut.WriteRecord(samHeader, samRecord[index]);
      }

      // lets the reads thread know that I am done writing the 
      // buffer that are destined to be my file
      if (sem_post(&semaphores_workers[windex]) == -1)
          error("sem_post: semaphores_workers");

      // time to stop variable is valid
      if (samrecord_data->stop) break;
    }

    // close the bamfile
    samOut.Close();
}


void process_file(int tindex, std::string filenameI1, String filenameR1, 
                   String filenameR2,  unsigned int barcode_length, 
                   unsigned int umi_length, 
                   const WHITE_LIST_DATA* white_list_data, 
                   SAM_RECORD_BINS *samrecord_data) {

    /// setting the shortest sequence allowed to be read
    FastQFile fastQFileI1(4, 4);
    FastQFile fastQFileR1(4, 4);
    FastQFile fastQFileR2(4, 4);

    // open the I1 file
    bool empty_I1_file_list = false;
    if (!filenameI1.empty()) {
        if (fastQFileI1.openFile(String(filenameI1.c_str()), BaseAsciiMap::UNKNOWN) != 
            FastQStatus::FASTQ_SUCCESS) {
            std::cerr << "Failed to open file: " <<  filenameI1.c_str();
            return;
        }
    } else {
        empty_I1_file_list = true;
    }

    // open the R1 file
    if (fastQFileR1.openFile(filenameR1, BaseAsciiMap::UNKNOWN) != 
        FastQStatus::FASTQ_SUCCESS) {
          std::cerr << "Failed to open file: " <<  filenameR1.c_str();
          return;
    }

    // open the R2 file
    if (fastQFileR2.openFile(filenameR2, BaseAsciiMap::UNKNOWN) != 
       FastQStatus::FASTQ_SUCCESS) {
          std::cerr << "Failed to open file: " <<  filenameR2.c_str();
          return;
    }

    // point to the array of records already allocated for this reader
    SamRecord *samRecord  = samrecord_data->samrecords[tindex];
    int block_size = samrecord_data->block_size;

    // Keep reading the file until there are no more fastq sequences to process.
    int i = 0;
    int n_barcode_errors = 0;
    int n_barcode_corrected = 0;
    int n_barcode_correct = 0;
    int r = 0;
    printf("Opening the thread in %d\n", tindex);


    while (fastQFileR1.keepReadingFile()) {
       if ((empty_I1_file_list == true || 
             (
             empty_I1_file_list == false && 
             fastQFileI1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
             )
           ) &&
          fastQFileR1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS &&
          fastQFileR2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS) {
          i = i + 1;

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
          samRecord[r].resetRecord();

          // add read group and the sam flag
          samRecord[r].addTag("RG", 'Z', "A");
          samRecord[r].setFlag(4);

          // add idenfifier, sequence and quality score of the alignments
          samRecord[r].setReadName(fastQFileR2.mySequenceIdentifier.c_str());
          samRecord[r].setSequence(fastQFileR2.myRawSequence.c_str());
          samRecord[r].setQuality(fastQFileR2.myQualityString.c_str());

          // add barcode and quality 
          samRecord[r].addTag("CR", 'Z', barcode.c_str());
          samRecord[r].addTag("CY", 'Z', barcodeQString.c_str());

          // add UMI
          samRecord[r].addTag("UR", 'Z', UMI.c_str());
          samRecord[r].addTag("UY", 'Z', UMIQString.c_str());

          // add raw sequence and quality sequence for the index 
          if (empty_I1_file_list == false) {
              std::string indexseq = std::string(fastQFileI1.myRawSequence.c_str());
              std::string indexSeqQual = std::string(fastQFileI1.myQualityString.c_str());
              samRecord[r].addTag("SR", 'Z', indexseq.c_str());
              samRecord[r].addTag("SY", 'Z', indexSeqQual.c_str());
          }

          std::string correct_barcode;
          // bucket barcode is used to pick the target bam file 
          // This is done because in the case of incorrigible barcodes
          // we need a mechanism to uniformly distribute the alignments
          // so that no bam is oversized to putting all such barcode less
          // sequences into one particular
          std::string bucket_barcode;
          if (white_list_data->mutations.find(barcode) != 
               white_list_data->mutations.end()) {

             // if -1 then the raw barcode is the correct barcode
             if (white_list_data->mutations.at(barcode) == -1) {
                  correct_barcode = barcode;
                  n_barcode_correct += 1;
             } else {
               // it is a 1-mutation of some whitelist barcode so get the 
               // barcode by indexing into the vector of whitelist barcodes 
               correct_barcode = 
               white_list_data->barcodes.at(white_list_data->mutations.at(barcode));
               n_barcode_corrected += 1;
             }
             // is used for computing the file index
             bucket_barcode = correct_barcode;
             samRecord[r].addTag("CB", 'Z', correct_barcode.c_str());
          } else { // not possible to correct the raw barcode
             n_barcode_errors += 1;
             bucket_barcode = barcode;
          }
          samrecord_data->num_records[tindex]++;
          // destination bam file index computed based on the bucket_barcode
          int32_t bucket = std::hash<std::string>{}(bucket_barcode.c_str()) % 
                       samrecord_data->num_files;
          // store the index of the record in the right vector that is done 
          // to serve as a bucket B to hold the indices of samrecords that are 
          // going to bamfile B 
          samrecord_data->file_index[tindex][bucket].push_back(r);
          // write a block of samrecords
          r = r + 1;
 
          // Once block_size amount of samrecord is read, or there 
          // is no more sequence to be read from the file then it is time to 
          // signal the writer threads to clear the buffer to the bam files 
          // only one reader should be successfull in doing so. 

          // However, if a block of data is not read or end of the FASTQ file is
          // seen then continue reading the FASTQ files and keep creating the 
          // sam records in its buffer. This is the same behavior across all 
          // reader threads
          if (r == block_size || !fastQFileR1.keepReadingFile()) {
                // mutex makes sure only one reader thread can lock
                mtx.lock();
              
                // it sets itself as the active thread who wants the 
                // readers to clear the data
                samrecord_data->active_thread_no = tindex;

                // send a signal to every writer thread, i.e., write out data 
                // data to any file where the samheader should be written to
                for (int32_t j = 0; j < samrecord_data->num_files; j++) {
                   if (sem_post(&semaphores[j]) == -1)
                      error("sem_post: semaphores");
                }

                // there is where I wait white the writers are writing
                for (int32_t j = 0; j < samrecord_data->num_files; j++) {
                   if (sem_wait(&semaphores_workers[j]) == -1)
                      error("sem_wait: semaphores_workers");
                }

                // they are done writing
                for (int j = 0; j < samrecord_data->num_files; j++) {
                      samrecord_data->file_index[tindex][j].clear();
                }
                // start reading a new block of FASTQ sequences 
                r = 0;
                // not records to write in the current  
                samrecord_data->num_records[tindex] = 0;
                 
                // release the mutex, so that other readers might want to write 
                mtx.unlock();
          }
          if (i % 10000000 == 0) {
              printf("%d\n", i);
              std::string a = std::string(fastQFileR1.myRawSequence.c_str());
              //printf("%s\n", fastQFileI1.mySequenceIdLine.c_str());
              printf("%s\n", fastQFileR1.mySequenceIdLine.c_str());
              printf("%s\n", fastQFileR2.mySequenceIdLine.c_str());
          }
         /*
          if(i== 5000000) {
              for(int j = 0; j < samrecord_data->num_files; j++) {
                  sem_destroy(&semaphores[j]);
              }
             break;
          }
          */
       }  //  if successful read of a sequence
    }

    // Finished processing all of the sequences in the file.
    // Close the input files.
    if(empty_I1_file_list == false) fastQFileI1.closeFile();

    fastQFileR1.closeFile();
    fastQFileR2.closeFile();
    printf("Total barcodes:%d\n correct:%d\ncorrected:%d\nuncorrectible"
           ":%d\nuncorrected:%lf\n", 
           i, n_barcode_correct, n_barcode_corrected, n_barcode_errors, 
           n_barcode_errors/static_cast<double>(i) *100);
}
