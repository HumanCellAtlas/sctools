/*

*/

#include "fastqprocess.h"
#include "utilities.h"
#include <cstdint>

std::mutex mtx;
sem_t *semaphores = 0;
sem_t *semaphores_workers = 0;

using namespace std;


SAM_RECORD_BINS * create_samrecord_holders(int16_t nthreads, \
                                           int32_t block_size, \
                                           const string sample_id, \
                                           int16_t num_files) {
    SAM_RECORD_BINS *samrecord_data = new SAM_RECORD_BINS;
    if ((samrecord_data->samrecords = new SamRecord *[nthreads]) == 0) {
        std::cerr << "Failed to allocate memory for the " \
                     "samRecords pointer arrays" << std::endl;
        return 0;
    }
    for (int i = 0; i < nthreads; i++) {
       if ((samrecord_data->samrecords[i] = new SamRecord[block_size]) == 0) {
          std::cerr << "Failed to allocate memory for the " \
                       "samRecords" << std::endl;
          return 0;
         }
    }
    if ((samrecord_data->num_records = new int[nthreads]) == 0) {
       std::cerr << "Failed to allocate memory for the num " \
                    "records array" << std::endl;
       return 0;
    }
    if ((samrecord_data->file_index = new vector<int> *[nthreads]) == 0) {
       std::cerr << "Failed to allocate memory for the pointer for " \
                    "array of vectors" << std::endl;
       return 0;
    }
    for (int i = 0; i < nthreads; i++) {
      if ((samrecord_data->file_index[i] = new vector<int>[num_files]) == 0) {
         std::cerr << "Failed to allocate memory for the vectors for " \
                      "index of file" << std::endl;
         return 0;
      }
    }
    samrecord_data->block_size = block_size;
    samrecord_data->sample_id = sample_id;
    samrecord_data->num_files = num_files;
    samrecord_data->stop = false;
    return samrecord_data;
}

int process_inputs(const INPUT_OPTIONS &options, \
                   const WHITE_LIST_DATA *white_list_data) {
     int num_files = 200;
     int block_size = 100000;
     num_files = get_num_blocks(options);
     // create the data for the threads
     SAM_RECORD_BINS *samrecord_data = \
               create_samrecord_holders(options.R1s.size(), block_size, \
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
        readers[i] = std::thread(process_file, i,  options.I1s[i].c_str(), \
                       options.R1s[i].c_str(), options.R2s[i].c_str(), \
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
    delete [] readers;
    delete [] writers;
}

void bam_writers(int windex, SAM_RECORD_BINS *samrecord_data) {
    std::thread::id  thread_id = std::this_thread::get_id();
    SamFile samOut;
    string outputfile;
    char buf[200];
    sprintf(buf, "subfile_%d.bam", windex);
    outputfile = buf;
    // open to write the outputfile
    samOut.OpenForWrite(outputfile.c_str());
    // Write the sam header.
    SamFileHeader samHeader;
    // add the HD tags
    samHeader.setHDTag("VN", "1.6");
    samHeader.setHDTag("SO", "unsorted");
    // add the RG group tags
    SamHeaderRG *headerRG = new SamHeaderRG;
    headerRG->setTag("ID", "A");
    headerRG->setTag("SM", samrecord_data->sample_id.c_str());
    samHeader.addRG(headerRG);
    samOut.WriteHeader(samHeader);

    while (true) {
      if (sem_wait(&semaphores[windex]) == -1)
         error("sem_wait:semaphores");
      SamRecord *samRecord  = \
         samrecord_data->samrecords[samrecord_data->active_thread_no];
      for (auto index : samrecord_data->file_index[samrecord_data->active_thread_no][windex]) {
           samOut.WriteRecord(samHeader, samRecord[index]);
      }

      if (sem_post(&semaphores_workers[windex]) == -1)
          error("sem_post: semaphores_workers");

      if (samrecord_data->stop) break;
    }

    // printf("closing %d\n", windex);
    samOut.Close();
}


void process_file(int tindex, String filename1, String filename2, \
                   String filename3, \
                   unsigned int barcode_length, unsigned int umi_length, \
                   const WHITE_LIST_DATA* white_list_data, \
                   SAM_RECORD_BINS *samrecord_data) {
    FastQFile fastQFile1(4, 4);
    FastQFile fastQFile2(4, 4);
    FastQFile fastQFile3(4, 4);

    if (fastQFile1.openFile(filename1, BaseAsciiMap::UNKNOWN) != \
        FastQStatus::FASTQ_SUCCESS) {
       std::cerr << "Failed to open file: " <<  filename1.c_str();
       return;
    }

    if (fastQFile2.openFile(filename2, BaseAsciiMap::UNKNOWN) != \
        FastQStatus::FASTQ_SUCCESS) {
       std::cerr << "Failed to open file: " <<  filename2.c_str();
       return;
    }

    if (fastQFile3.openFile(filename3, BaseAsciiMap::UNKNOWN) != \
        FastQStatus::FASTQ_SUCCESS) {
       std::cerr << "Failed to open file: " <<  filename3.c_str();
       return;
    }

    std::thread::id  thread_id = std::this_thread::get_id();
    SamRecord *samRecord  = samrecord_data->samrecords[tindex];
    int block_size = samrecord_data->block_size;

    // Keep reading the file until there are no more fastq sequences to process.
    int i = 0;
    int n_barcode_errors = 0;
    int n_barcode_corrected = 0;
    int n_barcode_correct = 0;
    int r = 0;
    printf("Opening the thread in %d\n", tindex);
    while (fastQFile1.keepReadingFile()) {
       if (fastQFile1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS &&
          fastQFile2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS &&
          fastQFile3.readFastQSequence() == FastQStatus::FASTQ_SUCCESS) {
          i = i + 1;
          // check the sequence names matching
          std::string a = std::string(fastQFile2.myRawSequence.c_str());
          std::string b = std::string(fastQFile2.myQualityString.c_str());
          // extract the barcode and UMI
          string barcode = a.substr(0, barcode_length);
          string UMI  = a.substr(barcode_length, umi_length);
          // extract barcode and UMI quality string
          string barcodeQString = b.substr(0, barcode_length);
          string UMIQString  = b.substr(barcode_length, umi_length);
          // reset the samrecord 
          samRecord[r].resetRecord();
          samRecord[r].setReadName(fastQFile3.mySequenceIdentifier.c_str());
          samRecord[r].setSequence(fastQFile3.myRawSequence.c_str());
          samRecord[r].setQuality(fastQFile3.myQualityString.c_str());
          // Barcode
          samRecord[r].addTag("CR", 'Z', barcode.c_str());
          samRecord[r].addTag("CY", 'Z', barcodeQString.c_str());
          // UMI
          samRecord[r].addTag("UR", 'Z', UMI.c_str());
          samRecord[r].addTag("UY", 'Z', UMIQString.c_str());

          std::string indexseq = std::string(fastQFile1.myRawSequence.c_str());
          std::string indexSeqQual = std::string(fastQFile1.myQualityString.c_str());
          samRecord[r].addTag("SR", 'Z', indexseq.c_str());
          samRecord[r].addTag("SY", 'Z', indexSeqQual.c_str());
          samRecord[r].addTag("RG", 'Z', "A");
          samRecord[r].setFlag(4);
          
          string correct_barcode;
          string bucket_barcode;
          if (white_list_data->mutations.find(barcode) != \
               white_list_data->mutations.end()) {
             if (white_list_data->mutations.at(barcode) == -1) {
                  correct_barcode = barcode;
                  n_barcode_correct += 1;
             } else {
               correct_barcode = \
               white_list_data->barcodes.at(white_list_data->mutations.at(barcode));
               n_barcode_corrected += 1;
             }
             // is used for computing the file index
             bucket_barcode = correct_barcode;
             samRecord[r].addTag("CB", 'Z', correct_barcode.c_str());
          } else {
             n_barcode_errors += 1;
             bucket_barcode = barcode;
          }
          samrecord_data->num_records[tindex]++;
          int32_t bucket = std::hash<std::string>{}(bucket_barcode.c_str()) % \
                       samrecord_data->num_files;
          samrecord_data->file_index[tindex][bucket].push_back(r);
          // write a block of samrecords
          r = r + 1;
          if (r == block_size || !fastQFile1.keepReadingFile()) {
                mtx.lock();
                samrecord_data->active_thread_no = tindex;
                //  printf("Ready to write %d\n", tindex);
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
                r = 0;
                samrecord_data->num_records[tindex] = 0;
                mtx.unlock();
          }
          if (i % 1000000 == 0) {
              printf("%d\n", i);
              std::string a = std::string(fastQFile2.myRawSequence.c_str());
              printf("%s\n", fastQFile1.mySequenceIdLine.c_str());
              printf("%s\n", fastQFile2.mySequenceIdLine.c_str());
              printf("%s\n", fastQFile3.mySequenceIdLine.c_str());
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
    // Close the input file.
    fastQFile1.closeFile();
    fastQFile2.closeFile();
    fastQFile3.closeFile();
    printf("Total barcodes:%d\n correct:%d\ncorrected:%d\nuncorrectible"
           ":%d\nuncorrected:%lf\n", \
           i, n_barcode_correct, n_barcode_corrected, n_barcode_errors, \
           n_barcode_errors/static_cast<double>(i) *100);
}

