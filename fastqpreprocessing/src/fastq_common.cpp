#include <gzstream.h>
#include <iostream>
#include <fstream>
#include <cstdint>

// number of samrecords per buffer in each reader
constexpr size_t kSamRecordBufferSize = 10000;

#include "input_options.h"
#include "utilities.h"

#include "FastQFile.h"
#include "FastQStatus.h"
#include "BaseAsciiMap.h"
#include "SamFile.h"
#include "SamValidation.h"

#include <thread>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <condition_variable>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <functional>
#include <mutex>
#include <stack>

// Overview of multithreading:
// * There are reader threads and writer threads. (Writers are either fastq or
//   bam, depending on how the program was run).
// * Each {reader, writer} has its own {input, output} file.
// * Each reader has an entry in g_read_arenas, and each writer has an entry in
//   g_write_queues.
// * Readers load each chunk of their processed results into SamRecord pointers
//   loaned out by their arena. They put the pointer in the correct write queue.
// * When a write queue finishes writing a SamRecord to the file, it notifies
//   the record pointer's arena that the record's memory is no longer in use.
//   The arena can then give that pointer to its reader for a new read.

// A pointer to a valid SamRecord waiting to be written to disk, and the index
// of the g_read_arenas that pointer should be released to after the write.
using PendingWrite = std::pair<SamRecord*, int>;

class WriteQueue
{
public:
  PendingWrite dequeueWrite()
  {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&] { return !queue_.empty(); });
    auto pair = queue_.front();
    queue_.pop();
    return pair;
  }
  void enqueueWrite(PendingWrite write)
  {
    mutex_.lock();
    queue_.push(write);
    mutex_.unlock();
    cv_.notify_one();
  }
  void shutdown()
  {
    mutex_.lock();
    queue_.push(std::make_pair(nullptr, -1));
    mutex_.unlock();
    cv_.notify_one();
  }
private:
  std::mutex mutex_;
  std::condition_variable cv_;
  std::queue<PendingWrite> queue_;
};

std::vector<std::unique_ptr<WriteQueue>> g_write_queues;

// I wrote this class to stay close to the performance characteristics of the
// original code, but I suspect the large buffers might not be necessary.
// If it doesn't slow things down noticeably, it would be cleaner to just delete
// this class, and have the WriteQueue accept unique_ptr<SamRecord> (with the
// addition of some reasonable bound on how much WriteQueue can have
// outstanding; maybe kSamRecordBufferSize items), and let them be directly
// destroyed after writing rather than be reused with this arena approach.
class SamRecordArena
{
public:
  SamRecordArena()
  {
    for (int i = 0; i < kSamRecordBufferSize; i++)
      samrecords_memory_.push_back(std::make_unique<SamRecord>());

    for (int i = samrecords_memory_.size() - 1; i >= 0; i--)
      available_samrecords_.push(samrecords_memory_[i].get());
  }

  SamRecord* acquireSamRecordMemory()
  {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&] { return !available_samrecords_.empty(); });
    SamRecord* sam = available_samrecords_.top();
    available_samrecords_.pop();
    return sam;
  }
  void releaseSamRecordMemory(SamRecord* sam)
  {
    mutex_.lock();
    available_samrecords_.push(sam);
    mutex_.unlock();
    cv_.notify_one();
  }
private:
  std::vector<std::unique_ptr<SamRecord>> samrecords_memory_;
  std::mutex mutex_;
  std::condition_variable cv_;
  // Reusing most-recently-used memory first ought to be more cache friendly.
  std::stack<SamRecord*> available_samrecords_;
};

std::vector<std::unique_ptr<SamRecordArena>> g_read_arenas;



void writeFastqRecord(ogzstream& r1_out, ogzstream& r2_out, SamRecord* sam)
{
  r1_out << "@" << sam->getReadName() << "\n" << sam->getString("CR").c_str()
         << sam->getString("UR") << "\n+\n" << sam->getString("CY") << sam->getString("UY") << "\n";
  r2_out << "@" << sam->getReadName() << "\n" << sam->getSequence() << "\n+\n"
         << sam->getQuality() << "\n";
}

void fastqWriterThread(int write_thread_index)
{
  std::string r1_output_fname = "fastq_R1_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r1_out(r1_output_fname.c_str());
  if (!r1_out)
    crash("ERROR: Failed to open R1 fastq file " + r1_output_fname + " for writing");

  std::string r2_output_fname = "fastq_R2_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r2_out(r2_output_fname.c_str());
  if (!r2_out)
    crash("ERROR: Failed to open R2 fastq file " + r2_output_fname + " for writing");

  while (true)
  {
    auto [sam, source_reader_index] = g_write_queues[write_thread_index]->dequeueWrite();
    if (source_reader_index == -1)
      break;

    writeFastqRecord(r1_out, r2_out, sam);
    g_read_arenas[source_reader_index]->releaseSamRecordMemory(sam);
  }

  // close the fastq files
  r1_out.close();
  r2_out.close();
}

void bamWriterThread(int write_thread_index, std::string sample_id)
{
  std::string bam_out_fname = "subfile_" + std::to_string(write_thread_index) + ".bam";
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
  headerRG->setTag("SM", sample_id.c_str());
  samHeader.addRG(headerRG);

  // add the header to the output bam
  samOut.WriteHeader(samHeader);

  while (true)
  {
    auto [sam, source_reader_index] = g_write_queues[write_thread_index]->dequeueWrite();
    if (source_reader_index == -1)
      break;

    samOut.WriteRecord(samHeader, *sam);
    g_read_arenas[source_reader_index]->releaseSamRecordMemory(sam);
  }

  // close the bamfile
  samOut.Close();
}

void fillSamRecordCommon(SamRecord* samRecord, FastQFile* fastQFileI1,
                         FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                         bool has_I1_file_list,
                         std::string const& barcode_seq, std::string const& barcode_quality,
                         std::string const& umi_seq, std::string const& umi_quality)
{
  // reset the samrecord
  samRecord->resetRecord();
  // add read group and the sam flag
  samRecord->addTag("RG", 'Z', "A");
  samRecord->setFlag(4);
  // add identifier, sequence and quality score of the alignments
  samRecord->setReadName(fastQFileR2->mySequenceIdentifier.c_str());
  samRecord->setSequence(fastQFileR2->myRawSequence.c_str());
  samRecord->setQuality(fastQFileR2->myQualityString.c_str());
  // add barcode and quality
  samRecord->addTag("CR", 'Z', barcode_seq.c_str());
  samRecord->addTag("CY", 'Z', barcode_quality.c_str());
  // add UMI
  samRecord->addTag("UR", 'Z', umi_seq.c_str());
  samRecord->addTag("UY", 'Z', umi_quality.c_str());
  // add raw sequence and quality sequence for the index
  if (has_I1_file_list)
  {
    samRecord->addTag("SR", 'Z', fastQFileI1->myRawSequence.c_str());
    samRecord->addTag("SY", 'Z', fastQFileI1->myQualityString.c_str());
  }
}

// Computes the whitelist-corrected barcode and adds it to sam_record.
// Returns the index of the bamfile bucket / writer thread where sam_record
// should be sent.
int32_t correctBarcodeToWhitelist(
    const std::string& barcode, SamRecord* sam_record, const WhiteListData* white_list_data,
    int* n_barcode_corrected, int* n_barcode_correct, int* n_barcode_errors, int num_writer_threads)
{
  std::string correct_barcode;
  // bucket barcode is used to pick the target bam file
  // This is done because in the case of incorrigible barcodes
  // we need a mechanism to uniformly distribute the alignments
  // so that no bam is oversized to putting all such barcode less
  // sequences into one particular. Incorregible barcodes are simply
  // added withouth the CB tag
  std::string bucket_barcode;
  if (auto it = white_list_data->mutations.find(barcode) ; it != white_list_data->mutations.end())
  {
    int64_t mutation_index = it->second;
    if (mutation_index == -1) // -1 means raw barcode is correct
    {
      correct_barcode = barcode;
      *n_barcode_correct += 1;
    }
    else
    {
      // it is a 1-mutation of some whitelist barcode so get the
      // barcode by indexing into the vector of whitelist barcodes
      correct_barcode = white_list_data->barcodes[mutation_index];
      *n_barcode_corrected += 1;
    }
    // is used for computing the file index
    bucket_barcode = correct_barcode;

    // corrected barcode should be added to the samrecord
    sam_record->addTag("CB", 'Z', correct_barcode.c_str());
  }
  else     // not possible to correct the raw barcode
  {
    *n_barcode_errors += 1;
    bucket_barcode = barcode;
  }
  // destination bam file index computed based on the bucket_barcode
  return std::hash<std::string> {}(bucket_barcode) % num_writer_threads;
}

// Returns true if successfully read a sequence.
bool readOneItem(FastQFile& fastQFileI1, bool has_I1_file_list,
                   FastQFile& fastQFileR1, FastQFile& fastQFileR2)
{
  return (!has_I1_file_list ||
      (
        has_I1_file_list &&
        fastQFileI1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
      )
  )
  && fastQFileR1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
  && fastQFileR2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS;
}

void fastQFileReaderThread(
    int reader_thread_index, std::string filenameI1, String filenameR1,
    String filenameR2, const WhiteListData* white_list_data,
    std::function <void(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> sam_record_filler,
    std::function <std::string(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> barcode_getter)
{
  /// setting the shortest sequence allowed to be read
  FastQFile fastQFileI1(4, 4);
  FastQFile fastQFileR1(4, 4);
  FastQFile fastQFileR2(4, 4);

  bool has_I1_file_list = true;
  if (!filenameI1.empty())
  {
    if (fastQFileI1.openFile(String(filenameI1.c_str()), BaseAsciiMap::UNKNOWN) !=
        FastQStatus::FASTQ_SUCCESS)
    {
      crash(std::string("Failed to open file: ") + filenameI1);
    }
  }
  else
    has_I1_file_list = false;

  if (fastQFileR1.openFile(filenameR1, BaseAsciiMap::UNKNOWN) !=
      FastQStatus::FASTQ_SUCCESS)
  {
    crash(std::string("Failed to open file: ") + filenameR1.c_str());
  }
  if (fastQFileR2.openFile(filenameR2, BaseAsciiMap::UNKNOWN) !=
      FastQStatus::FASTQ_SUCCESS)
  {
    crash(std::string("Failed to open file: ") + filenameR2.c_str());
  }

  // Keep reading the file until there are no more fastq sequences to process.
  int total_reads = 0;
  int n_barcode_errors = 0;
  int n_barcode_corrected = 0;
  int n_barcode_correct = 0;
  printf("Opening the thread in %d\n", reader_thread_index);

  while (fastQFileR1.keepReadingFile())
  {
    if (readOneItem(fastQFileI1, has_I1_file_list, fastQFileR1, fastQFileR2))
    {
      total_reads++;

      SamRecord* samrec = g_read_arenas[reader_thread_index]->acquireSamRecordMemory();

      // prepare the samrecord with the sequence, barcode, UMI, and their quality sequences
      sam_record_filler(samrec, &fastQFileI1, &fastQFileR1, &fastQFileR2, has_I1_file_list);
      std::string barcode = barcode_getter(samrec, &fastQFileI1, &fastQFileR1, &fastQFileR2, has_I1_file_list);

      // bucket barcode is used to pick the target bam file
      // This is done because in the case of incorrigible barcodes
      // we need a mechanism to uniformly distribute the alignments
      // so that no bam is oversized to putting all such barcode less
      // sequences into one particular. Incorregible barcodes are simply
      // added withouth the CB tag
      int32_t bam_bucket = correctBarcodeToWhitelist(
          barcode, samrec, white_list_data, &n_barcode_corrected, &n_barcode_correct,
          &n_barcode_errors, g_write_queues.size());

      g_write_queues[bam_bucket]->enqueueWrite(std::make_pair(samrec, reader_thread_index));

      if (total_reads % 10000000 == 0)
      {
        printf("%d\n", total_reads);
        std::string a = std::string(fastQFileR1.myRawSequence.c_str());
        printf("%s\n", fastQFileR1.mySequenceIdLine.c_str());
        printf("%s\n", fastQFileR2.mySequenceIdLine.c_str());
      }
    }
  }

  // Finished processing all of the sequences in the file.
  // Close the input files.
  if (has_I1_file_list)
    fastQFileI1.closeFile();
  fastQFileR1.closeFile();
  fastQFileR2.closeFile();
  printf("Total barcodes:%d\n correct:%d\ncorrected:%d\nuncorrectible"
         ":%d\nuncorrected:%lf\n",
         total_reads, n_barcode_correct, n_barcode_corrected, n_barcode_errors,
         n_barcode_errors/static_cast<double>(total_reads) * 100);
}

void mainCommon(
    std::string white_list_file, int num_writer_threads, std::string output_format,
    std::vector<std::string> I1s, std::vector<std::string> R1s, std::vector<std::string> R2s,
    std::string sample_id,
    std::function <void(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> sam_record_filler,
    std::function <std::string(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> barcode_getter)
{
  std::cout << "reading whitelist file " << white_list_file << "...";
  // stores barcode correction map and vector of correct barcodes
  WhiteListData white_list_data = readWhiteList(white_list_file);
  std::cout << "done" << std::endl;


  for (int i = 0; i < num_writer_threads; i++)
    g_write_queues.push_back(std::make_unique<WriteQueue>());

  // execute the bam file writers threads
  std::vector<std::thread> writers;
  if (output_format == "BAM")
    for (int i = 0; i < num_writer_threads; i++)
      writers.emplace_back(bamWriterThread, i, sample_id);
  else if (output_format == "FASTQ")
    for (int i = 0; i < num_writer_threads; i++)
      writers.emplace_back(fastqWriterThread, i);
  else
    crash("ERROR: Output-format must be either FASTQ or BAM");

  // execute the fastq readers threads
  std::vector<std::thread> readers;

  for (unsigned int i = 0; i < R1s.size(); i++)
  {
    assert(I1s.empty() || I1s.size() == R1s.size());
    // if there is no I1 file then send an empty file name
    std::string I1 = I1s.empty() ? "" : I1s[i];

    g_read_arenas.push_back(std::make_unique<SamRecordArena>());
    readers.emplace_back(fastQFileReaderThread, i, I1.c_str(), R1s[i].c_str(),
                         R2s[i].c_str(), &white_list_data, sam_record_filler, barcode_getter);
  }

  for (auto& reader : readers)
    reader.join();

  // Now that there's nothing left to read, we can safely append a shutdown
  // signal to all the write queues.
  for (auto& write_queue : g_write_queues)
    write_queue.shutdown();

  for (auto& writer : writers)
    writer.join();
}
