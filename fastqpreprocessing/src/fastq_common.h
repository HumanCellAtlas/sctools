#ifndef __SCTOOLS_FASTQPREPROCESSING_FASTQ_COMMON_H_
#define __SCTOOLS_FASTQPREPROCESSING_FASTQ_COMMON_H_

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <string>
#include <vector>

#include "FastQFile.h"
#include "FastQStatus.h"
#include "SamFile.h"
#include "SamValidation.h"

// A pointer to a valid SamRecord waiting to be written to disk, and the index
// of the g_read_arenas that pointer should be released to after the write.
using PendingWrite = std::pair<SamRecord*, int>;

class WriteQueue
{
public:
  static constexpr int kShutdown = -1;
  PendingWrite dequeueWrite();
  void enqueueWrite(PendingWrite write);
  void enqueueShutdownSignal();
private:
  std::mutex mutex_;
  std::condition_variable cv_;
  std::queue<PendingWrite> queue_;
};

void fillSamRecordCommon(SamRecord* samRecord, FastQFile* fastQFileI1,
                         FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                         bool has_I1_file_list,
                         std::string const& barcode_seq, std::string const& barcode_quality,
                         std::string const& umi_seq, std::string const& umi_quality);

void mainCommon(
    std::string white_list_file, int num_writer_threads, std::string output_format,
    std::vector<std::string> I1s, std::vector<std::string> R1s, std::vector<std::string> R2s,
    std::string sample_id,
    std::function<void(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> sam_record_filler,
    std::function<std::string(SamRecord*, FastQFile*, FastQFile*, FastQFile*, bool)> barcode_getter,
    std::function<void(WriteQueue*, SamRecord*, int)> output_handler);

#endif // __SCTOOLS_FASTQPREPROCESSING_FASTQ_COMMON_H_
