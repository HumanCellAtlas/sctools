#include "fastq_common.h"
#include "input_options.h"

std::vector<std::pair<char, int>> parseReadStructure(std::string const& read_structure)
{
  std::vector<std::pair<char, int>> ret;
  int next_ind = 0;
  while (next_ind < read_structure.size())
  {
    int type_ind = read_structure.find_first_not_of("0123456789", next_ind);
    assert(type_ind != std::string::npos);
    char type = read_structure[type_ind];
    int len = std::stoi(read_structure.substr(next_ind, type_ind - next_ind));
    ret.emplace_back(type, len);
    next_ind = type_ind + 1;
  }
  return ret;
}

std::vector<std::pair<char, int>> g_parsed_read_structure;

void fillSamRecordWithReadStructure(SamRecord* sam, FastQFile* fastQFileI1,
                                    FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                                    bool has_I1_file_list)
{
  // check the sequence names matching
  std::string a = std::string(fastQFileR1->myRawSequence.c_str());
  std::string b = std::string(fastQFileR1->myQualityString.c_str());
  // extract the raw barcode and UMI 8C18X6C9M1X and raw barcode and UMI quality string

  std::string barcode_seq, barcode_quality, umi_seq, umi_quality;
  int cur_ind = 0;
  for (auto [tag, length] : g_parsed_read_structure)
  {
    switch (tag)
    {
    case 'C':
      barcode_seq += a.substr(cur_ind, length);
      barcode_quality += b.substr(cur_ind, length);
      break;
    case 'M':
      umi_seq += a.substr(cur_ind, length);
      umi_quality += b.substr(cur_ind, length);
      break;
    default:
      break;
    }
    cur_ind += length;
  }
  fillSamRecordCommon(sam, fastQFileI1, fastQFileR1, fastQFileR2, has_I1_file_list,
                      barcode_seq, barcode_quality, umi_seq, umi_quality);
}

std::string slideseqBarcodeGetter(SamRecord* sam, FastQFile* fastQFileI1,
                                  FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                                  bool has_I1_file_list)
{
  return std::string(sam->getString("CR").c_str());
}

void outputHandler(WriteQueue* cur_write_queue, SamRecord* samrec, int reader_thread_index)
{
  cur_write_queue->enqueueWrite(std::make_pair(samrec, reader_thread_index));
}

int main(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options = readOptionsFastqSlideseq(argc, argv);
  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);

  g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.sample_id,
             fillSamRecordWithReadStructure, slideseqBarcodeGetter, outputHandler);
  return 0;
}
