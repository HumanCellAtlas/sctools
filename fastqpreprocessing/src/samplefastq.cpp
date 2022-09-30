#include "fastq_common.h"
#include "input_options.h"
#include <fstream>

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

void noOpFillSamRecordWithReadStructure(SamRecord* sam, FastQFile* fastQFileI1,
                                    FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                                    bool has_I1_file_list)
{

}

std::string slideseqBarcodeGetter(SamRecord* sam, FastQFile* fastQFileI1,
                                  FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                                  bool has_I1_file_list)
{
  // check the sequence names matching
  std::string a = std::string(fastQFileR1->myRawSequence.c_str());
  // extract the raw barcode 8C18X6C9M1X

  std::string barcode_seq;
  int cur_ind = 0;
  for (auto [tag, length] : g_parsed_read_structure)
  {
    switch (tag)
    {
    case 'C':
      barcode_seq += a.substr(cur_ind, length);
      break;
    default:
      break;
    }
    cur_ind += length;
  }
  return barcode_seq;
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

  std::ofstream outfile_r1("sampled_down.R1");
  if (!outfile_r1)
    crash("Failed to open output file sampled_down.R1");
  std::ofstream outfile_r2("sampled_down.R2");
  if (!outfile_r2)
    crash("Failed to open output file sampled_down.R2");

  g_parsed_read_structure = parseReadStructure(options.read_structure);
  mainCommon(options.white_list_file, /*num_writer_threads=*/1, options.output_format,
             options.I1s, options.R1s, options.R2s, options.sample_id,
             noOpFillSamRecordWithReadStructure, slideseqBarcodeGetter,
             [&outfile_r1, &outfile_r2](WriteQueue* ignored1, SamRecord* sam, int reader_thread_index)
             {
               // Assumed read structure of 8C18X6C9M1X with a fixed spacer sequence
               const char* barcode = sam->getString("CR").c_str();
               const char* quality_score = sam->getString("CY").c_str();
               outfile_r1 << "@" << sam->getReadName() << "\n"
                      << std::string_view(barcode, 8) << "CTTCAGCGTTCCCGAGAG" << std::string_view(barcode+8, 6) << sam->getString("UR") <<"T\n"
                      << "+\n"
                      << std::string_view(quality_score, 8)<<"FFFFFFFFFFFFFFFFFF" << std::string_view(quality_score+8, 6) << sam->getString("UY") <<"F"<< "\n";

               outfile_r2 << "@" << sam->getReadName() << "\n"
                      << sam->getSequence() << "\n"
                      << "+\n"
                      << sam->getQuality() << "\n";
	       releaseReaderThreadMemory(reader_thread_index,sam);
             });
  return 0;
}
