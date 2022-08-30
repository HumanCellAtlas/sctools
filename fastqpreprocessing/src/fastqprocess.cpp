/**
 *  @file   fastqprocess.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/

#include "fastq_common.h"
#include "input_options.h"

unsigned int g_barcode_length;
unsigned int g_umi_length;

void fillSamRecord(SamRecord* samRecord, FastQFile* fastQFileI1,
                   FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                   bool has_I1_file_list)
{
  // check the sequence names matching
  std::string a = std::string(fastQFileR1->myRawSequence.c_str());
  std::string b = std::string(fastQFileR1->myQualityString.c_str());

  // extract the raw barcode and UMI
  std::string barcode_seq = a.substr(0, g_barcode_length);
  std::string umi_seq = a.substr(g_barcode_length, g_umi_length);

  // extract raw barcode and UMI quality string
  std::string barcode_quality = b.substr(0, g_barcode_length);
  std::string umi_quality = b.substr(g_barcode_length, g_umi_length);

  fillSamRecordCommon(samRecord, fastQFileI1, fastQFileR1, fastQFileR2, has_I1_file_list,
                      barcode_seq, barcode_quality, umi_seq, umi_quality);
}

std::string barcodeGetter(SamRecord* samRecord, FastQFile* fastQFileI1,
                          FastQFile* fastQFileR1, FastQFile* fastQFileR2,
                          bool has_I1_file_list)
{
  return std::string(fastQFileR1->myRawSequence.c_str()).substr(0, g_barcode_length);
}

int main(int argc, char** argv)
{
  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);
  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);

  g_barcode_length = options.barcode_length;
  g_umi_length = options.umi_length;

  mainCommon(options.white_list_file, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.sample_id,
             fillSamRecord, barcodeGetter);
  return 0;
}
