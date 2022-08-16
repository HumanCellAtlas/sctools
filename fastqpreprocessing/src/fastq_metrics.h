#ifndef __FASTQ_METRICS_H__
#define __FASTQ_METRICS_H__
/**
 *  @file   fastq_metrics.h
 *  @brief  computes metrics based on the read structure
 *  @author Farzaneh Khajouei and Fred Douglas
 *  @date   2022-05-25
 ***********************************************/
#include <string>
#include <unordered_map>
#include <vector>
#include <thread>
#include "BaseAsciiMap.h"
#include "utilities.h"
#include "input_options.h"
#include "FastQFile.h"
#include "FastQStatus.h"

class PositionWeightMatrix
{
public:
  PositionWeightMatrix(int length): A(length), C(length), G(length), T(length), N(length) {}
  void recordChunk(std::string s);
  PositionWeightMatrix& operator+=(PositionWeightMatrix const& rhs);
  void writeToFile(std::string filename);

  std::vector<int> A;
  std::vector<int> C;
  std::vector<int> G;
  std::vector<int> T;
  std::vector<int> N;
};

class FastQMetricsShard
{
public:
  FastQMetricsShard(std::string read_structure);
  void ingestBarcodeAndUMI(std::string_view raw_seq);
  void processShard(String filenameR1, std::string read_structure,
                    const WhiteListData* white_list_data);
  static void mergeMetricsShardsToFile(std::string filename_prefix,
                                       std::vector<FastQMetricsShard> shards,
                                       int umi_length, int CB_length);
  FastQMetricsShard& operator+=(FastQMetricsShard const& rhs);


private:
  std::string read_structure_;
  int barcode_length_;
  int umi_length_;
  std::vector<std::pair<char, int>> tagged_lengths_;
  std::unordered_map<std::string, int> barcode_counts_;
  std::unordered_map<std::string, int> umi_counts_;
  PositionWeightMatrix barcode_;
  PositionWeightMatrix umi_;
};

#endif // __FASTQ_METRICS_H__
