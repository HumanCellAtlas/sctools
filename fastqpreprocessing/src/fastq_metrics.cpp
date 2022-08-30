/**
 *  @file   fastq_metrics.cpp
 *  @brief  functions for computing metrics
 *  @author Farzaneh Khajouei and Fred Douglas
 *  @date   2022-05-25
 ***********************************************/
#include "FastQFile.h"
#include "FastQStatus.h"
#include "fastq_metrics.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cassert>

using std::string;

std::vector<std::pair<char, int>> parseReadStructure(std::string read_structure)
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

int getLengthOfType(string read_structure,char type)
{
  int total_length = 0;
  for (auto [curr_type, length] : parseReadStructure(read_structure))
    if (curr_type == type)
      total_length += length;
  return total_length;
}

void PositionWeightMatrix::recordChunk(string s)
{
  for (int index = 0; index < s.size(); index++)
  {
    switch (s[index])
    {
    case 'A':
    case 'a':
      A[index]++;
      break;
    case 'C':
    case 'c':
      C[index]++;
      break;
    case 'G':
    case 'g':
      G[index]++;
      break;
    case 'T':
    case 't':
      T[index]++;
      break;
    case 'N':
    case 'n':
      N[index]++;
      break;
    default:
      std::cerr<<"Unknown character:"<<s[index]<<std::endl;
    }
  }
}

FastQMetricsShard::FastQMetricsShard(std::string read_structure)
  : read_structure_(read_structure),
    barcode_length_(getLengthOfType(read_structure_,'C')),
    umi_length_(getLengthOfType(read_structure_,'M')),
    tagged_lengths_(parseReadStructure(read_structure_)),
    barcode_(barcode_length_),
    umi_(umi_length_) {}

// Read a chunk from a fastq r1 and get UMI and Cellbarcode filled
void FastQMetricsShard::ingestBarcodeAndUMI(std::string_view raw_seq)
{
  // extract the raw barcode and UMI 8C18X6C9M1X and raw barcode and UMI quality string
  std::string barcode_seq, umi_seq;
  int cur_ind = 0;
  for (auto [tag, length] : tagged_lengths_)
  {
    switch (tag)
    {
    case 'C':
      barcode_seq += raw_seq.substr(cur_ind, length);
      break;
    case 'M':
      umi_seq += raw_seq.substr(cur_ind, length);
      break;
    default:
      break;
    }
    cur_ind += length;
  }

  barcode_counts_[barcode_seq]++;
  umi_counts_[umi_seq]++;
  barcode_.recordChunk(barcode_seq);
  umi_.recordChunk(umi_seq);
}


// This is a wrapper to use std thread
void processShard(FastQMetricsShard* fastq_metrics_shard, String filenameR1,
                  std::string read_structure, const WhiteListData* white_list_data)
{
  fastq_metrics_shard->processShard(filenameR1, read_structure, white_list_data);
}
void FastQMetricsShard::processShard(String filenameR1, std::string read_structure,
                                     const WhiteListData* white_list_data)
{
  /// setting the shortest sequence allowed to be read
  FastQFile fastQFileR1(4, 4);
  // open the R1 file
  if (fastQFileR1.openFile(filenameR1, BaseAsciiMap::UNKNOWN) != FastQStatus::FASTQ_SUCCESS)
    crash("Failed to open R1 file");

  // Keep reading the file until there are no more fastq sequences to process.
  int n_lines_read = 0;
  while (fastQFileR1.keepReadingFile())
  {
    if (fastQFileR1.readFastQSequence() != FastQStatus::FASTQ_SUCCESS)
      break;

    ingestBarcodeAndUMI(std::string_view(fastQFileR1.myRawSequence.c_str(),fastQFileR1.myRawSequence.Length()));

    n_lines_read++;
    if (n_lines_read % 10000000 == 0)
    {
      printf("%d\n", n_lines_read);
      std::string a = std::string(fastQFileR1.myRawSequence.c_str());
      printf("%s\n", fastQFileR1.mySequenceIdLine.c_str());
    }
  }
  // Finished processing all of the sequences in the file.
  // Close the input files.
  fastQFileR1.closeFile();
}

PositionWeightMatrix& PositionWeightMatrix::operator+=(const PositionWeightMatrix& rhs)
{
  for (int i=0; i < A.size(); i++)
  {
    A[i] += rhs.A[i];
    C[i] += rhs.C[i];
    G[i] += rhs.G[i];
    T[i] += rhs.T[i];
    N[i] += rhs.N[i];
  }
  return *this;
}

FastQMetricsShard& FastQMetricsShard::operator+=(const FastQMetricsShard& rhs)
{
  for (auto [key,value] : rhs.barcode_counts_)
    barcode_counts_[key] += value;
  for (auto [key,value] : rhs.umi_counts_)
    umi_counts_[key] += value;

  barcode_+=rhs.barcode_;
  umi_+=rhs.umi_;
  return *this;
}

/** @copydoc process_inputs */
void process_inputs(const INPUT_OPTIONS_FASTQ_READ_STRUCTURE& options,
                    const WhiteListData* white_list_data)
{
  // number of files based on the input size
  int num_files = options.R1s.size();

  // compute UMI and cell_barcode lengths

  int umi_length = getLengthOfType(options.read_structure,'M');
  int CB_length = getLengthOfType(options.read_structure,'C');

  // create the data for the threads
  vector <FastQMetricsShard> fastqMetrics;
  for (int i = 0; i < num_files; i++)
    fastqMetrics.emplace_back(options.read_structure);

  // execute the fastq readers threads
  vector<std::thread> readers;
  for (unsigned int i = 0; i < options.R1s.size(); i++)
  {
    readers.emplace_back(processShard,
                         &fastqMetrics[i],
                         options.R1s[i].c_str(),
                         options.read_structure.c_str(),
                         white_list_data);

  }

  // every reader thread joins.
  for (unsigned int i = 0; i < options.R1s.size(); i++)
    readers[i].join();

  std::cout << "Done reading all shards. Will now aggregate and write to file; "
            << "this will take a few minutes." << std::endl;
  FastQMetricsShard::mergeMetricsShardsToFile(options.sample_id, fastqMetrics, umi_length, CB_length);
}

void writeCountsFile(std::unordered_map<string,int> counts, std::string filename)
{
  std::ofstream out(filename, std::ofstream::out);
  std::vector<std::pair<std::string,int>> sorted_counts;
  for (auto [str, count] : counts)
    sorted_counts.emplace_back(str, count);
  std::sort(sorted_counts.begin(), sorted_counts.end(), //sort counts from most to fewest!
            [](std::pair<std::string,int> const& a, std::pair<std::string,int> const& b)
  {
    return a.second > b.second;
  });
  for (auto [str, count] : sorted_counts)
    out << count << "\t" << str << "\n";
}
void PositionWeightMatrix::writeToFile(std::string filename)
{
  std::ofstream out(filename, std::ofstream::out);
  out << "position\tA\tC\tG\tT\tN\n";
  for (int i = 0; i < A.size(); i++)
    out << (i + 1) << "\t" << A[i] << "\t" << C[i] << "\t" << G[i] << "\t" << T[i] << "\t" << N[i] << "\n";
}
void FastQMetricsShard::mergeMetricsShardsToFile(std::string filename_prefix, vector<FastQMetricsShard> shards, int umi_length, int CB_length)
{
  FastQMetricsShard total(shards[0].read_structure_);
  for (FastQMetricsShard const& shard : shards)
    total += shard;

  writeCountsFile(total.umi_counts_, filename_prefix + ".numReads_perCell_XM.txt");
  writeCountsFile(total.barcode_counts_, filename_prefix + ".numReads_perCell_XC.txt");
  total.barcode_.writeToFile(filename_prefix + ".barcode_distribution_XC.txt");
  total.umi_.writeToFile(filename_prefix + ".barcode_distribution_XM.txt");
}

int main(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options = readOptionsFastqMetrics(argc, argv);
  std::cout << "reading whitelist file " << options.white_list_file << "...";
  WhiteListData white_list_data = readWhiteList(options.white_list_file);
  std::cout << "done" << std::endl;

  process_inputs(options, &white_list_data);
  return 0;
}
