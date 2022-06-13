#ifndef __FASTQ_METRICS_H__
#define __FASTQ_METRICS_H__
/**
 *  @file   fastq_metrics.h
 *  @brief  functions for computing metrics
 *  @author Farzaneh Khajouei and Fred Douglas
 *  @date   2022-05-25
 ***********************************************/

#include <string>
#include <unordered_map>
#include <vector>
#include <thread>

#include "utilities.h"
#include "input_options.h"

class PositionWeightMatrix
{
public:
    PositionWeightMatrix(int length): A(length), C(length), G(length), T(length), N(length){}
    void recordChunk(std::string s);
    PositionWeightMatrix& operator+=(const PositionWeightMatrix& rhs);
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
    FastQMetricsShard(std::string read_structure) : read_structure_(read_structure),
                     barcode_length_(getLengthOfType(read_structure_,'C')),
                     umi_length_(getLengthOfType(read_structure_,'M')),
                     tagged_lengths_(parseReadStructure(read_structure_)),
                     barcode_(barcode_length),
                     umi_(umi_length){}
    void ingestBarcodeAndUMI(std::string raw_seq);
    void processShard( std::string filenameR1, std::string read_structure, const WHITE_LIST_DATA* white_list_data);
    static void mergeMetricsShardsToFile(std::string filename_prefix, vector<FastQMetricsShard> shards, int umi_length, int CB_length);
    FastQMetricsShard& operator+=(const FastQMetricsShard& rhs);


private:
    std::string read_structure_;
    int barcode_length_;
    int umi_length_;
    std::vector<std::pair<char, int>> tagged_lengths_;
    std::unordered_map<string,int> barcode_counts_;
    std::unordered_map<string,int> umi_counts_;
    PositionWeightMatrix barcode_;
    PositionWeightMatrix umi_;
};

#endif // __FASTQ_METRICS_H__