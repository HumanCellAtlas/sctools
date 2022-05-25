#ifndef __FASTQ_METRICS_H__
#define __FASTQ_METRICS_H__
/**
 *  @file   fastq_metrics.h
 *  @brief  functions for computing metrics
 *  @author Farzaneh Khajouei
 *  @date   2022-05-25
 ***********************************************/
#include "SamFile.h"

#include <string>
#include <unordered_map>
#include <vector>

#include "utilities.h"
#include "input_options.h"

struct PositionWeightMatrix
{
public:
    PositionWeightMatrix(int length): A(length), C(length), G(length), T(length), N(length){}
    std::vector<int> A;
    std::vector<int> C;
    std::vector<int> G;
    std::vector<int> T;
    std::vector<int> N;
};

class FastQMetricsShard
{
public:
    FastQMetricsShard(int barcode_length, int umi_length) : barcode_(barcode_length), umi_(umi_length){}
    void ingestSamRecord(const SamRecord* sam_record);
private:
    std::unordered_map<string,int> barcode_counts_;
    std::unordered_map<string,int> umi_counts_;
    PositionWeightMatrix barcode_;
    PositionWeightMatrix umi_;
};

void mergeMetricsShardsToFile(std::string filename, vector<FastQMetricsShard> metrics_shards );
#endif // __FASTQ_METRICS_H__