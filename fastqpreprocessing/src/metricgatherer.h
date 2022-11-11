#ifndef __METRIC_GATHERER__
#define __METRIC_GATHERER__

#include <unordered_map>
#include <string>
#include <regex>
#include <iostream>
#include <vector>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <unordered_set>

#include "input_options.h"

class OnlineGaussianSufficientStatistic
{
public:
  // incorporates new_value into the online estimate of mean and variance
  void update(double new_value)
  {
    _count += 1.0;
    _sum += new_value;
    sum_EX2 += (new_value*new_value);
  }

  // return the mean value
  double getMean()
  {
    _mean = _sum/_count;
    return _mean;
  }

  // calculate and return the variance
  double calculateVariance()
  {
    if (_count < 2)
      return -1.0;
    return sum_EX2 / (_count - 1) - (_sum/_count) * (_sum / (_count - 1));
  }

  void clear()
  {
    _mean_squared_error = 0.0;
    _mean = 0.0;
    _count = 0;
    _sum = 0;
    sum_EX2 = 0.0;
  }

private:
  double _mean_squared_error = 0.0;
  double sum_EX2 = 0.0;
  double _mean = 0.0;
  double _sum = 0.0;
  double _count = 0.0;
};

// TODO better name
class LineFields
{
public:
  // TODO merge with splitStringToFields?
  LineFields(std::string const& s);
  std::string_view first_tag; // 0
  std::string_view second_tag; // 1
  std::string_view third_tag;  // 2
  std::string_view reference; // 3
  std::string_view alignment_location; // 4
  std::string_view position_str; // 5
  int is_strand; // 6
  // 7 unused
  float cell_barcode_base_above_30; // 8
  float genomic_read_quality; // 9
  float genomic_reads_base_quality_above_30; // 10
  int number_mappings; // 11
  int perfect_molecule_barcode; // 12
  // cigar N field (3) indicates a read is spliced if the value is non-zero
  int read_spliced; // 13
  int read_is_duplicate; // 14
  int cell_barcode_perfect; // 15
  float molecule_barcode_base_above_30; // 16

private:
  int getNextAsInt();
  float getNextAsFloat();
  std::string_view getNextField();

  std::string const& s_;
  size_t cur_start_ = 0;
  size_t cur_tab_ = std::string::npos;
  int fields_gotten_ = 0;
};

class MetricGatherer
{
public:
  MetricGatherer(MetricType metric_type);
  virtual ~MetricGatherer();
  //  get the headers
  virtual std::string getHeader() = 0;

  void ingestLine(std::string const& str, std::ofstream& fmetric_out);

  void output_metrics(std::ofstream& fmetric_out);
  virtual void output_metrics_extra(std::ofstream& fmetric_out) = 0;
  virtual void parse_extra_fields(LineFields const& fields) = 0;
  virtual void clear();

protected:
  const std::string kCommonHeaders[24] =
  {
    "n_reads",
    "noise_reads",
    "perfect_molecule_barcodes",
    "reads_mapped_exonic",
    "reads_mapped_intronic",
    "reads_mapped_utr",
    "reads_mapped_uniquely",
    "reads_mapped_multiple",
    "duplicate_reads",
    "spliced_reads",
    "antisense_reads",
    "molecule_barcode_fraction_bases_above_30_mean",
    "molecule_barcode_fraction_bases_above_30_variance",
    "genomic_reads_fraction_bases_quality_above_30_mean",
    "genomic_reads_fraction_bases_quality_above_30_variance",
    "genomic_read_quality_mean",
    "genomic_read_quality_variance",
    "n_molecules",
    "n_fragments",
    "reads_per_molecule",
    "reads_per_fragment",
    "fragments_per_molecule",
    "fragments_with_single_read_evidence",
    "molecules_with_single_read_evidence"
  };

  void parseAlignedReadFields(LineFields const& fields, std::string hyphenated_tags);

private:
  const MetricType metric_type_;

  // count information
  int n_reads_ = 0;
  const int noise_reads = 0; //# long polymers, N-sequences; NotImplemented

  std::unordered_map<std::string, int> fragment_histogram_;
  std::unordered_map<std::string, int> molecule_histogram_;

  // molecule information
  OnlineGaussianSufficientStatistic molecule_barcode_fraction_bases_above_30_;

  int perfect_molecule_barcodes_ = 0;

  OnlineGaussianSufficientStatistic genomic_reads_fraction_bases_quality_above_30_;

  OnlineGaussianSufficientStatistic genomic_read_quality_;

  // alignment location information
  int reads_mapped_exonic_ = 0;
  int reads_mapped_intronic_ = 0;
  int reads_mapped_utr_ = 0;

  // in future we can implement this when we have a gene model
  // self.reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
  // self._read_distance_from_termination_site = OnlineGaussianSufficientStatistic()

  // alignment uniqueness information
  int reads_mapped_uniquely_ = 0;
  int reads_mapped_multiple_ = 0;
  int duplicate_reads_ = 0;

  // alignment splicing information
  int spliced_reads_ = 0;
  const int kAntisenseReads = 0; // TODO is never changed from 0
  // int plus_strand_reads_ = 0;  // strand balance (currently unused)




  std::string prev_tag_;
};

class CellMetricGatherer: public MetricGatherer
{
public:
  CellMetricGatherer(std::string gtf_file,
                     std::string mitochondrial_gene_names_filename);

  std::string getHeader() override;
  void output_metrics_extra(std::ofstream& fmetric_out) override;
  void parse_extra_fields(LineFields const& fields) override;

  void clear();

private:
  std::unordered_set<std::string> mitochondrial_genes_;

  int perfect_cell_barcodes_; // The number of reads whose cell barcodes contain no errors (tag ``CB`` == ``CR``)
  int reads_mapped_intergenic_; // The number of reads mapped to an intergenic region for this cell

  int reads_unmapped_;
  //  The number of reads that were mapped to too many loci across the genome and as a
  //  consequence, are reported unmapped by the aligner
  const int kReadsMappedTooManyLoci = 0; // TODO is never changed from 0

  // The variance of the fraction of Illumina base calls for the cell barcode sequence that
  // are greater than 30, across molecules
  float cell_barcode_fraction_bases_above_30_variance;

  // The average fraction of Illumina base calls for the cell barcode sequence that
  // are greater than 30, across molecules
  float cell_barcode_fraction_bases_above_30_mean;

  OnlineGaussianSufficientStatistic cell_barcode_fraction_bases_above_30_;
  std::unordered_map<std::string, int> genes_histogram_;

  std::string cell_specific_headers[11] =
  {
    "perfect_cell_barcodes",
    "reads_mapped_intergenic",
    "reads_unmapped",
    "reads_mapped_too_many_loci",
    "cell_barcode_fraction_bases_above_30_variance",
    "cell_barcode_fraction_bases_above_30_mean",
    "n_genes",
    "genes_detected_multiple_observations",
    "n_mitochondrial_genes",
    "n_mitochondrial_molecules",
    "pct_mitochondrial_molecules"
  };
};


class GeneMetricGatherer: public MetricGatherer
{
public:
  GeneMetricGatherer() : MetricGatherer(MetricType::Gene) {}

  std::string getHeader() override;
  void output_metrics_extra(std::ofstream& fmetric_out) override;
  void parse_extra_fields(LineFields const& fields) override;

  void clear();

private:
  std::unordered_map<std::string, int> cells_histogram_;
  std::string gene_specific_headers[2] =
  {
    "number_cells_detected_multiple",
    "number_cells_expressing"
  };
};

class UmiMetricGatherer: public MetricGatherer // TODO TODO
{
public:
  UmiMetricGatherer() : MetricGatherer(MetricType::Umi) {}

  std::string getHeader() override {return"";}
  void output_metrics_extra(std::ofstream& fmetric_out) override {}
  void parse_extra_fields(LineFields const& fields) override {}
  void clear(){}
};

#endif
