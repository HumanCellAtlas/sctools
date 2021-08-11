#ifndef __METRIC_GATHERER__
#define __METRIC_GATHERER__
#include <unordered_map>
#include <map>
#include <string>
#include <regex>
#include <iostream>
#include <vector>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "datatypes.h"

using namespace std;

//    Implementation of Welford's online mean and variance algorithm
//
//    Methods
//    -------
//    update(new_value: float)
//        incorporate new_value into the online estimate of mean and variance
//    getMean()
//        return the mean value
//    calculate_variance()
//        calculate and return the variance
//    mean_and_variance()
//        return both mean and variance
class OnlineGaussianSufficientStatistic {
    private:
        double _mean_squared_error = 0.0;
        double sum_EX2 = 0.0;
        double _mean = 0.0;
        double _sum = 0.0;
        double _count = 0.0;

    public:
      void update(double new_value) {
         _count += 1.0;
         _sum += new_value;
         sum_EX2 += (new_value*new_value);
      }

      // return the mean value
      double getMean() {
         _mean = _sum/_count;
         return _mean;
      }

      // calculate and return the variance
      double calculate_variance() {
        if (_count < 2) {
           //return double("nan")
          return -1.0;
        }
        return  sum_EX2/(_count -1) - (_sum/_count)*(_sum/(_count -1));
      }

      void clear() {
        _mean_squared_error = 0.0;
        _mean = 0.0;
        _count = 0;
        _sum = 0;
        sum_EX2 = 0.0;
      }
};

class Metrics {
  private:
    // count information
    int n_reads = 0;
    int noise_reads = 0; //# long polymers, N-sequences; NotImplemented

    //std::unordered_map<std::string, int> _fragment_histogram;
    std::map<std::string, int> _fragment_histogram;

    //self._molecule_histogram: Counter[str] = Counter()
    std::unordered_map<std::string, int> _molecule_histogram;

    // molecule information
    OnlineGaussianSufficientStatistic _molecule_barcode_fraction_bases_above_30;

    int perfect_molecule_barcodes = 0;

    OnlineGaussianSufficientStatistic _genomic_reads_fraction_bases_quality_above_30;

    OnlineGaussianSufficientStatistic _genomic_read_quality;

    // alignment location information
    int reads_mapped_exonic = 0;
    int reads_mapped_intronic = 0;
    int reads_mapped_utr = 0;

    // todo implement this once we have a gene model
    // self.reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
    // self._read_distance_from_termination_site = OnlineGaussianSufficientStatistic()

    // alignment uniqueness information
    int reads_mapped_uniquely = 0;
    int reads_mapped_multiple = 0;
    int duplicate_reads = 0;

    // alignment splicing information
    int spliced_reads = 0;
    int antisense_reads = 0;
    int plus_strand_reads = 0;  // strand balance  # todo implement property here

    // higher-order methods, filled in by finalize() when all data is extracted
    float molecule_barcode_fraction_bases_above_30_mean = -1;
    float molecule_barcode_fraction_bases_above_30_variance = -1;
    float genomic_reads_fraction_bases_quality_above_30_mean = -1;
    float genomic_reads_fraction_bases_quality_above_30_variance = -1;
    float genomic_read_quality_mean = -1;
    float genomic_read_quality_variance  = -1;
    float n_molecules = -1;
    float n_fragments = -1;
    float reads_per_molecule = -1;
    float reads_per_fragment = -1;
    float fragments_per_molecule = -1;
    int fragments_with_single_read_evidence = -1;
    int molecules_with_single_read_evidence = -1;
  private:
    std::regex rgx;
    std::sregex_token_iterator end;
    std::string prev_tag;
    char *record[20]; 

  protected:
    std::string headers[24] = {
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


  public:
    Metrics() {
      rgx.assign("\t");
      prev_tag = std::string("");
    }
    virtual ~Metrics() {}
    //  get the headers
    virtual std::string getHeader() { return std::string(""); };

    void parse_line(std::string &str, ofstream &fmetric_out, 
            std::set<std::string> &mitochondrial_genes,
            METRIC_TYPE metric_type);

    void output_metrics(ofstream &fmetric_out);
    virtual void output_metrics_extra(ofstream &fmetric_out) {}
    virtual void parse_extra_fields(const std::string &first_tag, 
                            const std::string &second_tag,
                            const std::string &third_tag,
                            char** record) {}
    virtual void finalize(std::set<std::string> &mitochondrial_genes);
    virtual void clear();
};

class CellMetrics: public Metrics {
  private:
    int perfect_cell_barcodes; // The number of reads whose cell barcodes contain no errors (tag ``CB`` == ``CR``)
    int reads_mapped_intergenic; // The number of reads mapped to an intergenic region for this cell

    // reads unmapped
    int reads_unmapped;
    //  The number of reads that were mapped to too many loci across the genome and as a
    //  consequence, are reported unmapped by the aligner
    int reads_mapped_too_many_loci;

    // The variance of the fraction of Illumina base calls for the cell barcode sequence that
    // are greater than 30, across molecules
    float cell_barcode_fraction_bases_above_30_variance;

    // The average fraction of Illumina base calls for the cell barcode sequence that
    // are greater than 30, across molecules
    float cell_barcode_fraction_bases_above_30_mean;

    int n_genes;  //The number of genes detected by this cell

    int genes_detected_multiple_observations; // The number of genes that are observed by more than one read in this cell
    int n_mitochondrial_genes; // The number of mitochondrial genes detected by this cell
    int n_mitochondrial_molecules; // The number of molecules from mitochondrial genes detected for this cell
    int pct_mitochondrial_molecules; // The percentage of molecules from mitoc
    
    OnlineGaussianSufficientStatistic _cell_barcode_fraction_bases_above_30;
    std::unordered_map<std::string, int> _genes_histogram;

    std::string specific_headers[11] = {
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

  public:
    std::string getHeader();
    void output_metrics_extra(ofstream &fmetric_out);
    void parse_extra_fields(const std::string &first_tag, 
                            const std::string &second_tag,
                            const std::string &third_tag,
                            char** record);

    void finalize(std::set<std::string> &mitochondrial_genes);

    void clear();
};


class GeneMetrics: public Metrics {
  private:
    int number_cells_detected_multiple;
    int number_cells_expressing;

    std::unordered_map<std::string, int> _cells_histogram;
    std::string specific_headers[2] = {
                 "number_cells_detected_multiple",  
                 "number_cells_expressing" 
               };
 
  public:
    GeneMetrics() {
       number_cells_detected_multiple = 0;
       number_cells_expressing = 0;
    }

  public:
    std::string getHeader();
    void output_metrics_extra(ofstream &fmetric_out);
    void parse_extra_fields(const std::string &first_tag, 
                            const std::string &second_tag,
                            const std::string &third_tag,
                            char** record);

    void finalize(std::set<std::string> &mitochondrial_genes);
    void clear();
};




#endif 
