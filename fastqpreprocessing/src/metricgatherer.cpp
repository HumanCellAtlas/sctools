/**
 *  @file   metricgatherer.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/
#include "metricgatherer.h"

#define ASSERT(condition) { if(!(condition)){ std::cerr << "ASSERT FAILED: " << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; } }
using namespace std;

namespace consts
{
const std::string INTRONIC_ALIGNMENT_LOCATION_TAG_VALUE = "INTRONIC";
const std::string CODING_ALIGNMENT_LOCATION_TAG_VALUE = "CODING";
const std::string UTR_ALIGNMENT_LOCATION_TAG_VALUE = "UTR";
const std::string INTERGENIC_ALIGNMENT_LOCATION_TAG_VALUE = "INTERGENIC";
}

template<typename T>
inline void freeContainer(T& p_container)
{
  T empty;
  using std::swap;
  swap(p_container, empty);
}


std::string to_nan(float x)
{
  stringstream s;
  s << std::setprecision(10) << x;
  return x==-1 ? "nan" : s.str();
}

void Metrics::clear()
{
  n_reads = 0;
  noise_reads = 0; //# long polymers, N-sequences; NotImplemented
  freeContainer(_fragment_histogram);
  freeContainer(_molecule_histogram);

  freeContainer(_molecule_barcode_fraction_bases_above_30);
  perfect_molecule_barcodes = 0;
  freeContainer(_genomic_reads_fraction_bases_quality_above_30);
  freeContainer(_genomic_read_quality);

  reads_mapped_exonic = 0;
  reads_mapped_intronic = 0;
  reads_mapped_utr = 0;

  // alignment uniqueness information
  reads_mapped_uniquely = 0;
  reads_mapped_multiple = 0;
  duplicate_reads = 0;

  // alignment splicing information
  spliced_reads = 0;
  antisense_reads = 0;
  plus_strand_reads = 0;  // strand balance

  // higher-order methods, filled in by finalize() when all data is extracted
  molecule_barcode_fraction_bases_above_30_mean = -1;
  molecule_barcode_fraction_bases_above_30_variance = -1;
  genomic_reads_fraction_bases_quality_above_30_mean = -1;
  genomic_reads_fraction_bases_quality_above_30_variance = -1;
  genomic_read_quality_mean = -1;
  genomic_read_quality_variance  = -1;
  n_molecules = -1;
  n_fragments = -1;
  reads_per_molecule = -1;
  reads_per_fragment = -1;
  fragments_per_molecule = -1;
  fragments_with_single_read_evidence = -1;
  molecules_with_single_read_evidence = -1;
}


void Metrics::output_metrics(ofstream& fmetric_out)
{
  fmetric_out << std::setprecision(10) <<  prev_tag << ","
              << this->n_reads << ","
              << this->noise_reads << ","
              << this->perfect_molecule_barcodes << ","
              << this->reads_mapped_exonic << ","
              << this->reads_mapped_intronic << ","
              << this->reads_mapped_utr << ","
              << this->reads_mapped_uniquely << ","
              << this->reads_mapped_multiple << ","
              << this->duplicate_reads << ","
              << this->spliced_reads << ","
              << this->antisense_reads << ","
              << this->molecule_barcode_fraction_bases_above_30_mean << ","
              << to_nan(this->molecule_barcode_fraction_bases_above_30_variance) << ","
              << this->genomic_reads_fraction_bases_quality_above_30_mean << ","
              << to_nan(this->genomic_reads_fraction_bases_quality_above_30_variance) << ","
              << to_nan(this->genomic_read_quality_mean) << ","
              << to_nan(this->genomic_read_quality_variance) << ","
              << this->n_molecules << ","
              << this->n_fragments << ","
              << this->reads_per_molecule << ","
              << this->reads_per_fragment << ","
              << this->fragments_per_molecule << ","
              << this->fragments_with_single_read_evidence << ","
              << this->molecules_with_single_read_evidence;
}


unsigned int offset = 3;
void Metrics::parse_line(std::string& str, ofstream& fmetric_out,
                         std::set<std::string>& mitochondrial_genes,
                         METRIC_TYPE metric_type)
{
  char line[1000];
  std::string NONE("None");
  std::size_t len = str.copy(line, str.size(), 0);
  line[str.size()]='\0';

  assert(len < 1000);

  char* c = line;

  unsigned int k = 0;
  record[k] = c;
  while (*c!='\0')
  {
    if (*c == '\t')
    {
      *c='\0';
      record[++k] = c + 1;
    }
    c++;
  }

  assert(k==16);

  std::string current_tag = record[0];

  // ignore the None gene
  if (current_tag.compare(NONE)==0)
    return;

  // load the tags
  std::string first_tag, second_tag, third_tag;
  first_tag = record[0];
  second_tag = record[1];
  third_tag = record[2];
  if (metric_type == GENE && first_tag.find(",")!=std::string::npos)
    return;

  std::string tags = first_tag + std::string("-") + second_tag + std::string("-") + third_tag;
  if (prev_tag.compare(current_tag)!=0 && prev_tag.size()!=0)
  {
    this->finalize(mitochondrial_genes);
    this->output_metrics(fmetric_out);
    this->output_metrics_extra(fmetric_out);
    this->clear();
  }

  this->parse_extra_fields(first_tag, second_tag, third_tag, record);

  this->n_reads += 1;

  // the tags passed to this function define a molecule, this increments the counter,
  // identifying a new molecule only if a new tag combination is observed

  /* updating the molecule histogram with tags */
  if (this->_molecule_histogram.find(tags)==this->_molecule_histogram.end())
    this->_molecule_histogram[tags] = 0;
  this->_molecule_histogram[tags] += 1;

  this->_molecule_barcode_fraction_bases_above_30.update(std::stof(record[offset + 13]));
  this->perfect_molecule_barcodes += std::stoi(record[offset + 9]);

  this->_genomic_reads_fraction_bases_quality_above_30.update(std::stof(record[offset + 7]));
  this->_genomic_read_quality.update(std::stof(record[offset + 6]));

  // the remaining portions deal with aligned reads, so if the read is not mapped, we are
  // done with it
  if (std::string(record[offset + 0]).compare("*")==0)
    return;

  // get components that define a unique sequence fragment and increment the histogram
  std::string position_str = record[offset + 2];
  std::string strand = std::stoi(std::string(record[offset + 3]))==1 ? "true" : "false";
  string reference = record[offset + 0];

  std::string _ref_pos_str_tags = reference + std::string("\t") +
                                  position_str + std::string("\t") +
                                  strand + std::string("\t") + tags;
  std::string ref_pos_str_tags = std::to_string(std::hash<std::string> {}(_ref_pos_str_tags));

  /* updating the fragment histogram with tag, strand and pos */
  if (this->_fragment_histogram.find(ref_pos_str_tags)==this->_fragment_histogram.end())
    this->_fragment_histogram[ref_pos_str_tags] = 0;
  this->_fragment_histogram[ref_pos_str_tags] += 1;

  std::string alignment_location = std::string(record[offset + 1]);
  if (alignment_location.compare(consts::CODING_ALIGNMENT_LOCATION_TAG_VALUE)==0)
    reads_mapped_exonic += 1;
  else if (alignment_location.compare(consts::INTRONIC_ALIGNMENT_LOCATION_TAG_VALUE)==0)
    reads_mapped_intronic += 1;
  else if (alignment_location.compare(consts::UTR_ALIGNMENT_LOCATION_TAG_VALUE)==0)
    reads_mapped_utr += 1;

  // in futher check if read maps outside window (when we add a  gene model)
  // and  create distances from terminate side (needs gene model) uniqueness
  int number_mappings = std::stoi(std::string(record[offset + 8]));

  if (number_mappings==1)
    this->reads_mapped_uniquely += 1;
  else
    this->reads_mapped_multiple += 1;  // without multi-mapping, this number is zero!

  this->duplicate_reads += std::stoi(std::string(record[offset + 11]));

  // cigar N field (3) indicates a read is spliced if the value is non-zero
  this->spliced_reads += std::stoi(std::string(record[offset + 10]));

  prev_tag = current_tag;
}


//  Calculate metrics that require information from all molecules of an entity
//  ``finalize()`` replaces attributes in-place that were initialized by the constructor as
//  ``None`` with a value calculated across all molecule data that has been aggregated.

void Metrics::finalize(std::set<std::string>& mitochondrial_genes)
{
  this->molecule_barcode_fraction_bases_above_30_mean =
      this->_molecule_barcode_fraction_bases_above_30.getMean();

  this->molecule_barcode_fraction_bases_above_30_variance =
      this->_molecule_barcode_fraction_bases_above_30.calculate_variance();

  this->genomic_reads_fraction_bases_quality_above_30_mean =
      this->_genomic_reads_fraction_bases_quality_above_30.getMean();

  this->genomic_reads_fraction_bases_quality_above_30_variance =
      this->_genomic_reads_fraction_bases_quality_above_30.calculate_variance();

  this->genomic_read_quality_mean = this->_genomic_read_quality.getMean();

  this->genomic_read_quality_variance = this->_genomic_read_quality.calculate_variance();

  this->n_molecules = this->_molecule_histogram.size();

  this->n_fragments = this->_fragment_histogram.size();

  try
  {
    this->reads_per_molecule = this->n_reads / this->n_molecules;
  }
  catch (runtime_error& e)
  {
    this->reads_per_molecule = -1;   // float("nan")
  }

  try
  {
    this->reads_per_fragment = this->n_reads / this->n_fragments;
  }
  catch (runtime_error& e)
  {
    this->reads_per_fragment = -1; //float("nan")
  }

  try
  {
    this->fragments_per_molecule = this->n_fragments / this->n_molecules;
  }
  catch (runtime_error& e)
  {
    this->fragments_per_molecule = -1; // float("nan")
  }

  this->fragments_with_single_read_evidence = 0;
  for (auto local_cit = this->_fragment_histogram.cbegin();
       local_cit != this->_fragment_histogram.cend();
       local_cit++)
  {
    if (local_cit->second==1)
      this->fragments_with_single_read_evidence++;
  }

  this->molecules_with_single_read_evidence = 0;
  for (auto local_cit = this->_molecule_histogram.cbegin();
       local_cit != this->_molecule_histogram.cend();
       local_cit++)
  {
    if (local_cit->second==1)
      this->molecules_with_single_read_evidence++;
  }
}

////////////////  CellMetrics ////////////////////////
std::string CellMetrics::getHeader()
{
  std::string s;
  for (int i=0; i<24; i++)
    s += std::string(",") + headers[i];
  for (int i=0; i<11; i++)
    s += std::string(",") + specific_headers[i];
  return s;
}

// Parses a record to extract gene-specific information
void CellMetrics::parse_extra_fields(const std::string& first_tag,
                                     const std::string& second_tag,
                                     const std::string& third_tag,
                                     char** record)
{
  this->_cell_barcode_fraction_bases_above_30.update(std::stof(record[offset + 5]));
  this->perfect_cell_barcodes += std::stoi(record[offset + 12]);

  if (std::string(record[offset + 1]).size() != 0)
  {
    if (std::string(record[offset + 1]).compare(consts::INTERGENIC_ALIGNMENT_LOCATION_TAG_VALUE)==0)
      this->reads_mapped_intergenic += 1;
  }
  else
    this->reads_unmapped += 1;

  /* updating the genes histogram with tags */
  if (this->_genes_histogram.find(third_tag) == this->_genes_histogram.end())
    this->_genes_histogram[third_tag] = 0;
  this->_genes_histogram[third_tag] += 1;
}

void CellMetrics::output_metrics_extra(ofstream& fmetric_out)
{
  fmetric_out << std::setprecision(10)
              << "," << this->perfect_cell_barcodes
              << "," << this->reads_mapped_intergenic
              << "," << this->reads_unmapped
              << "," << this->reads_mapped_too_many_loci
              << std::setprecision(10)
              << "," << to_nan(this->cell_barcode_fraction_bases_above_30_variance)
              << "," << this->cell_barcode_fraction_bases_above_30_mean
              << "," << this->n_genes
              << "," << this->genes_detected_multiple_observations
              << "," << this->n_mitochondrial_genes
              << "," << this->n_mitochondrial_molecules
              << "," << this->pct_mitochondrial_molecules
              << std::endl;
}

void CellMetrics::finalize(std::set<std::string>& mitochondrial_genes)
{
  // call the finalize function in the parent class
  Metrics::finalize(mitochondrial_genes);

  this->cell_barcode_fraction_bases_above_30_mean =
      this->_cell_barcode_fraction_bases_above_30.getMean();

  this->cell_barcode_fraction_bases_above_30_variance =
      this->_cell_barcode_fraction_bases_above_30.calculate_variance();

  this->n_genes = this->_genes_histogram.size();

  this->genes_detected_multiple_observations = 0;
  for (auto local_cit = this->_genes_histogram.cbegin();
       local_cit != this->_genes_histogram.cend();
       local_cit++)
  {
    if (local_cit->second > 1)
      this->genes_detected_multiple_observations++;
  }

  this->n_mitochondrial_genes = 0;
  for (auto local_cit = this->_genes_histogram.cbegin();
       local_cit!=this->_genes_histogram.cend();
       local_cit++)
  {
    if (mitochondrial_genes.find(local_cit->first) != mitochondrial_genes.end())
      this->n_mitochondrial_genes++;
  }

  this->n_mitochondrial_molecules = 0;
  for (auto local_cit = this->_genes_histogram.cbegin();
       local_cit != this->_genes_histogram.cend();
       local_cit++)
  {
    if (mitochondrial_genes.find(local_cit->first)!=mitochondrial_genes.end())
      this->n_mitochondrial_molecules += local_cit->second;
  }

  if (this->n_mitochondrial_molecules > 0)
  {
    int tot_molecules = 0;
    for (auto local_cit = this->_genes_histogram.cbegin();
         local_cit != this->_genes_histogram.cend();
         local_cit++)
    {
      tot_molecules += local_cit->second;
    }
    this->pct_mitochondrial_molecules = (this->n_mitochondrial_molecules/tot_molecules * 100.0);
  }
  else
    this->pct_mitochondrial_molecules = 0.0;
}

void CellMetrics::clear()
{
  // call the clear function in the parent class
  Metrics::clear();

  this->_cell_barcode_fraction_bases_above_30.clear();
  this->perfect_cell_barcodes = 0;
  this->reads_mapped_intergenic = 0;
  this->reads_unmapped = 0;
  this->reads_mapped_too_many_loci = 0;
  this->_genes_histogram.clear();

  this->n_genes = 0;
  this->genes_detected_multiple_observations = 0;
  this->n_mitochondrial_genes = 0;
  this->n_mitochondrial_molecules = 0;
  this->pct_mitochondrial_molecules = 0;
}


////////////////  GeneMetrics ////////////////////////
std::string GeneMetrics::getHeader()
{
  std::string s;
  for (int i=0; i<24; i++)
    s += std::string(",") + headers[i];
  for (int i=0; i<2; i++)
    s += std::string(",") + specific_headers[i];
  return s;
}

void GeneMetrics::parse_extra_fields(const std::string& first_tag,
                                     const std::string& second_tag,
                                     const std::string& third_tag,
                                     char** record)
{
  // updating the cell histogram with tags
  if (this->_cells_histogram.find(second_tag)==this->_cells_histogram.end())
    this->_cells_histogram[second_tag] = 0;

  this->_cells_histogram[second_tag] += 1;
}

void GeneMetrics::output_metrics_extra(ofstream& fmetric_out)
{
  fmetric_out <<  ","  << number_cells_detected_multiple
              <<  ","  << number_cells_expressing
              << std::endl;
}

void GeneMetrics::finalize(std::set<std::string>& mitochondrial_genes)
{
  // call the finalize function in the parent class
  Metrics::finalize(mitochondrial_genes);
  this->number_cells_expressing = this->_cells_histogram.size();

  this->number_cells_detected_multiple = 0;
  for (auto local_cit = this->_cells_histogram.cbegin();
       local_cit != this->_cells_histogram.cend();
       local_cit++)
  {
    if (local_cit->second > 1)
      this->number_cells_detected_multiple++;
  }
}

void GeneMetrics::clear()
{
  // call the clear function in the parent class
  Metrics::clear();
  number_cells_detected_multiple = 0;
  number_cells_expressing = 0;
  freeContainer(_cells_histogram);
}
