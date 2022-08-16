/**
 *  @file   metricgatherer.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/
#include "metricgatherer.h"

template<typename T>
inline void freeContainer(T& p_container)
{
  T empty;
  using std::swap;
  swap(p_container, empty);
}


std::string to_nan(float x)
{
  std::stringstream s;
  s << std::setprecision(10) << x;
  return x==-1 ? "nan" : s.str();
}

void Metrics::clear()
{
  n_reads = 0;
  // noise_reads = 0; //# long polymers, N-sequences; NotImplemented
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


void Metrics::output_metrics(std::ofstream& fmetric_out)
{
  fmetric_out << std::setprecision(10) <<  prev_tag << ","
              << n_reads << ","
              << noise_reads << ","
              << perfect_molecule_barcodes << ","
              << reads_mapped_exonic << ","
              << reads_mapped_intronic << ","
              << reads_mapped_utr << ","
              << reads_mapped_uniquely << ","
              << reads_mapped_multiple << ","
              << duplicate_reads << ","
              << spliced_reads << ","
              << antisense_reads << ","
              << molecule_barcode_fraction_bases_above_30_mean << ","
              << to_nan(molecule_barcode_fraction_bases_above_30_variance) << ","
              << genomic_reads_fraction_bases_quality_above_30_mean << ","
              << to_nan(genomic_reads_fraction_bases_quality_above_30_variance) << ","
              << to_nan(genomic_read_quality_mean) << ","
              << to_nan(genomic_read_quality_variance) << ","
              << n_molecules << ","
              << n_fragments << ","
              << reads_per_molecule << ","
              << reads_per_fragment << ","
              << fragments_per_molecule << ","
              << fragments_with_single_read_evidence << ","
              << molecules_with_single_read_evidence;
}


unsigned int offset = 3;
void Metrics::parse_line(std::string& str, std::ofstream& fmetric_out,
                         std::unordered_set<std::string>& mitochondrial_genes,
                         MetricType metric_type)
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
  if (metric_type == MetricType::Gene && first_tag.find(",")!=std::string::npos)
    return;

  std::string tags = first_tag + std::string("-") + second_tag + std::string("-") + third_tag;
  if (prev_tag.compare(current_tag)!=0 && prev_tag.size()!=0)
  {
    finalize(mitochondrial_genes);
    output_metrics(fmetric_out);
    output_metrics_extra(fmetric_out);
    clear();
  }

  parse_extra_fields(first_tag, second_tag, third_tag, record);

  n_reads += 1;

  // the tags passed to this function define a molecule, this increments the counter,
  // identifying a new molecule only if a new tag combination is observed

  /* updating the molecule histogram with tags */
  if (_molecule_histogram.find(tags)==_molecule_histogram.end())
    _molecule_histogram[tags] = 0;
  _molecule_histogram[tags] += 1;

  _molecule_barcode_fraction_bases_above_30.update(std::stof(record[offset + 13]));
  perfect_molecule_barcodes += std::stoi(record[offset + 9]);

  _genomic_reads_fraction_bases_quality_above_30.update(std::stof(record[offset + 7]));
  _genomic_read_quality.update(std::stof(record[offset + 6]));

  // the remaining portions deal with aligned reads, so if the read is not mapped, we are
  // done with it
  if (std::string(record[offset + 0]).compare("*")==0)
    return;

  // get components that define a unique sequence fragment and increment the histogram
  std::string position_str = record[offset + 2];
  std::string strand = std::stoi(std::string(record[offset + 3]))==1 ? "true" : "false";
  std::string reference = record[offset + 0];

  std::string _ref_pos_str_tags = reference + std::string("\t") +
                                  position_str + std::string("\t") +
                                  strand + std::string("\t") + tags;
  std::string ref_pos_str_tags = std::to_string(std::hash<std::string> {}(_ref_pos_str_tags));

  /* updating the fragment histogram with tag, strand and pos */
  if (_fragment_histogram.find(ref_pos_str_tags)==_fragment_histogram.end())
    _fragment_histogram[ref_pos_str_tags] = 0;
  _fragment_histogram[ref_pos_str_tags] += 1;

  std::string alignment_location = std::string(record[offset + 1]);
  if (alignment_location == "CODING")
    reads_mapped_exonic += 1;
  else if (alignment_location == "INTRONIC")
    reads_mapped_intronic += 1;
  else if (alignment_location == "UTR")
    reads_mapped_utr += 1;

  // in futher check if read maps outside window (when we add a  gene model)
  // and  create distances from terminate side (needs gene model) uniqueness
  int number_mappings = std::stoi(std::string(record[offset + 8]));

  if (number_mappings==1)
    reads_mapped_uniquely += 1;
  else
    reads_mapped_multiple += 1;  // without multi-mapping, this number is zero!

  duplicate_reads += std::stoi(std::string(record[offset + 11]));

  // cigar N field (3) indicates a read is spliced if the value is non-zero
  spliced_reads += std::stoi(std::string(record[offset + 10]));

  prev_tag = current_tag;
}


//  Calculate metrics that require information from all molecules of an entity
//  ``finalize()`` replaces attributes in-place that were initialized by the constructor as
//  ``None`` with a value calculated across all molecule data that has been aggregated.

void Metrics::finalize(std::unordered_set<std::string>& mitochondrial_genes)
{
  molecule_barcode_fraction_bases_above_30_mean =
      _molecule_barcode_fraction_bases_above_30.getMean();

  molecule_barcode_fraction_bases_above_30_variance =
      _molecule_barcode_fraction_bases_above_30.calculate_variance();

  genomic_reads_fraction_bases_quality_above_30_mean =
      _genomic_reads_fraction_bases_quality_above_30.getMean();

  genomic_reads_fraction_bases_quality_above_30_variance =
      _genomic_reads_fraction_bases_quality_above_30.calculate_variance();

  genomic_read_quality_mean = _genomic_read_quality.getMean();

  genomic_read_quality_variance = _genomic_read_quality.calculate_variance();

  n_molecules = _molecule_histogram.size();

  n_fragments = _fragment_histogram.size();

  reads_per_molecule = -1;   // float("nan")
  if (n_molecules != 0)
    reads_per_molecule = n_reads / n_molecules;

  reads_per_fragment = -1; //float("nan")
  if (n_fragments != 0)
    reads_per_fragment = n_reads / n_fragments;

  fragments_per_molecule = -1; // float("nan")
  if (n_molecules != 0)
    fragments_per_molecule = n_fragments / n_molecules;

  fragments_with_single_read_evidence = 0;
  for (auto const& [key, val] : _fragment_histogram)
    if (val == 1)
      fragments_with_single_read_evidence++;

  molecules_with_single_read_evidence = 0;
  for (auto const& [key, val] : _molecule_histogram)
    if (val == 1)
      molecules_with_single_read_evidence++;
}

////////////////  CellMetrics ////////////////////////
std::string CellMetrics::getHeader()
{
  std::string s;
  for (int i=0; i<24; i++)
    s += std::string(",") + common_headers[i]; // TODO ok to start with ,?
  for (int i=0; i<11; i++)
    s += std::string(",") + cell_specific_headers[i];
  return s;
}

// Parses a record to extract gene-specific information
void CellMetrics::parse_extra_fields(std::string const& first_tag,
                                     std::string const& second_tag,
                                     std::string const& third_tag,
                                     char** record)
{
  _cell_barcode_fraction_bases_above_30.update(std::stof(record[offset + 5]));
  perfect_cell_barcodes += std::stoi(record[offset + 12]);

  std::string record_str(record[offset + 1]);
  if (!record_str.empty()) // TODO can the empty check be skipped, or is there a non-empty non-unmapped case?
  {
    if (record_str == "INTERGENIC")
      reads_mapped_intergenic += 1;
  }
  else
    reads_unmapped += 1;

  /* updating the genes histogram with tags */
  if (_genes_histogram.find(third_tag) == _genes_histogram.end())
    _genes_histogram[third_tag] = 0;
  _genes_histogram[third_tag] += 1;
}

void CellMetrics::output_metrics_extra(std::ofstream& fmetric_out)
{
  fmetric_out << std::setprecision(10)
              << "," << perfect_cell_barcodes
              << "," << reads_mapped_intergenic
              << "," << reads_unmapped
              << "," << reads_mapped_too_many_loci
              << std::setprecision(10)
              << "," << to_nan(cell_barcode_fraction_bases_above_30_variance)
              << "," << cell_barcode_fraction_bases_above_30_mean
              << "," << n_genes
              << "," << genes_detected_multiple_observations
              << "," << n_mitochondrial_genes
              << "," << n_mitochondrial_molecules
              << "," << pct_mitochondrial_molecules
              << std::endl;
}

void CellMetrics::finalize(std::unordered_set<std::string>& mitochondrial_genes)
{
  // call the finalize function in the parent class
  Metrics::finalize(mitochondrial_genes);

  cell_barcode_fraction_bases_above_30_mean =
      _cell_barcode_fraction_bases_above_30.getMean();

  cell_barcode_fraction_bases_above_30_variance =
      _cell_barcode_fraction_bases_above_30.calculate_variance();

  n_genes = _genes_histogram.size();

  genes_detected_multiple_observations = 0;
  n_mitochondrial_genes = 0;
  n_mitochondrial_molecules = 0;
  for (auto const& [gene, count] : _genes_histogram)
  {
    if (count > 1)
      genes_detected_multiple_observations++;
    if (mitochondrial_genes.find(gene) != mitochondrial_genes.end())
    {
      n_mitochondrial_genes++;
      n_mitochondrial_molecules += count;
    }
  }

  if (n_mitochondrial_molecules > 0)
  {
    int tot_molecules = 0;
    for (auto const& [gene, count] : _genes_histogram)
      tot_molecules += count;

    // TODO BUG associativity and integer division combine to make this always 0
    pct_mitochondrial_molecules = (n_mitochondrial_molecules/tot_molecules * 100.0);
  }
  else
    pct_mitochondrial_molecules = 0.0;
}

void CellMetrics::clear()
{
  // call the clear function in the parent class
  Metrics::clear();

  _cell_barcode_fraction_bases_above_30.clear();
  perfect_cell_barcodes = 0;
  reads_mapped_intergenic = 0;
  reads_unmapped = 0;
  reads_mapped_too_many_loci = 0;
  _genes_histogram.clear();

  n_genes = 0;
  genes_detected_multiple_observations = 0;
  n_mitochondrial_genes = 0;
  n_mitochondrial_molecules = 0;
  pct_mitochondrial_molecules = 0;
}


////////////////  GeneMetrics ////////////////////////
std::string GeneMetrics::getHeader()
{
  std::string s;
  for (int i=0; i<24; i++)
    s += std::string(",") + common_headers[i]; // TODO ok to start with ,?
  for (int i=0; i<2; i++)
    s += std::string(",") + gene_specific_headers[i];
  return s;
}

void GeneMetrics::parse_extra_fields(std::string const& first_tag,
                                     std::string const& second_tag,
                                     std::string const& third_tag,
                                     char** record)
{
  // updating the cell histogram with tags
  if (_cells_histogram.find(second_tag)==_cells_histogram.end())
    _cells_histogram[second_tag] = 0;

  _cells_histogram[second_tag] += 1;
}

void GeneMetrics::output_metrics_extra(std::ofstream& fmetric_out)
{
  fmetric_out <<  ","  << number_cells_detected_multiple
              <<  ","  << number_cells_expressing
              << std::endl;
}

void GeneMetrics::finalize(std::unordered_set<std::string>& mitochondrial_genes)
{
  // call the finalize function in the parent class
  Metrics::finalize(mitochondrial_genes);

  number_cells_expressing = _cells_histogram.size();
  number_cells_detected_multiple = 0;
  for (auto const& [cell, count] : _cells_histogram)
    if (count > 1)
      number_cells_detected_multiple++;
}

void GeneMetrics::clear()
{
  // call the clear function in the parent class
  Metrics::clear();
  number_cells_detected_multiple = 0;
  number_cells_expressing = 0;
  freeContainer(_cells_histogram);
}
