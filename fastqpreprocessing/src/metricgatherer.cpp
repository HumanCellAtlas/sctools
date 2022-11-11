#include "metricgatherer.h"

#include <charconv>

#include "mitochondrial_gene_selector.h"

std::string to_nan(float x)
{
  std::stringstream s;
  s << std::setprecision(10) << x;
  return x==-1 ? "nan" : s.str();
}

LineFields::LineFields(std::string const& s) : s_(s)
{
  first_tag =                           getNextField(); // 0
  second_tag =                          getNextField(); // 1
  third_tag =                           getNextField(); // 2
  reference =                           getNextField(); // 3
  alignment_location =                  getNextField(); // 4
  position_str =                        getNextField(); // 5
  is_strand =                           getNextAsInt(); // 6
  /*std::string_view junk = */          getNextField(); // 7 unused
  cell_barcode_base_above_30 =          getNextAsFloat(); // 8
  genomic_read_quality =                getNextAsFloat(); // 9
  genomic_reads_base_quality_above_30 = getNextAsFloat(); // 10
  number_mappings =                     getNextAsInt(); // 11
  perfect_molecule_barcode =            getNextAsInt(); // 12
  read_spliced =                        getNextAsInt(); // 13
  read_is_duplicate =                   getNextAsInt(); // 14
  cell_barcode_perfect =                getNextAsInt(); // 15
  molecule_barcode_base_above_30 =      getNextAsFloat(); // 16
}

int LineFields::getNextAsInt()
{
  std::string_view sv = getNextField();
  int ret = 0;
  auto [junk, err] { std::from_chars(sv.data(), sv.data() + sv.size(), ret) };
  if (err != std::errc())
    crash(std::string("Not a valid int: ") + std::string(sv));
  return ret;
}
float LineFields::getNextAsFloat()
{
  // TODO replace with from_chars as in getNextAsInt when our compiler has a good
  // from_chars<float> implementation
  return std::stof(std::string(getNextField()));
}
std::string_view LineFields::getNextField()
{
  cur_tab_ = s_.find('\t', cur_start_);
  if (cur_tab_ == std::string::npos)
  {
    if (fields_gotten_ != 16)
      crash("Found " + std::to_string(fields_gotten_+1) + " fields in line; expected 17");

    cur_tab_ = s_.length();
    std::string_view ret(s_.data() + cur_start_, cur_tab_ - cur_start_);
    fields_gotten_++;
    return ret;
  }
  else
  {
    std::string_view ret(s_.data() + cur_start_, cur_tab_ - cur_start_);
    cur_start_ = cur_tab_ + 1;
    fields_gotten_++;
    return ret;
  }
}

MetricGatherer::MetricGatherer(MetricType metric_type) : metric_type_(metric_type) {}

MetricGatherer::~MetricGatherer() {}

void MetricGatherer::clear()
{
  n_reads_ = 0;
  // noise_reads = 0; //# long polymers, N-sequences; NotImplemented
  fragment_histogram_.clear();
  molecule_histogram_.clear();

  molecule_barcode_fraction_bases_above_30_.clear();
  perfect_molecule_barcodes_ = 0;
  genomic_reads_fraction_bases_quality_above_30_.clear();
  genomic_read_quality_.clear();

  reads_mapped_exonic_ = 0;
  reads_mapped_intronic_ = 0;
  reads_mapped_utr_ = 0;

  // alignment uniqueness information
  reads_mapped_uniquely_ = 0;
  reads_mapped_multiple_ = 0;
  duplicate_reads_ = 0;

  // alignment splicing information
  spliced_reads_ = 0;
}

void MetricGatherer::output_metrics(std::ofstream& fmetric_out)
{
  float reads_per_molecule = -1.0f;   // float("nan")
  if (molecule_histogram_.size() != 0)
    reads_per_molecule = n_reads_ / (float)molecule_histogram_.size();

  float reads_per_fragment = -1.0f; //float("nan")
  if (fragment_histogram_.size() != 0)
    reads_per_fragment = n_reads_ / (float)fragment_histogram_.size();

  float fragments_per_molecule = -1.0f; // float("nan")
  if (molecule_histogram_.size() != 0)
    fragments_per_molecule = fragment_histogram_.size() / (float)molecule_histogram_.size();

  int fragments_with_single_read_evidence = 0;
  for (auto const& [key, val] : fragment_histogram_)
    if (val == 1)
      fragments_with_single_read_evidence++;

  int molecules_with_single_read_evidence = 0;
  for (auto const& [key, val] : molecule_histogram_)
    if (val == 1)
      molecules_with_single_read_evidence++;

  fmetric_out << std::setprecision(10) <<  prev_tag_ << ","
              << n_reads_ << ","
              << noise_reads << ","
              << perfect_molecule_barcodes_ << ","
              << reads_mapped_exonic_ << ","
              << reads_mapped_intronic_ << ","
              << reads_mapped_utr_ << ","
              << reads_mapped_uniquely_ << ","
              << reads_mapped_multiple_ << ","
              << duplicate_reads_ << ","
              << spliced_reads_ << ","
              << kAntisenseReads << ","
              << molecule_barcode_fraction_bases_above_30_.getMean() << ","
              << to_nan(molecule_barcode_fraction_bases_above_30_.calculateVariance()) << ","
              << genomic_reads_fraction_bases_quality_above_30_.getMean() << ","
              << to_nan(genomic_reads_fraction_bases_quality_above_30_.calculateVariance()) << ","
              << to_nan(genomic_read_quality_.getMean()) << ","
              << to_nan(genomic_read_quality_.calculateVariance()) << ","
              << molecule_histogram_.size() << ","
              << fragment_histogram_.size() << ","
              << reads_per_molecule << ","
              << reads_per_fragment << ","
              << fragments_per_molecule << ","
              << fragments_with_single_read_evidence << ","
              << molecules_with_single_read_evidence;
}

void MetricGatherer::ingestLine(std::string const& str, std::ofstream& fmetric_out)
{
  LineFields fields(str);

  std::string current_tag(fields.first_tag);
  if (current_tag == "None")
    return; // ignore the None gene

  if (metric_type_ == MetricType::Gene && fields.first_tag.find(",")!=std::string::npos)
    return;

  std::string hyphenated_tags = current_tag + "-";
  hyphenated_tags.append(fields.second_tag);
  hyphenated_tags.append("-");
  hyphenated_tags.append(fields.third_tag);

  // TODO why? this needs a comment explaining why this case makes us start over
  if (prev_tag_ != current_tag && !prev_tag_.empty())
  {
    output_metrics(fmetric_out);
    output_metrics_extra(fmetric_out);
    clear();
  }

  n_reads_++;

  // the tags passed to this function define a molecule, this increments the counter,
  // identifying a new molecule only if a new tag combination is observed

  molecule_histogram_[hyphenated_tags] += 1;

  molecule_barcode_fraction_bases_above_30_.update(fields.molecule_barcode_base_above_30);
  perfect_molecule_barcodes_ += fields.perfect_molecule_barcode;

  genomic_reads_fraction_bases_quality_above_30_.update(fields.genomic_reads_base_quality_above_30);
  genomic_read_quality_.update(fields.genomic_read_quality);

  parse_extra_fields(fields);

  // the remaining portions deal with aligned reads, so if the read is not mapped, we are
  // done with it
  if (fields.reference == "*")
    return;

  parseAlignedReadFields(fields, hyphenated_tags);

  // TODO TODO is it ok to not update this when the aligned thing is '*'?
  prev_tag_ = current_tag;
}

void MetricGatherer::parseAlignedReadFields(LineFields const& fields, std::string hyphenated_tags)
{
  // get components that define a unique sequence fragment and increment the histogram
  std::string is_strand = (fields.is_strand == 1 ? "true" : "false");
  std::string ref_pos_str_tags(fields.reference);
  ref_pos_str_tags.append("\t");
  ref_pos_str_tags.append(fields.position_str);
  ref_pos_str_tags.append("\t");
  ref_pos_str_tags.append(is_strand);
  ref_pos_str_tags.append("\t");
  ref_pos_str_tags.append(hyphenated_tags);
  fragment_histogram_[ref_pos_str_tags] += 1;

  if (fields.alignment_location == "CODING")
    reads_mapped_exonic_ += 1;
  else if (fields.alignment_location == "INTRONIC")
    reads_mapped_intronic_ += 1;
  else if (fields.alignment_location == "UTR")
    reads_mapped_utr_ += 1;

  // in futher check if read maps outside window (when we add a  gene model)
  // and  create distances from terminate side (needs gene model) uniqueness
  if (fields.number_mappings == 1)
    reads_mapped_uniquely_ += 1;
  else
    reads_mapped_multiple_ += 1;  // without multi-mapping, this number is zero!

  duplicate_reads_ += fields.read_is_duplicate;
  spliced_reads_ += fields.read_spliced;
}





////////////////  CellMetricGatherer ////////////////////////

CellMetricGatherer::CellMetricGatherer(std::string gtf_file,
                                       std::string mitochondrial_gene_names_filename)
  : MetricGatherer(MetricType::Cell)
{
  if (gtf_file.empty())
    crash("CellMetricGatherer needs a non-empty gtf_file name!");
  if (mitochondrial_gene_names_filename.empty())
    crash("CellMetricGatherer needs a non-empty mitochondrial_gene_names_filename!");
  mitochondrial_genes_ = getInterestingMitochondrialGenes(
      gtf_file, mitochondrial_gene_names_filename);
}

std::string CellMetricGatherer::getHeader()
{
  std::string s;
  for (int i=0; i<24; i++)
    s += std::string(",") + kCommonHeaders[i]; // TODO ok to start with ,?
  for (int i=0; i<11; i++)
    s += std::string(",") + cell_specific_headers[i];
  return s;
}

// Parses a record to extract gene-specific information
void CellMetricGatherer::parse_extra_fields(LineFields const& fields)
{
  cell_barcode_fraction_bases_above_30_.update(fields.cell_barcode_base_above_30);
  perfect_cell_barcodes_ += fields.cell_barcode_perfect;

  if (!fields.alignment_location.empty())
  {
    if (fields.alignment_location == "INTERGENIC")
      reads_mapped_intergenic_ += 1;
  }
  else
    reads_unmapped_ += 1;

  genes_histogram_[std::string(fields.third_tag)] += 1;
}

void CellMetricGatherer::output_metrics_extra(std::ofstream& fmetric_out)
{
  // The number of genes that are observed by more than one read in this cell
  int genes_detected_multiple_observations = 0;
  // The number of mitochondrial genes detected by this cell
  int n_mitochondrial_genes = 0;
  // The number of molecules from mitochondrial genes detected for this cell
  int n_mitochondrial_molecules = 0;
  float pct_mitochondrial_molecules = 0.0f;

  for (auto const& [gene, count] : genes_histogram_)
  {
    if (count > 1)
      genes_detected_multiple_observations++;
    if (mitochondrial_genes_.find(gene) != mitochondrial_genes_.end())
    {
      n_mitochondrial_genes++;
      n_mitochondrial_molecules += count;
    }
  }

  if (n_mitochondrial_molecules > 0)
  {
    int tot_molecules = 0;
    for (auto const& [gene, count] : genes_histogram_)
      tot_molecules += count;

    // TODO BUG associativity and integer division combine to make this always 0
    pct_mitochondrial_molecules = (n_mitochondrial_molecules/tot_molecules * 100.0f);
  }
  else
    pct_mitochondrial_molecules = 0.0f;

  fmetric_out << std::setprecision(10)
              << "," << perfect_cell_barcodes_
              << "," << reads_mapped_intergenic_
              << "," << reads_unmapped_
              << "," << kReadsMappedTooManyLoci
              << std::setprecision(10)
              << "," << to_nan(cell_barcode_fraction_bases_above_30_.calculateVariance())
              << "," << cell_barcode_fraction_bases_above_30_.getMean()
              << "," << genes_histogram_.size()
              << "," << genes_detected_multiple_observations
              << "," << n_mitochondrial_genes
              << "," << n_mitochondrial_molecules
              << "," << pct_mitochondrial_molecules
              << std::endl;
}

void CellMetricGatherer::clear()
{
  MetricGatherer::clear();

  cell_barcode_fraction_bases_above_30_.clear();
  perfect_cell_barcodes_ = 0;
  reads_mapped_intergenic_ = 0;
  reads_unmapped_ = 0;
  genes_histogram_.clear();
}


////////////////  GeneMetricGatherer ////////////////////////
std::string GeneMetricGatherer::getHeader()
{
  std::string s;
  for (int i=0; i<24; i++)
    s += std::string(",") + kCommonHeaders[i]; // TODO ok to start with ,?
  for (int i=0; i<2; i++)
    s += std::string(",") + gene_specific_headers[i];
  return s;
}

void GeneMetricGatherer::parse_extra_fields(LineFields const& fields)
{
  cells_histogram_[std::string(fields.second_tag)] += 1;
}

void GeneMetricGatherer::output_metrics_extra(std::ofstream& fmetric_out)
{
  int number_cells_detected_multiple = 0;
  for (auto const& [cell, count] : cells_histogram_)
    if (count > 1)
      number_cells_detected_multiple++;

  fmetric_out <<  ","  << number_cells_detected_multiple
              <<  ","  << cells_histogram_.size()
              << std::endl;
}

void GeneMetricGatherer::clear()
{
  MetricGatherer::clear();
  cells_histogram_.clear();
}
