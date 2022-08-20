/**
 *  @file   tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori M. Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <fstream>
#include <queue>
#include <regex>
#include <string>
#include <sstream>
#include <unordered_set>

#include "htslib_tagsort.h"
#include "metricgatherer.h"

constexpr int kDataBufferSize = 1000;

struct Context
{
  std::vector<std::vector<std::string>> data;

  std::vector<long int> file_offset;
  std::vector<int> data_size;
  std::vector<int> ptrs;
  std::vector<bool> isempty;
  int index_ = -1;
  int num_active_files = 0;
  const int num_parts_;

  Context(unsigned int num_parts) : num_parts_(num_parts)
  {
    // set the file offsets to 0
    for (int i=0; i < num_parts_; i++)
      file_offset.push_back(0);

    // set the isempty for each file to false
    for (int i=0; i < num_parts_; i++)
      isempty.push_back(false);

    // set a vector of vectors of data for each file
    for (int i=0; i < num_parts_; i++)
      data.push_back(std::vector<std::string>());

    // set the data_size of the buffer for each file to 0
    for (int i=0; i < num_parts_; i++)
      data_size.push_back(0);

    // set the pointer to f each buffer to kDataBufferSize
    for (int i=0; i < num_parts_; i++)
      ptrs.push_back(kDataBufferSize);
  }

  void print_status()
  {
    std::cout << "Contx status " << std::endl;
    for (int i=0; i < num_parts_; i++)
    {
      index_ = i;
      std::cout << "\t" << index_ << "\t" << data[index_].size() << "\t"
                << data_size[index_] << "\t" << ptrs[index_] << std::endl;
    }
  }

  void clear()
  {
    data_size.clear();
    ptrs.clear();
    isempty.clear();
  }
};

using QUEUETUPLE = std::tuple<std::string, int, int>;

inline std::string ltrim(std::string& s)
{
  auto it = find_if_not(s.begin(), s.end(), [](int c) { return isspace(c); });
  s.erase(s.begin(), it);
  return s;
}

// remove the " (quotes) from the beginning and end of the string
// (TODO and the middle; hopefully nobody is trying to use escaped quotes).
std::string removeQuotes(std::string& s)
{
  s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c)
  {
    return c=='\"';
  }), s.end());
  return s;
}

std::vector<std::string> splitStringToFields(std::string const& str, char delim)
{
  std::stringstream splitter(str);
  std::vector<std::string> ret;
  for (std::string field; std::getline(splitter, field, delim); )
    ret.push_back(field);
  return ret;
}

class MitochondrialGeneSelector
{
public:
  MitochondrialGeneSelector(std::string const& mitochondrial_gene_names_filename)
  {
    if (mitochondrial_gene_names_filename.empty())
    {
      default_old_behavior_ = true;
      return;
    }

    std::ifstream input_file(mitochondrial_gene_names_filename);
    if (!input_file)
    {
      crash("ERROR failed to open the mitochondrial gene names file named: " +
            mitochondrial_gene_names_filename);
    }
    for (std::string line; std::getline(input_file, line);)
    {
      if (line.empty() || line[0] == '#') // skip comment lines
        continue;
      mito_genes_.insert(line);
    }
  }

  bool interestedInGeneName(std::string const& gene_name)
  {
    if (default_old_behavior_)
      return std::regex_search(gene_name, std::regex("^mt-", std::regex_constants::icase));
    else
      return mito_genes_.find(gene_name) != mito_genes_.end();
  }

private:
  bool default_old_behavior_ = false;
  std::unordered_set<std::string> mito_genes_;
};

// TODO function is named "get gene names", and there is something in there called
//      "gene name", but it instead returns a set of "gene id"s. correct?
//
// The file at gtf_filename should be unzipped.
std::unordered_set<std::string> get_mitochondrial_gene_names(
    std::string const& gtf_filename, std::string const& mitochondrial_gene_names_filename)
{
  std::unordered_set<std::string> mitochondrial_gene_ids;

  MitochondrialGeneSelector gene_selector(mitochondrial_gene_names_filename);

  std::ifstream input_file(gtf_filename);
  if (!input_file)
    crash("ERROR failed to open the GTF file named: " + gtf_filename);

  for (std::string line; std::getline(input_file, line);)
  {
    if (line.empty() || line[0] == '#') // skip comment lines
      continue;

    std::vector<std::string> tabbed_fields = splitStringToFields(line, '\t');
    assert(tabbed_fields.size() > 8);
    if (tabbed_fields[2] != "gene") // skip the line unless it is a gene
      continue;
    // split the semicolon-separated attributes field
    std::vector<std::string> attribs = splitStringToFields(tabbed_fields[8], ';');

    std::string gene_name;
    std::string gene_id;
    // now examine each of the attribute name-value pairs
    for (std::string attrib : attribs)
    {
      // each attribute is a space-separated key-value pair
      std::vector<std::string> key_and_val = splitStringToFields(ltrim(attrib), ' ');
      if (key_and_val.size() != 2)
        crash("Expected 2 fields, found " + std::to_string(key_and_val.size()) + " fields");

      // the second element in the pair is the value string
      std::string& key = key_and_val[0];
      std::string value = removeQuotes(key_and_val[1]);

      if (key == "gene_id")
        gene_id = value;
      if (key == "gene_name")
        gene_name = value;
    }
    if (gene_name.empty())
    {
      crash("Malformed GTF file detected. Record is of type gene but does not "
            "have a gene_name in line:\n" + line);
    }

    if (gene_selector.interestedInGeneName(gene_name))
      mitochondrial_gene_ids.insert(gene_id); // TODO what if gene_id is empty?
  }
  std::cout << "Number of mitochondrial genes found " << mitochondrial_gene_ids.size() << std::endl;
  return mitochondrial_gene_ids;
}


/*
 * @brief fills the buffer for the files
 *
 * @param contx is the context of the file
 * @return int number of alignments processed
*/
int fill_buffer(Context& contx, std::vector<std::string> const& partial_files)
{
  contx.data[contx.index_].clear();
  int k = 0;
  int filling_counter = 0;

  std::ifstream input_file(partial_files[contx.index_]);
  if (!input_file)
    crash("ERROR failed to open the file " + partial_files[contx.index_]);

  input_file.seekg(contx.file_offset[contx.index_]);

  // the order of the loop condition is iportant first make sure if you can accomodate then try to read,
  // otherwise it might create a read but never processed
  for (std::string line; k < kDataBufferSize && std::getline(input_file, line); k++)
  {
    contx.data[contx.index_].push_back(line);
    filling_counter++;
  }
  assert(contx.data[contx.index_].size() <= kDataBufferSize);

  contx.file_offset[contx.index_] = input_file.tellg();

  contx.data_size[contx.index_] = contx.data[contx.index_].size();

  if (contx.data_size[contx.index_] != 0)
  {
    contx.ptrs[contx.index_] = 0;
    contx.isempty[contx.index_] = false;
  }
  else
  {
    contx.ptrs[contx.index_] = kDataBufferSize;
    contx.isempty[contx.index_] = true;
  }

#ifdef DEBUG
  std::cout << "-->" << std::endl;
  for (int m = 0; m < contx.num_parts_; m++)
    std::cout << "\t" << m << " : " << contx.data_size[m] << " : " << contx.ptrs[m] << std::endl;
#endif

  return filling_counter;
}

// TODO if after other refactoring this ends up being the only regex use, then
//      probably would be worth switching away from regex here.
// From e.g. "A\tB\tC\tD\tE", extract "A\tB\tC"
std::string extractCompTag(std::string& s)
{
  const std::regex rgx("\t");
  const std::sregex_token_iterator end;
  std::sregex_token_iterator iter(s.begin(), s.end(), rgx, -1);
  std::stringstream comp_tag;
  for (auto k = 0; k < 3 && iter != end; ++iter, k++)
  {
    if (k > 0)
      comp_tag << "\t";
    comp_tag << *iter;
  }
  return comp_tag.str();
}

// returns number of alignments processed
int mergeSortedPartialFiles(INPUT_OPTIONS_TAGSORT const& options,
                            std::vector<std::string> const& partial_files)
{
  const std::string& sorted_output_file = options.sorted_output_file;
  const std::string& metric_type  = options.metric_type;
  const std::string& metric_output_file = options.metric_output_file;
  int filling_counter = 0;

  std::unordered_set<std::string> mitochondrial_genes;
  if (!options.gtf_file.empty())
  {
    mitochondrial_genes = get_mitochondrial_gene_names(
        options.gtf_file, options.mitochondrial_gene_names_filename);
  }

  // input the buffer size and partial files
  Context contx(partial_files.size());
  auto cmp = [](const QUEUETUPLE &a, const  QUEUETUPLE &b)
  {
    return std::get<0>(a) > std::get<0>(b);
  };
  std::priority_queue<QUEUETUPLE, std::vector<QUEUETUPLE>,  decltype(cmp) > heap(cmp);

  for (int i=0; i < contx.num_parts_; i++)
  {
    contx.index_ = i;
    filling_counter += fill_buffer(contx, partial_files);
  }

  // create the heap from the first batch loaded data
  contx.num_active_files = 0;
  for (int i=0; i< contx.num_parts_; i++)
  {
    contx.index_ = i;
    if (contx.ptrs[i] != kDataBufferSize)
    {
      heap.push(QUEUETUPLE(extractCompTag(contx.data[i][contx.ptrs[i]]), i, contx.ptrs[i]));
      contx.ptrs[i]++;
      contx.num_active_files += 1;
    }
  }

  //  now merge by pop an push
  std::ofstream fout;
  if (options.compute_metric) // TODO i think this is a mistake, and should actually be options.output_sorted_info
    fout.open(sorted_output_file);

  // pop and push from the heap
  int num_alignments = 0;
  int i, j;

  Metrics* metric_gatherer = nullptr;
  MetricType metric_type_enum = MetricType::Cell;
  if (metric_type.compare("cell")==0)
  {
    metric_gatherer = new CellMetrics;
    metric_type_enum = MetricType::Cell;
  }
  else if (metric_type.compare("gene")==0)
  {
    metric_gatherer = new GeneMetrics;
    metric_type_enum = MetricType::Gene;
  }
  else
    crash("Expected metric_type 'cell' or 'gene', got: " + metric_type);

  metric_gatherer->clear();

  std::ofstream fmetric_out;
  if (options.compute_metric)
  {
    fmetric_out.open(metric_output_file.c_str());
    fmetric_out << metric_gatherer->getHeader() << std::endl;
  }

  // TODO just write directly to fout... don't think the 'binary' is significant,
  //      but not 100% sure
  std::stringstream str(std::stringstream::out | std::stringstream::binary);
  std::string prev_comp_tag = "";
  while (!heap.empty())
  {
    // read the top
    QUEUETUPLE qtuple = heap.top();
    std::string curr_comp_tag = std::get<0>(qtuple);
    assert(prev_comp_tag.compare(curr_comp_tag) <= 0);

#ifdef DEBUG
    contx.print_status();
    if (prev_comp_tag.compare(curr_comp_tag) <= 0)
      std::cout << "Expected " << prev_comp_tag << "\n\t\t" << curr_comp_tag << std::endl;
    else
      crash("Anomaly " + prev_comp_tag + "\n\t\t" + curr_comp_tag);
#endif
    i = std::get<1>(qtuple);  //buffer no
    j = std::get<2>(qtuple);  //the pointer into the ith buffer array

    heap.pop();

    // start writing in chunks from the stream buffer
    if (num_alignments%kDataBufferSize==0)
    {
      if (options.output_sorted_info)
      {
        fout.write(str.str().c_str(), str.str().length());
        str.clear();
        str.str("");
      }
    }

    // load into stream buffer
    std::string field  = contx.data[i][j];
    if (options.output_sorted_info)
      str << field << std::endl;

    if (options.compute_metric)
      metric_gatherer->parse_line(field, fmetric_out, mitochondrial_genes, metric_type_enum);
    num_alignments += 1;

    // if ismpty is true means the file has been fully read
    if (!contx.isempty[i] && contx.ptrs[i] == contx.data_size[i])
    {
      contx.index_ = i;
      filling_counter += fill_buffer(contx, partial_files);
    }

    // make sure it is not empty
    if (contx.data_size[i] > 0)
    {
      heap.push(QUEUETUPLE(extractCompTag(contx.data[i][contx.ptrs[i]]), i, contx.ptrs[i]));
      contx.ptrs[i]++;
    }
    else     // one more file is fully read
      contx.num_active_files -= 1;

    if (num_alignments % 1000000 == 0)
      std::cout << "num alns read " << num_alignments << std::endl;

    prev_comp_tag = curr_comp_tag;
  }

  // process the final line
  metric_gatherer->finalize(mitochondrial_genes);
  metric_gatherer->output_metrics(fmetric_out);
  metric_gatherer->output_metrics_extra(fmetric_out);
  delete metric_gatherer;

  // close the metric file
  if (options.compute_metric)
    fmetric_out.close();

  // write out the remaining data
  if (options.output_sorted_info)
  {
    fout.write(str.str().c_str(), str.str().length());
    str.str("");
    str.clear();
  }

  // close output files as there is no more to write
  if (options.output_sorted_info)
    fout.close();

  std::cout << "Written "<< num_alignments << " alignments in total" << std::endl;
  contx.clear();
  return filling_counter;
}

void warnIfNo_mitochondrial_gene_names_filename(INPUT_OPTIONS_TAGSORT const& options)
{
  if (options.mitochondrial_gene_names_filename.empty())
  {
    std::string msg =
"*** WARNING! You did not specify --mitochondrial_gene_names_filename.\n"
"Therefore, we fell back to selecting only genes beginning with 'mt-' (case\n"
"insensitive). Please write a list of all gene names you're interested in into\n"
"a file, and pass the filename with --mitochondrial_gene_names_filename.";
    std::cout << msg << std::endl;
    std::cerr << msg << std::endl;
  }
}

/* Flag set by ‘--verbose’. */
int main(int argc, char** argv)
{
  INPUT_OPTIONS_TAGSORT options = readOptionsTagsort(argc, argv);
  warnIfNo_mitochondrial_gene_names_filename(options);

  std::cout << "bam input " << options.bam_input << std::endl;
  std::cout << "temp folder " << options.temp_folder << std::endl;
  std::cout << "sorted output file " <<  options.sorted_output_file << std::endl;
  std::cout << "metric output file " <<  options.metric_output_file << std::endl;
  std::cout << "temp folder " << options.alignments_per_thread << std::endl;
  std::cout << "tags:" << std::endl;

  for (auto const& [tag, tag_order_num] : options.tag_order)
    std::cout << "\t" << tag << "\t" << tag_order_num << std::endl;

  /* first create a list of sorted, and simplified sorted files */
  std::vector<std::string> partial_files = create_sorted_file_splits_htslib(options);

  /* now merge the sorted files to create one giant sorted file by using
    a head to compare the values based on the tags used  */
  std::cout << "Merging " <<  partial_files.size() << " sorted files!"<< std::endl;

  int filling_counter = mergeSortedPartialFiles(options, partial_files);

  // we no longer need the partial files
  for (unsigned int i=0; i < partial_files.size(); i++)
    if (remove(partial_files[i].c_str()) != 0)
      std::cerr << "Warning: error deleting file " << partial_files[i] << std::endl;

  partial_files.clear();
  std::cout << "Aligments " << filling_counter << " loaded to buffer " << std::endl;

  warnIfNo_mitochondrial_gene_names_filename(options);
  return 0;
}
