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
  int filling_counter = 0;

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
    fout.open(options.sorted_output_file);

  // pop and push from the heap
  int num_alignments = 0;
  int i, j;

  std::unique_ptr<MetricGatherer> metric_gatherer;
  if (options.metric_type == MetricType::Cell)
  {
    metric_gatherer = std::make_unique<CellMetricGatherer>(
        options.gtf_file, options.mitochondrial_gene_names_filename);
  }
  else if (options.metric_type == MetricType::Gene)
    metric_gatherer = std::make_unique<GeneMetricGatherer>();
  else if (options.metric_type == MetricType::Umi)
    metric_gatherer = std::make_unique<UmiMetricGatherer>();
  else
    crash("new MetricType enum value is not yet handled by MetricGatherer!");
  metric_gatherer->clear();

  std::ofstream fmetric_out;
  if (options.compute_metric)
  {
    fmetric_out.open(options.metric_output_file.c_str());
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
    i = std::get<1>(qtuple);  //buffer no
    j = std::get<2>(qtuple);  //the pointer into the ith buffer array

#ifdef DEBUG
    assert(prev_comp_tag.compare(curr_comp_tag) <= 0);
    contx.print_status();
    if (prev_comp_tag.compare(curr_comp_tag) <= 0)
      std::cout << "Expected " << prev_comp_tag << "\n\t\t" << curr_comp_tag << std::endl;
    else
      crash("Anomaly " + prev_comp_tag + "\n\t\t" + curr_comp_tag);
#endif

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
    std::string field = contx.data[i][j];
    if (options.output_sorted_info)
      str << field << std::endl;

    if (options.compute_metric)
      metric_gatherer->ingestLine(field, fmetric_out);
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
  metric_gatherer->output_metrics(fmetric_out);
  metric_gatherer->output_metrics_extra(fmetric_out);

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
  std::cout << "temp folder " << options.alignments_per_batch << std::endl;
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
