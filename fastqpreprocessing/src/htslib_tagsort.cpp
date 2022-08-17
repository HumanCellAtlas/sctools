/**
 *  @file   htslib_tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

constexpr int kThreshold = 30; // qual score threshold

#include "htslib_tagsort.h"

#include <tuple>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unordered_map>
#include <memory>
#include <semaphore.h>
#include <set>

sem_t g_tagsort_semaphore;
std::mutex g_mtx;
std::set<unsigned int> g_threads_to_join;

#define SEM_PRINT_VAL(X,Y)                               \
  ({                                                     \
     sem_getvalue(&X, &Y);                               \
     std::cout << "SEM_VALUE " << Y << std::endl;        \
  })                                                     \

#define SEM_INIT(X, n)  (sem_init(&X, 0, n) )

#define SEM_DESTROY(X)  (sem_destroy(&X) )

#define SEM_WAIT(X)                                      \
 ({                                                      \
     if (sem_wait(&X) == -1)                             \
         crashWithPerror("sem_wait: g_tagsort_semaphore");   \
 })

#define SEM_POST(X)                                      \
 ({                                                      \
     if (sem_post(&X) == -1)                             \
        crashWithPerror("sem_post: g_tagsort_semaphore");    \
 })


extern "C" {
  bam_hdr_t* sam_hdr_read(samFile*);   //read header
  htsFile* hts_open(const char* fn, const char* mode);
}


/*
  @brief get the int tag or -1

*/
inline int get_itag_or_default(bam1_t* aln, const char* tagname, int default_value)
{
  uint8_t* p;
  int tag_value = -1;
  if ((p = bam_aux_get(aln, tagname)) == NULL)
    tag_value = default_value;
  else
    tag_value = bam_aux2i(p);

  return tag_value;
}

/*
  @brief get the string tag or the default

*/
inline char* get_Ztag_or_default(bam1_t* aln, const char* tagname, char* default_value)
{
  uint8_t* p;
  char* tag_value = NULL;
  if ((p = bam_aux_get(aln, tagname)) == NULL)
    tag_value = default_value;
  else
  {
    tag_value = bam_aux2Z(p);
    if (strcmp(tag_value, "-") == 0)
      tag_value = default_value;
  }
  return tag_value;
}

const char kNone[] = "None";
const char kNoChr[] = "*";
void parseOneAlignment(std::vector<TAGTUPLE>* tuple_records, const bam1_t* aln,
                       InputOptionsTagsort const& options, const bam_hdr_t* bamHdr,
                       std::unordered_map<std::string, std::string*>& string_map)
{
  // extract the barcodes corrected and  corrected
  char* barcode = get_Ztag_or_default(aln, options.barcode_tag.c_str(), kNone);
  char* barcode_raw = get_Ztag_or_default(aln, "CR", "");

  // to be called perfect, the corrected and raw barcodes should match
  int perfect_cell_barcode = (strcmp(barcode, barcode_raw) == 0) ? 1 : 0;

  // barcode quality score
  char* barcode_qual = get_Ztag_or_default(aln, "CY", "");

  //average barcode across the query and the fraction of barcodes above threshold
  float sum_barcode_qual = 0;
  float num_bp_above_threshold = 0;
  size_t len = strlen(barcode_qual);
  for (unsigned int k = 0; k < len; k++)
  {
    // barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
    uint8_t qual_score = (((uint8_t)barcode_qual[k]) - 33);
    sum_barcode_qual += qual_score;
    if (qual_score > kThreshold)
      num_bp_above_threshold += 1;
  }
  int avg_cell_barcode_qual = sum_barcode_qual / (float)len; // TODO truncation intended?
  float cell_barcode_qual_above_threshold = (float)num_bp_above_threshold / (float)len;

  // corrected molecule barcodes (UMIs)
  char* umi = get_Ztag_or_default(aln, options.umi_tag.c_str(), kNone);
  // raw molecule barcodes
  char* umi_raw = get_Ztag_or_default(aln, "UR", "");

  // to be called perfect, the corrected and raw molecular barcodes should match
  int perfect_molecule_barcode = (strcmp(umi, umi_raw) == 0) ? 1 : 0;

  // qual score for molecular barcodes
  char* umi_qual = get_Ztag_or_default(aln, "UY", "");

  float sum_umi_qual = 0;
  float num_umi_above_threshold = 0;
  len = strlen(umi_qual);
  for (unsigned int k = 0; k < len; k++)
  {
    // molecular barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
    sum_umi_qual += ((uint8_t)umi_qual[k] -33);
    if (((uint8_t)umi_qual[k] - 33) > kThreshold)
      num_umi_above_threshold += 1;
  }
  float frac_umi_qual_above_threshold = (float)num_umi_above_threshold / (float)len;

  char* gene_id = get_Ztag_or_default(aln, options.gene_tag.c_str(), kNone);
  char* location_tag = get_Ztag_or_default(aln, "XF", "");

  int nh_num = get_itag_or_default(aln, "NH", -1);

  char* chr = (aln->core.tid == -1) ? kNoChr : bamHdr->target_name[aln->core.tid];

  uint32_t pos = aln->core.pos; // position.
  uint32_t isrev = bam_is_rev(aln) ? 1 : 0;   // is reverse stand
  uint32_t is_duplicate = ((aln->core.flag & BAM_FDUP) != 0) ? 1 : 0;

  // sequence quality score
  float avg_sequence_qual = 0, sum_qual = 0;
  float qual_above_threshold = 0;
  uint8_t* qual_seq = bam_get_qual(aln);  // pointer to the qual data
  len = aln->core.l_qseq; //length of qual seq.
  for (unsigned int k = 0; k < len; k++)
  {
    // the qual string are already in phred scores
    sum_qual += qual_seq[k];
    if (qual_seq[k] > kThreshold)
      qual_above_threshold += 1;
  }
  avg_sequence_qual = sum_qual / (float)len;
  qual_above_threshold = qual_above_threshold / (float)len;

  uint32_t* cigar = bam_get_cigar(aln);
  // see if it is spliced, i.e., N appears in the CIGAR string
  uint32_t spliced_read = 0;
  for (unsigned int k = 0; k < aln->core.n_cigar; k++)
  {
    uint32_t op = cigar[k] & BAM_CIGAR_MASK;
    if (op == 3 && (cigar[k] >> BAM_CIGAR_SHIFT) != 0)
    {
      spliced_read = 1;
      break;
    }
  }

  // the order of the three tags are define by the order of the supplied input arguments
  // tag.order [tag_name] -> order map
  std::string tags[3];
  tags[options.tag_order[options.barcode_tag]] = barcode;
  tags[options.tag_order[options.umi_tag]] = umi;
  tags[options.tag_order[options.gene_tag]] = gene_id;
  // TODO how long are these strings? the indirection might not be worth it
  if (string_map.find(barcode)==string_map.end())
    string_map[std::string(barcode)] = new std::string(barcode);
  if (string_map.find(umi)==string_map.end())
    string_map[std::string(umi)] = new std::string(umi);
  if (string_map.find(gene_id)==string_map.end())
    string_map[gene_id] = new std::string(gene_id);
  TRIPLET* triplet = new TRIPLET(string_map[tags[0]], string_map[tags[1]], string_map[tags[2]]);

  tuple_records->emplace_back(
      triplet, /* triplet of tags pointers */
      std::string(chr),  /* record[0] */
      std::string(location_tag), /* record[1] */
      pos,   /* record [2] */
      isrev, /* record[3] */
      avg_cell_barcode_qual,  /* record[4] */
      cell_barcode_qual_above_threshold, /* record[5] */
      avg_sequence_qual, /* record[6] */
      qual_above_threshold, /* record[7]  */
      nh_num, /* record[8] */
      perfect_molecule_barcode, /* record[9] */
      spliced_read, /* record[10] */
      is_duplicate, /* record[11] */
      perfect_cell_barcode,  /* record[12] */
      frac_umi_qual_above_threshold /* record[13] */);
}

inline bool sortbyfirst(std::pair<TRIPLET*, int> const& a,
                        std::pair<TRIPLET*, int> const& b)
{
  using std::get;
  if ((*get<0>(*a.first)).compare(*get<0>(*b.first)) != 0)
    return ((*get<0>(*a.first)).compare(*get<0>(*b.first)) < 0);
  if ((*get<1>(*a.first)).compare(*get<1>(*b.first)) != 0)
    return ((*get<1>(*a.first)).compare(*get<1>(*b.first)) < 0);
  return ((*get<2>(*a.first)).compare(*get<2>(*b.first)) < 0);
}

// Generates a random alphanumeric string (AZaz09) of a fixed length.
constexpr int kStringLen = 40;
std::string randomString()
{
  auto randchar = []() -> char
  {
    const char charset[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(kStringLen, 0);
  std::generate_n(str.begin(), kStringLen, randchar);
  return str;
}

using TRIPLET = std::tuple<std::string*, std::string*, std::string*>;

using TAGTUPLE = std::tuple<
    TRIPLET* /*  tuple<std::string *, std::string *, std::string *>*/,
    std::string /* reference */,
    std::string /* biotype */,
    int /* pos */,
    int /*rev strand   1 for yes, 0 otherwise*/,
    float /*avg barcode qual score */,
    float /* frac of barcode qual score >30 */,
    float /*avg qual seq */,
    float /*fract of >30 score qual seq*/,
    int /*NH*/,
    int /*perfect molecule barcode, 1 is yes, 0 otherwise*/,
    int /*spliced reads 1 yes, 0 otherwise*/,
    int /*is duplicate */,
    int /*perfect cell barcode 1 is yes, 0 otherwise*/,
    float /* fraction of umi qual score > 30 */
    >;

/**
 * @brief This function takes a vector of tuples of the tags, sorts them
 * in the dictionary order of the tags and then writes these in the same
 * order to a txt file
 *
 * @details
 * The function take the vector of tags tuples, writes the sorted tuples into
 * a file. The filename is generated randomly (enought to avoid collision with other files)
 * in the temp folder specified.
 *
 * @param tuple_records: vector<TAGTUPLE> &, reference to a vector of TAGTUPLES
 * @return a string for the random file name
*/
void write_out_partial_txt_file(std::vector<TAGTUPLE> const& tuple_records,
                                std::string const& tmp_folder,
                                std::vector<std::string>* partial_files)
{
  using std::get;

  std::string tempfile = tmp_folder + "/" + randomString() + ".txt";
  std::ofstream outfile(tempfile);

  // Sort by triplet, maintaining each triplet's pre-sorted index...
  std::vector<std::pair<TRIPLET*, int>> index_pairs;
  for (int i = 0; i < tuple_records.size(); i++)
    index_pairs.emplace_back(get<0>(tuple_records[i]), i);
  std::sort(index_pairs.begin(), index_pairs.end(), sortbyfirst);

  // ...then write the triplets in sorted order, linked up with the tuple_record
  // located at its pre-sorted index, i.e. the record for the triplet.
  for (auto [triplet_ptr, record_index] : index_pairs)
  {
    // TODO?
    // what if you ran out of disk space ???? NEED TO add logic
    outfile << *get<0>(*triplet_ptr) /*  first tag */ << "\t"
            << *get<1>(*triplet_ptr) /*  second tag */ << "\t"
            << *get<2>(*triplet_ptr) /*  third tag */ << "\t"
            << get<1>(tuple_records[record_index]) /* record[0] */  << "\t"
            << get<2>(tuple_records[record_index]) /* record[1] */  << "\t"
            << get<3>(tuple_records[record_index]) /* record[2] */  << "\t"
            << get<4>(tuple_records[record_index]) /* record[3] */  << "\t"
            << get<5>(tuple_records[record_index]) /* record[4] */  << "\t"
            << get<6>(tuple_records[record_index]) /* record[5] */  << "\t"
            << get<7>(tuple_records[record_index]) /* record[6] */  << "\t"
            << get<8>(tuple_records[record_index]) /* record[7] */  << "\t"
            << get<9>(tuple_records[record_index]) /* record[8] */  << "\t"
            << get<10>(tuple_records[record_index]) /* record[9] */ << "\t"
            << get<11>(tuple_records[record_index]) /* record[10] */ << "\t"
            << get<12>(tuple_records[record_index]) /* record[11] */ << "\t"
            << get<13>(tuple_records[record_index]) /* record[12] */ << "\t"
            << get<14>(tuple_records[record_index]) /* record[13] */ << "\n";
  }

  outfile.close(); // TODO is closing specifically here relevant to the mutex?
  g_mtx.lock();
  partial_files->push_back(tempfile);
  g_mtx.unlock();
}

unsigned int g_num_thread_deallocations = 0;

void process_alignments(InputOptionsTagsort const& options, bam1_t** alns,
                        bam_hdr_t* bamHdr, unsigned int buf_no, unsigned int n,
                        std::vector<std::string>* partial_files)
{
  std::vector<TAGTUPLE> tuple_records;
  // TODO this maps from a string to a heap-allocated copy of that string.
  std::unordered_map<std::string, std::string*> string_map;

  for (unsigned int i = 0; i < n; i++)
    parseOneAlignment(&tuple_records, alns[i], options, bamHdr, string_map);

  write_out_partial_txt_file(tuple_records, options.temp_folder, partial_files);

  // delete the triplet
  for (auto it=tuple_records.begin(); it != tuple_records.end(); ++it)
    delete std::get<0>(*it);

  //  free the memory for the strings
  for (auto it=string_map.begin(); it != string_map.end(); ++it)
    delete it->second;

  g_mtx.lock();
  g_threads_to_join.insert(buf_no);
  g_mtx.unlock();

  g_num_thread_deallocations += 1;
  SEM_POST(g_tagsort_semaphore);
}

/**
 * @brief From the input bam create a list of txt files with the records (lines)
 * sorted according to the * tags
 *
 * @details
 * The input bam file is read chunk by chunk, sorted by the tags and the written
 * out as a text file in the sorted manner.
 *
 * @param options: InputOptionsTagsort the inputs to the program
 * @return a vector containing the file paths of the partial files
*/
std::vector<std::string> create_sorted_file_splits_htslib(InputOptionsTagsort& options)
{
  // size of individual chunks to sort in memory as an approx 20 mil alignments makes 1 GB bam
  //open bam file
  samFile* fp_in=0;
  if ((fp_in = hts_open(options.bam_input.c_str(),"r"))==0)
    crash(options.bam_input + ": cannot open file.");

  bam_hdr_t* bamHdr = sam_hdr_read(fp_in); //read header
  bam1_t*** aln_arr;

  // allocated memory for the alignments for the various threads
  if ((aln_arr = (bam1_t***)malloc(sizeof(bam1_t**)*options.nthreads))==NULL)
    crash("ERROR Failed to allocate memory");

  for (unsigned int i = 0; i < options.nthreads; i++)
  {
    if ((aln_arr[i] = (bam1_t**)malloc(sizeof(bam1_t*)*options.alignments_per_thread))==NULL)
      crash("ERROR Failed to allocate memory for alignments");
    for (unsigned int k = 0; k < options.alignments_per_thread; k++)
      aln_arr[i][k]= bam_init1(); //initialize an alignment
  }

  std::cout << "Running htslib" << std::endl;

  unsigned int batch = 0;
  long int num_alignments = 0;
  std::string tags[3];

  std::thread thread_id[kMaxTagsortThreads];

  std::set<unsigned int> busy_buffers, idle_buffers;
  for (unsigned int j = 0; j < options.nthreads; j++)
    idle_buffers.insert(j);

  SEM_INIT(g_tagsort_semaphore, options.nthreads);
  SEM_WAIT(g_tagsort_semaphore);

  g_mtx.lock();
  unsigned int buf_no = *idle_buffers.begin();
  idle_buffers.erase(buf_no);
  busy_buffers.insert(buf_no);
  g_mtx.unlock();

  std::vector<std::string> partial_files;

  unsigned int ind_in_buf = 0;
  while (sam_read1(fp_in, bamHdr, aln_arr[buf_no][ind_in_buf]) > 0)
  {
    ind_in_buf++;
    num_alignments++;

    if (ind_in_buf == options.alignments_per_thread)
    {
      std::cout << "Batch number : " << batch << std::endl;
      batch++;

      // thread begins with the new buf no
      thread_id[buf_no] = std::thread(process_alignments, std::ref(options),
                                      aln_arr[buf_no], bamHdr, buf_no, ind_in_buf,
                                      &partial_files);

      // move forward only if there is room for new buffer
      SEM_WAIT(g_tagsort_semaphore);

      // join any completed thread
      g_mtx.lock();
      for (unsigned int thread_ind : g_threads_to_join)
      {
        thread_id[thread_ind].join();
        busy_buffers.erase(thread_ind);
        idle_buffers.insert(thread_ind);
      }
      // clear the threads
      g_threads_to_join.clear();
      g_mtx.unlock();

      // now find an idle buffer and make it busy so that we can load new data
      g_mtx.lock();
      // this is the new buffer to fill
      buf_no = *idle_buffers.begin();
      idle_buffers.erase(buf_no);
      busy_buffers.insert(buf_no);
      //std::cout << buf_no << " busy" << std::endl;
      g_mtx.unlock();
      ind_in_buf = 0;
    }
  }

  thread_id[buf_no] = std::thread(process_alignments, std::ref(options),
                                  aln_arr[buf_no], bamHdr, buf_no, ind_in_buf,
                                  &partial_files);

  // make sure you consume all the counts, the remaining semaphore is the value
  int sem_value;
  sem_getvalue(&g_tagsort_semaphore, &sem_value);
  for (unsigned int k=0; k< options.nthreads; k++)
    SEM_WAIT(g_tagsort_semaphore);
  SEM_DESTROY(g_tagsort_semaphore);

  // join any completed thread that are ready to be joined
  g_mtx.lock();
  for (unsigned int thread_ind : g_threads_to_join)
    thread_id[thread_ind].join();

  g_threads_to_join.clear();
  g_mtx.unlock();

  // release the memory and clear data-structures
  for (unsigned int i = 0; i < options.nthreads; i++)
  {
    std::cout << "releasing memory for thread " << i << std::endl;
    for (unsigned int k = 0; k < options.alignments_per_thread; k++)
      bam_destroy1(aln_arr[i][k]);
    free(aln_arr[i]);
  }
  free(aln_arr);

  sam_hdr_destroy(bamHdr);
  hts_close(fp_in);

  std::cout << "Deallocate threads " << g_num_thread_deallocations << std::endl;
  std::cout << std::endl << "Read " << i << " records" << std::endl;
  std::cout << "Read " << num_alignments << " records as batches" << std::endl;

  return partial_files;
}  // function

