/**
 *  @file   htslib_tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include "globals.h"
#include "htslib_tagsort.h"

unsigned int num_thread_deallocations = 0;
extern sem_t semaphore;
extern std::mutex mtx;
extern std::set<unsigned int> busy_buffer, idle_buffer, threads_to_join;

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
         crashWithPerror("sem_wait: semaphore");                          \
 })

#define SEM_POST(X)                                      \
 ({                                                      \
     if (sem_post(&X) == -1)                             \
        crashWithPerror("sem_post: semaphore");                           \
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
  if ((p=bam_aux_get(aln, tagname))==NULL)
  {
    tag_value = default_value;
  }
  else
  {
    tag_value = bam_aux2i(p);
  }
  return  tag_value;
}

/*
  @brief get the string tag or the default

*/
inline char* get_Ztag_or_default(bam1_t* aln, const char* tagname, char* default_value)
{
  uint8_t* p;
  char* tag_value = NULL;
  if ((p=bam_aux_get(aln, tagname))==NULL)
  {
    tag_value = default_value;
  }
  else
  {
    tag_value = bam_aux2Z(p);
    if (strcmp(tag_value, "-")==0)
    {
      tag_value = default_value;
    }
  }
  return  tag_value;
}

void process_alignments(InputOptionsTagsort& options, bam1_t** aln, bam_hdr_t* bamHdr, unsigned int buf_no, unsigned int n)
{
  int threshold = THRESHOLD; //qual score threshold
  std::vector<TAGTUPLE>  tuple_records;
  std::unordered_map<std::string, std::string*>  string_map;
  std::string tags[3];

  char empty[] = {'\0'};
  char none[] = {'N', 'o', 'n', 'e', '\0'};
  char nochr[] = {'*', '\0'};

  char* barcode, *gene_id, *umi, *umi_raw, *umi_qual,  *barcode_qual, *barcode_raw;
  char* location_tag;
  int nh_num;
  char* chr;

  uint32_t len = 0; //length of qual seq.

  for (unsigned int i = 0; i < n; i++)
  {
    // extract the barcodes corrected and  corrected
    barcode = (char*)get_Ztag_or_default(aln[i], options.barcode_tag.c_str(), none);
    barcode_raw = (char*)get_Ztag_or_default(aln[i], "CR", empty);

    // to be called perfect the corrected and raw barcodes should match
    int perfect_cell_barcode = 1;
    unsigned int barcode_len = strlen(barcode);

    if (strcmp(barcode, "None")==0)
    {
      // empty barcodes are not perfect
      perfect_cell_barcode = 0;
    }
    else
    {
      // perfect barcodes have the same corrected and raw barcodes
      for (unsigned int k=0; k < barcode_len; k++)
        if (barcode[k] != barcode_raw[k])
        {
          perfect_cell_barcode = 0;
          break;
        }
    }

    // barcode quality score
    barcode_qual = (char*)get_Ztag_or_default(aln[i], "CY", empty);

    //average barcode across the query and the fraction of barcodes above threshold
    float sum_barcode_qual = 0;
    float num_bp_above_threshold = 0;
    len = strlen(barcode_qual);

    for (unsigned int k = 0; k < len; k++)
    {
      // barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
      sum_barcode_qual += ((uint8_t)barcode_qual[k] -33);
      if (((uint8_t)barcode_qual[k] - 33)  > threshold) num_bp_above_threshold += 1;
    }
    int avg_cell_barcode_qual = sum_barcode_qual/(float)len;
    float cell_barcode_qual_above_threshold = (float)num_bp_above_threshold/(float)len;


    // corrected molecule barcodes (UMIs)
    umi = (char*)get_Ztag_or_default(aln[i], options.umi_tag.c_str(), none);
    // raw molecule barcodes
    umi_raw = (char*)get_Ztag_or_default(aln[i], "UR", empty);

    // to be called perfect,  the corrected and raw molecular barcodes should match
    unsigned int umi_len = strlen(umi);
    int perfect_molecule_barcode = 1;
    if (strcmp(umi, "None")==0)    // not equal length case
    {
      perfect_molecule_barcode = 0;
    }
    else
    {
      // if equal then compare char by char
      for (unsigned int k=0; k < umi_len; k++)
        if (umi[k] != umi_raw[k])
        {
          perfect_molecule_barcode = 0;
          break;
        }
    }
    // qual score for molecular barcodes
    umi_qual = (char*)get_Ztag_or_default(aln[i], "UY", empty);

    float sum_umi_qual = 0;
    float num_umi_above_threshold = 0;
    len = strlen(umi_qual);
    for (unsigned int k = 0; k < len; k++)
    {
      // molecular barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
      sum_umi_qual += ((uint8_t)umi_qual[k] -33);
      if (((uint8_t)umi_qual[k] - 33)  > threshold) num_umi_above_threshold += 1;
    }
    float frac_umi_qual_above_threshold = (float)num_umi_above_threshold/(float)len;

    gene_id = (char*)get_Ztag_or_default(aln[i], options.gene_tag.c_str(), none);
    location_tag = (char*)get_Ztag_or_default(aln[i], "XF", empty);

    nh_num = get_itag_or_default(aln[i], "NH", -1);

    chr = bamHdr->target_name[aln[i]->core.tid];
    if (aln[i]->core.tid==-1)
    {
      chr = nochr;
    }
    else
    {
      chr = bamHdr->target_name[aln[i]->core.tid];
    }

    uint32_t pos = aln[i]->core.pos; // position.
    uint32_t isrev = bam_is_rev(aln[i]) ? 1 : 0;   // is reverse stand

    uint32_t is_duplicate = 0;  // is duplicate
    if ((aln[i]->core.flag & BAM_FDUP) !=0)
    {
      is_duplicate = 1;
    }

    // sequence quality score
    uint8_t* qual_seq = bam_get_qual(aln[i]);  // pointer to the qual data
    len = aln[i]->core.l_qseq; //length of qual seq.

    float avg_sequence_qual = 0, sum_qual = 0;
    float qual_above_threshold = 0;

    for (unsigned int k = 0; k < len; k++)
    {
      // the qual string are already in phred scores
      sum_qual += (qual_seq[k]);
      if ((qual_seq[k]) > threshold) qual_above_threshold += 1;
    }
    avg_sequence_qual = sum_qual/(float)len;
    qual_above_threshold = qual_above_threshold/(float)len;

    uint32_t* cigar = bam_get_cigar(aln[i]);
    // see if it is spliced, i.e., N appears in the CIGAR string
    uint32_t spliced_read = 0;
    for (unsigned int k = 0; k < aln[i]->core.n_cigar; ++k)
    {
      uint32_t op = cigar[k] & BAM_CIGAR_MASK;
      if (op==3 && (cigar[k] >> BAM_CIGAR_SHIFT)!=0)
      {
        spliced_read = 1;
        break;
      }
    }

    // the order of the three tags are define by the order of the supplied input arguments
    // tag.order [tag_name] -> order map
    tags[options.tag_order[options.barcode_tag]] = barcode;
    tags[options.tag_order[options.umi_tag]] = umi;
    tags[options.tag_order[options.gene_tag]] = gene_id;

    if (string_map.find(barcode)==string_map.end())
    {
      string_map[std::string(barcode)] = new std::string(barcode);
    }

    if (string_map.find(umi)==string_map.end())
    {
      string_map[std::string(umi)] = new std::string(umi);
    }

    if (string_map.find(gene_id)==string_map.end())
    {
      string_map[gene_id] = new std::string(gene_id);
    }

    TRIPLET* triplet = new TRIPLET(string_map[tags[0]],  string_map[tags[1]], string_map[tags[2]]);

    tuple_records.push_back(
      std::make_tuple(
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
        frac_umi_qual_above_threshold /* record[13] */
      )
    );
  }//for loop

  write_out_partial_txt_file(tuple_records, options.temp_folder);

  // delete the triplet
  for (auto it=tuple_records.begin(); it != tuple_records.end(); it++)
    delete std::get<0>(*it);

  //  free the memory for the strings
  for (auto it=string_map.begin(); it != string_map.end(); it++)
    delete it->second;

  mtx.lock();
  threads_to_join.insert(buf_no);
  mtx.unlock();

  num_thread_deallocations += 1;
  SEM_POST(semaphore);
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
void create_sorted_file_splits_htslib(InputOptionsTagsort& options)
{
  std::string input_bam = options.bam_input;
  std::string tmp_folder = options.temp_folder;

  // size of individual chunks to sort in memory as an approx 20 mil alignments makes 1 GB bam
  //open bam file
  samFile* fp_in=0;
  if ((fp_in = hts_open(options.bam_input.c_str(),"r"))==0)
  {
    std::cerr << "ERROR Cannot open the file " << options.bam_input << std::endl;
    exit(1);
  }

  bam_hdr_t* bamHdr = sam_hdr_read(fp_in); //read header
  bam1_t*** aln_arr;

  // allocated memory for the alignments for the various threads
  if ((aln_arr = (bam1_t***)malloc(sizeof(bam1_t**)*options.nthreads))==NULL)
  {
    std::cerr << "ERROR Failed to allocate memory " <<   std::endl;
    exit(1);
  }
  for (unsigned int i = 0; i < options.nthreads; i++)
  {
    if ((aln_arr[i] = (bam1_t**)malloc(sizeof(bam1_t*)*options.alignments_per_thread))==NULL)
    {
      std::cerr << "ERROR Failed to allocate memory for alignments " <<   std::endl;
    }
    for (unsigned int k = 0; k < options.alignments_per_thread; k++)
    {
      aln_arr[i][k]= bam_init1(); //initialize an alignment
    }
  }

  std::cout << "Running htslib" << std::endl;
  // Keep reading records until ReadRecord returns false.
  long int i =0;
  unsigned int batch = 0;
  long int num_alignments =0;
  std::string tags[3];

  std::thread thread_id[MAX_THREADS];

  unsigned int buf_i=0;

  for (unsigned int j = 0; j < options.nthreads; j++)
  {
    idle_buffers.insert(j);
  }

  std::set<unsigned int>::iterator buf_it;
  unsigned int buf_no;

  SEM_INIT(semaphore, options.nthreads);

  int sem_value;
  SEM_WAIT(semaphore);

  mtx.lock();
  buf_no = *idle_buffers.begin();
  idle_buffers.erase(buf_no);
  busy_buffers.insert(buf_no);
  mtx.unlock();

  while (sam_read1(fp_in, bamHdr,aln_arr[buf_no][buf_i]) > 0)
  {
    buf_i++;
    num_alignments++;

    if (buf_i!=0 && buf_i==options.alignments_per_thread)
    {
      std::cout << "Batch number : " << batch << std::endl;
      batch++;

      // thread begins with the new buf no
      thread_id[buf_no] = std::thread(process_alignments, std::ref(options), aln_arr[buf_no], bamHdr, buf_no, buf_i);

      // move forward only if there is room for new buffer
      SEM_WAIT(semaphore);

      // join any completed thread
      mtx.lock();
      for (buf_it=threads_to_join.begin(); buf_it!=threads_to_join.end(); buf_it++)
      {
        thread_id[*buf_it].join();
        busy_buffers.erase(*buf_it);
        idle_buffers.insert(*buf_it);
      }
      // clear the threads
      threads_to_join.clear();
      mtx.unlock();

      // now find an idle buffer and make it busy so that we can load new data
      mtx.lock();
      // this is the new buffer to fill
      buf_no = *idle_buffers.begin();
      idle_buffers.erase(buf_no);
      busy_buffers.insert(buf_no);
      //std::cout << buf_no << " busy" << std::endl;
      mtx.unlock();
      buf_i = 0;
    }
    i = i + 1;
  }  // while loop

  thread_id[buf_no] = std::thread(process_alignments, std::ref(options), aln_arr[buf_no], bamHdr, buf_no, buf_i);

  // make sure you consume all the counts, the remaining semaphore is the value
  sem_getvalue(&semaphore, &sem_value);
  for (unsigned int k=0; k< options.nthreads; k++)
  {
    SEM_WAIT(semaphore);
  }
  SEM_DESTROY(semaphore);

  // join any completed thread that are ready to be joined
  mtx.lock();
  for (buf_it=threads_to_join.begin(); buf_it!=threads_to_join.end(); buf_it++)
  {
    thread_id[*buf_it].join();
  }
  threads_to_join.clear();
  mtx.unlock();

  // release the memory and clear data-structures
  for (unsigned int i = 0; i < options.nthreads; i++)
  {
    std::cout << "releasing memory for thread " << i << std::endl;
    for (unsigned int k = 0; k < options.alignments_per_thread; k++)
    {
      bam_destroy1(aln_arr[i][k]);
    }
    free(aln_arr[i]);
  }
  free(aln_arr);

  sam_hdr_destroy(bamHdr);
  hts_close(fp_in);

  std::cout << "Deallocate threads " << num_thread_deallocations << std::endl;
  std::cout << std::endl << "Read " << i << " records" << std::endl;
  std::cout << "Read " << num_alignments << " records as batches" << std::endl;
}  // function

