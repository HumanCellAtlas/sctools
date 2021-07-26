/**
 *  @file   htslib_tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/
#include "htslib_tagsort.h"

#define NUM_ALIGNMENTS_PER_CHUNK 1000000
#define THRESHOLD 30.0

extern sem_t semaphore;

#define SEM_INIT(X)  (sem_init(&X, 0, 1) )
#define SEM_WAIT(X)                                      \
 ({                                                      \
     if (sem_wait(&X) == -1)                             \
         error("sem_wait: semaphore");                   \
 })


extern "C" {
    bam_hdr_t * sam_hdr_read(samFile *); //read header
    htsFile *hts_open(const char *fn, const char *mode);
}

/*
  @brief get the int tag or -1
  
*/
inline int get_itag_or_default(bam1_t *aln, const char *tagname, int default_value) {
    uint8_t *p;
    int tag_value = -1;
    if ((p=bam_aux_get(aln, tagname))==NULL) {
         tag_value = default_value;
    } else {
        tag_value = bam_aux2i(p);
    }
    return  tag_value;
}

/*
  @brief get the string tag or the default
  
*/
inline char *get_Ztag_or_default(bam1_t *aln, const char *tagname, char *default_value) {
    uint8_t *p;
    char *tag_value = NULL;
    if ((p=bam_aux_get(aln, tagname))==NULL) {
         tag_value = default_value;
    } else {
        tag_value = bam_aux2Z(p);
    }
    return  tag_value;
}

using namespace std;

/**
 * @brief From the input bam create a list of txt files with the records (lines) 
 * sorted according to the * tags 
 *
 * @details
 * The input bam file is read chunk by chunk, sorted by the tags and the written 
 * out as a text file in the sorted manner. 
 *
 * @param options: INPUT_OPTIONS_TAGSORT the inputs to the program
 * @return a vector containing the file paths of the partial files
*/
namespace htslib {
std::vector<std::string> create_sorted_file_splits_htslib(INPUT_OPTIONS_TAGSORT &options) {

    string input_bam = options.bam_input;
    string tmp_folder = options.temp_folder;
    
    // size of individual chunks to sort in memory as an approx 20 mil alignments makes 1 GB bam
    int num_align_per_file =  static_cast<int>(options.inmemory_chunk_size * NUM_ALIGNMENTS_PER_CHUNK);
    std::vector<string> partial_files;

    //open bam file

    samFile *fp_in=0;
    if ((fp_in = hts_open(options.bam_input.c_str(),"r"))==0) {
        std::cerr << "ERROR Cannot open the file " << options.bam_input << std::endl;
        exit(1);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment


    char empty[] = {'\0'};
    char none[] = {'N', 'o', 'n', 'e', '\0'};
    char nochr[] = {'*', '\0'};

    std::cout << "Running htslib" << std::endl;
    char *barcode, *gene_id, *umi, *umi_raw, *umi_qual,  *barcode_qual, *barcode_raw;
    char *location_tag;
    int nh_num;
    char *chr;

   // Keep reading records until ReadRecord returns false.
    long int i =0;
    unsigned int k = 0;
    unsigned int batch = 0;
    long int num_alignments =0;
    uint32_t len = 0; //length of qual seq.
    int threshold = THRESHOLD; //qual score threshold

    vector<TAGTUPLE>  tuple_records[2];
    std::unordered_map<std::string, std::string *>  string_map[2];

    std::string tags[3]; 

    SEM_INIT(semaphore);
    std::vector<std::thread> thread_ids;

    while(sam_read1(fp_in, bamHdr,aln) > 0) {
/*
       if (i > 25000) { 
           std::cout << "Remove me from file " << __FILE__ << " line no " << __LINE__ << std::endl;
           break;
       }
*/
       // extract the barcodes corrected and  corrected
       barcode = (char *)get_Ztag_or_default(aln, options.barcode_tag.c_str(), none);
       barcode_raw = (char *)get_Ztag_or_default(aln, "CR", empty);

       // to be called perfect the corrected and raw barcodes should match
       int perfect_cell_barcode = 1;
       unsigned int barcode_len = strlen(barcode);

       if (strcmp(barcode, "None")==0)  {
          // empty barcodes are not perfect
           perfect_cell_barcode = 0;     
       } else {
          // perfect barcodes have the same corrected and raw barcodes
          for (k=0; k < barcode_len; k++) 
            if (barcode[k] != barcode_raw[k]) {
               perfect_cell_barcode = 0;     
               break;
            }
       }

       // barcode quality score
       barcode_qual = (char *)get_Ztag_or_default(aln, "CY", empty);

       //average barcode across the query and the fraction of barcodes above threshold
       float sum_barcode_qual = 0;
       float num_bp_above_threshold = 0;
       len = strlen(barcode_qual);
       for (k = 0; k < len; k++) {
          // barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
          sum_barcode_qual += ((uint8_t)barcode_qual[k] -33);
          if (((uint8_t)barcode_qual[k] - 33)  > threshold) num_bp_above_threshold += 1;
       }
       int avg_cell_barcode_qual = sum_barcode_qual/(float)len;
       float cell_barcode_qual_above_threshold =  (float)num_bp_above_threshold/(float)len;


       // corrected molecule barcodes (UMIs)
       umi = (char *)get_Ztag_or_default(aln, options.umi_tag.c_str(), none);
       // raw molecule barcodes
       umi_raw = (char *)get_Ztag_or_default(aln, "UR", empty);

       // to be called perfect,  the corrected and raw molecular barcodes should match 
       unsigned int umi_len = strlen(umi);
       int perfect_molecule_barcode = 1;
       if (strcmp(umi, "None")==0) {  // not equal length case
          perfect_molecule_barcode = 0;     
       } else {   
          // if equal then compare char by char
          for (k=0; k < umi_len; k++) 
            if (umi[k] != umi_raw[k]) {
              perfect_molecule_barcode = 0;     
              break;
           }
       }

       // qual score for molecular barcodes
       umi_qual = (char *)get_Ztag_or_default(aln, "UY", empty);

       float sum_umi_qual = 0;
       float num_umi_above_threshold = 0;
       len = strlen(umi_qual);
       for (k = 0; k < len; k++) {
          // molecular barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
          sum_umi_qual += ((uint8_t)umi_qual[k] -33);
          if (((uint8_t)umi_qual[k] - 33)  > threshold) num_umi_above_threshold += 1;
       }
       float frac_umi_qual_above_threshold =  (float)num_umi_above_threshold/(float)len;

       gene_id = (char *)get_Ztag_or_default(aln, options.gene_tag.c_str(), none);
       location_tag = (char *)get_Ztag_or_default(aln, "XF", empty);

       nh_num = get_itag_or_default(aln, "NH", -1);

       chr = bamHdr->target_name[aln->core.tid];
       if (aln->core.tid==-1){
          chr = nochr;
       } else {
          chr = bamHdr->target_name[aln->core.tid];
       }

       uint32_t pos = aln->core.pos; // position.
       uint32_t isrev = bam_is_rev(aln) ? 1 : 0;   // is reverse stand 
       
       uint32_t is_duplicate = 0;  // is duplicate
       if ((aln->core.flag & BAM_FDUP) !=0) {
          is_duplicate = 1;
       }

       // sequence quality score
       uint8_t *qual_seq = bam_get_qual(aln);  // pointer to the qual data
       len = aln->core.l_qseq; //length of qual seq.

       float avg_sequence_qual = 0, sum_qual = 0;
       float qual_above_threshold = 0;

       for (k = 0; k < len; k++) {
          // the qual string are already in phred scores
          sum_qual += (qual_seq[k]);
          if ((qual_seq[k]) > threshold) qual_above_threshold += 1;
       }
       avg_sequence_qual = sum_qual/(float)len;
       qual_above_threshold = qual_above_threshold/(float)len;


       uint32_t *cigar = bam_get_cigar(aln);
       // see if it is spliced, i.e., N appears in the CIGAR string
       uint32_t spliced_read = 0;
       for (k = 0; k < aln->core.n_cigar; ++k) {
           uint32_t op = cigar[k] & BAM_CIGAR_MASK;
           if (op==3 && (cigar[k] >> BAM_CIGAR_SHIFT)!=0) {
              spliced_read = 1;
              break;
           } 
       }
   
       if (i!=0 && i%num_align_per_file==0) {
           std::cout << "Batch number : " << batch << std::endl;
           std::string split_file_path;
           
           SEM_WAIT(semaphore);
           std::thread thread_id =  std::thread(write_out_partial_txt_file, tuple_records[batch%2],  tmp_folder, std::ref(partial_files));

           thread_ids.push_back(std::move(thread_id));

           num_alignments += tuple_records[batch%2].size();
           batch++;

           for(auto it=tuple_records[batch%2].begin(); it!=tuple_records[batch%2].end(); it++) { 
              delete get<0>(*it);
           }
           tuple_records[batch%2].clear(); 

           for(auto it=string_map[batch%2].begin(); it!=string_map[batch%2].end(); it++) { 
              delete it->second; 
           }
           string_map[batch%2].clear();
       }

       // the order of the three tags are define by the order of the supplied input arguments 
       // tag.order [tag_name] -> order map  
       tags[options.tag_order[options.barcode_tag]] = barcode;
       tags[options.tag_order[options.umi_tag]] = umi;
       tags[options.tag_order[options.gene_tag]] = gene_id;
   
       if (string_map[batch%2].find(barcode)==string_map[batch%2].end()) {
           string_map[batch%2][barcode] = new std::string(barcode);
       }
       if (string_map[batch%2].find(umi)==string_map[batch%2].end()) {
           string_map[batch%2][umi] = new std::string(umi);
       }
       if (string_map[batch%2].find(gene_id)==string_map[batch%2].end()) {
           string_map[batch%2][gene_id] = new std::string(gene_id);
       }
        
       TRIPLET* triplet = new TRIPLET(string_map[batch%2][tags[0]],  string_map[batch%2][tags[1]], string_map[batch%2][tags[2]]);

       tuple_records[batch%2].push_back(
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
       i = i + 1;
    }

    // release the memory and clear data-structures
    bam_destroy1(aln);
    sam_hdr_destroy(bamHdr);
    hts_close(fp_in);

    if (tuple_records[batch%2].size()>0) {
        SEM_WAIT(semaphore);

        write_out_partial_txt_file(tuple_records[batch%2], tmp_folder, partial_files);
        num_alignments += tuple_records[batch%2].size();

        for(auto it=tuple_records[batch%2].begin(); it!=tuple_records[batch%2].end(); it++) { 
            delete get<0>(*it);
        }
        tuple_records[batch%2].clear(); 

        for(auto it=string_map[batch%2].begin(); it!=string_map[batch%2].end(); it++) { 
            delete it->second; 
        }
        string_map[batch%2].clear();
        //partial_files.push_back(split_file_path);
    }


    for(auto it=thread_ids.begin(); it!=thread_ids.end(); it++) { 
       it->join();
    }
    thread_ids.erase(thread_ids.begin(), thread_ids.end()); 
    
    std::cout << std::endl << "Read " << i << " records" << std::endl;
    std::cout << std::endl << "Read " << num_alignments << " records as batches" << std::endl;

    // one of them has data, which is from the last thread inside the while loop
    // but we can clean both anyway
    for (k = 0; k < 2; k++)  {
      for(auto it=tuple_records[k%2].begin(); it!=tuple_records[k%2].end(); it++) { 
         delete get<0>(*it);
      }
      tuple_records[k%2].clear(); 
    }

    for (k = 0; k < 2; k++)  {
      for(auto it=string_map[k%2].begin(); it!=string_map[k%2].end(); it++) { 
        delete it->second; 
      }
      string_map[k%2].clear();
    }

    return partial_files;
}  //while loop

} //namespace 
