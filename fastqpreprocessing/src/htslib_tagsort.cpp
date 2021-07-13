/**
 *  @file   fastqprocess.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2020-08-27
 ***********************************************/
#include <htslib/sam.h>
#include "htslib_tagsort.h"

#define NUM_ALIGNMENTS_PER_CHUNK 1000000

extern "C" {
    bam_hdr_t * sam_hdr_read(samFile *); //read header
    htsFile *hts_open(const char *fn, const char *mode);
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

    samFile *fp_in = hts_open(options.bam_input.c_str(),"r"); //open bam file
    if (fp_in==0) std::cout << "Failed to open bam";
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment


    vector<TAGTUPLE>  tuple_records;
    char empty[] = {'\0'};

    std::cout << "Running htslib" << std::endl;
    char *barcode, *gene_id, *umi;
    char *location_tag;
    int nh_num;
    char *chr;

   // Keep reading records until ReadRecord returns false.
    long int i =0;
    int k = 0;
    long int num_alignments =0;
    uint8_t *p;
    while(sam_read1(fp_in, bamHdr,aln) > 0) {

       if ((p=bam_aux_get(aln, "CB"))==NULL) {
          barcode = empty;
       } else {
          barcode = bam_aux2Z(p);
       }
       if ((p=bam_aux_get(aln, "UB"))==NULL) {
          umi = empty;
       } else {
          umi = bam_aux2Z(p);
       }

       if ((p=bam_aux_get(aln, "GE"))==NULL) {
          gene_id = empty;
       } else {
          gene_id = bam_aux2Z(p);
       }

       if ((p=bam_aux_get(aln, "XF"))==NULL) {
          location_tag = empty;
       } else {
          location_tag = bam_aux2Z(p);
       }

       if ((p=bam_aux_get(aln, "NH"))==NULL) {
          nh_num = 0;
       } else {
          nh_num = bam_aux2i(p);
       }

       chr = bamHdr->target_name[aln->core.tid];

       if (aln->core.tid==-1){
          chr = empty;
       } else {
          chr = bamHdr->target_name[aln->core.tid];
       }

       uint8_t *q = bam_get_qual(aln);  // pointer to the qual data
       uint32_t len = aln->core.l_qseq; //length of qual seq.
       int threshold = 30;
       float avg_qual = 0;


       float qual_above_threshold = 0;
       for (k = 0; k < len; k++) {
          avg_qual += q[k];
          if (q[k]>threshold) qual_above_threshold += 1;
       }
       avg_qual /= len;
       qual_above_threshold /= len;

       if (i!=0 && i%num_align_per_file==0) {
           std::cout << "Batch number : " << partial_files.size() << std::endl;
           std::string split_file_path;
           split_file_path = write_out_partial_txt_file(tuple_records, tmp_folder);
           num_alignments += tuple_records.size();
           partial_files.push_back(split_file_path);

           tuple_records.clear();
//         if (partial_files.size() ==7) break;
       }

       tuple_records.push_back(
            std::make_tuple(std::string(barcode), std::string(umi), std::string(gene_id), 
                            std::string(chr), std::string(location_tag), 
                            avg_qual, qual_above_threshold, nh_num)
       );
       i = i + 1;
    }

    if (tuple_records.size()>0) {
        std::string split_file_path = write_out_partial_txt_file(tuple_records, tmp_folder);
        num_alignments += tuple_records.size();
        tuple_records.clear();
        partial_files.push_back(split_file_path);
     }

    std::cout << std::endl << "Read " << i << " records" << std::endl;
    std::cout << std::endl << "Read " << num_alignments << " records as batches" << std::endl;
    return partial_files;
}

}
