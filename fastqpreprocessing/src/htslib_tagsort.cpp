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

    samFile *fp_in = hts_open(options.bam_input.c_str(),"r"); //open bam file
    if (fp_in==0) std::cout << "Failed to open bam";
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
    bam1_t *aln = bam_init1(); //initialize an alignment


    vector<TAGTUPLE>  tuple_records;
    char empty[] = {'\0'};
    char none[] = {'N', 'o', 'n', 'e', '\0'};
    char nochr[] = {'*', '\0'};

    std::cout << "Running htslib" << std::endl;
    char *barcode, *gene_id, *umi, *umi_raw, *umi_qual,  *barcode_qual, *barcode_raw;
    char *qual_seq;
    char *location_tag;
    int nh_num;
    char *chr;

   // Keep reading records until ReadRecord returns false.
    long int i =0;
    int k = 0;
    long int num_alignments =0;
    uint8_t *p;
    uint32_t len = 0; //length of qual seq.
    int threshold = 30.0; //qual score threshold

    std::string tags[3]; 
    while(sam_read1(fp_in, bamHdr,aln) > 0) {
       barcode = (char *)get_Ztag_or_default(aln, options.barcode_tag.c_str(), none);
       barcode_raw = (char *)get_Ztag_or_default(aln, "CR", empty);

        // corrected and raw barcodes should match to be perfect
       int perfect_cell_barcode = 1;
       
       int barcode_len = strlen(barcode);
       // empty barcodes are not perfect
       if (strcmp(barcode, "None")==0)  {
           perfect_cell_barcode = 0;     
       } else {
          for (k=0; k < barcode_len; k++) 
            if (barcode[k] != barcode_raw[k]) {
               perfect_cell_barcode = 0;     
               break;
            }
       }

       // barcode quality score
       barcode_qual = (char *)get_Ztag_or_default(aln, "CY", empty);

       float sum_barcode_qual = 0;
       float num_bp_above_threshold = 0;
       len = strlen(barcode_qual);

       for (k = 0; k < len; k++) {
          sum_barcode_qual += ((uint8_t)barcode_qual[k] -33);
          if (((uint8_t)barcode_qual[k] - 33)  > threshold) num_bp_above_threshold += 1;
       }

       int avg_cell_barcode_qual = sum_barcode_qual/(float)len;
       float cell_barcode_qual_above_threshold =  (float)num_bp_above_threshold/(float)len;

       // corrected molecule barcodes
       umi = (char *)get_Ztag_or_default(aln, options.umi_tag.c_str(), empty);

       // raw molecule barcodes
       umi_raw = (char *)get_Ztag_or_default(aln, "UR", empty);

       // corrected and raw barcodes should match to be perfect
       int perfect_molecule_barcode = 1;
       if (strlen(umi)!=strlen(umi_raw)) {  // not equal length case
          perfect_molecule_barcode = 0;     
       } else {   // if equal then compare char by char
          int umi_len = strlen(umi);
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
          sum_umi_qual += ((uint8_t)umi_qual[k] -33);
          if (((uint8_t)umi_qual[k] - 33)  > threshold) num_umi_above_threshold += 1;
       }
       float frac_umi_qual_above_threshold =  (float)num_umi_above_threshold/(float)len;

       gene_id = (char *)get_Ztag_or_default(aln, options.gene_tag.c_str(), empty);

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
           std::cout << "Batch number : " << partial_files.size() << std::endl;
           std::string split_file_path;
           
           split_file_path = write_out_partial_txt_file(tuple_records, tmp_folder);
           num_alignments += tuple_records.size();
           partial_files.push_back(split_file_path);
           tuple_records.clear();
//         if (partial_files.size() ==7) break;
       }

       tags[options.tag_order[options.barcode_tag]] = barcode;
       tags[options.tag_order[options.umi_tag]] = umi;
       tags[options.tag_order[options.gene_tag]] = gene_id;
   
       tuple_records.push_back(
            std::make_tuple( 
                             std::string(tags[0]) + std::string("\t") + std::string(tags[1]) + std::string("\t") + std::string(tags[2]), 
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
