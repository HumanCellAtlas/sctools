#include "sort_write.h"
#include <sstream>

#define STRING_LEN  40

inline bool sortbyfirst(const std::pair<std::string, int>& a, 
                  const std::pair<std::string, int>& b) {

            return a.first > b.first;
}


using namespace std;
/** @copydoc write_out_partial_txt_file */
std::string write_out_partial_txt_file(const vector<TAGTUPLE> &tuple_records, std::string const & tmp_folder) {
    std::string tempfile = tmp_folder + string("/") + random_string(STRING_LEN) + std::string(".txt");

    ofstream output_fp;
    output_fp.open (tempfile); 
    std::vector<std::pair<std::string, int>> index_pairs;
    int k  = 0;
     
    for( auto it=tuple_records.begin(); it!=tuple_records.end(); it++, k++) { 
       index_pairs.push_back(std::make_pair(get<0>(*it) , k));
    }

    // sort using the lambda function
    std::sort(index_pairs.begin(), index_pairs.end(), sortbyfirst);
    //return "hello";

    stringstream str(stringstream::out|stringstream::binary);

    for (auto it=index_pairs.begin(); it!=index_pairs.end(); it++) {
          str << get<0>(tuple_records[it->second]) /* three tags */ << "\t"
              << get<1>(tuple_records[it->second]) /* record[0] */  << "\t"
              << get<2>(tuple_records[it->second]) /* record[1] */  << "\t" 
              << get<3>(tuple_records[it->second]) /* record[2] */  << "\t" 
              << get<4>(tuple_records[it->second]) /* record[3] */  << "\t" 
              << get<5>(tuple_records[it->second]) /* record[4] */  << "\t" 
              << get<6>(tuple_records[it->second]) /* record[5] */  << "\t" 
              << get<7>(tuple_records[it->second]) /* record[6] */  << "\t" 
              << get<8>(tuple_records[it->second]) /* record[7] */  << "\t" 
              << get<9>(tuple_records[it->second]) /* record[8] */  << "\t" 
              << get<10>(tuple_records[it->second]) /* record[9] */ << "\t" 
              << get<11>(tuple_records[it->second]) /* record[10] */ << "\t"
              << get<12>(tuple_records[it->second]) /* record[11] */ << "\t"
              << get<13>(tuple_records[it->second]) /* record[12] */ << "\t"
              << get<14>(tuple_records[it->second]) /* record[13] */
              << std::endl; 
    }
    // what is you ran out of disk space ???? NEED TO add logic
    output_fp.write(str.str().c_str(), str.str().length());
    output_fp.close();
    return tempfile;
}


