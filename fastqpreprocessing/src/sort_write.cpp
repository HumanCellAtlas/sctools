#include "sort_write.h"
#include <sstream>

#define STRING_LEN  40

inline bool sortbyfirst(const std::pair<std::string, int>& a, 
                  const std::pair<std::string, int>& b) {
    return (a.first < b.first);
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
       // 'UB:CB:GX"   or "GX:CB:UB"
       index_pairs.push_back(std::make_pair(get<0>(*it) + std::string("\t") + get<1>(*it) + std::string("\t") + get<2>(*it) , k));
    }

    stringstream str(stringstream::out|stringstream::binary);
    // sort using the lambda function
    std::sort(index_pairs.begin(), index_pairs.end(), sortbyfirst);

    for (auto it=index_pairs.begin(); it!=index_pairs.end(); it++) {
       str << get<0>(tuple_records[it->second]) << "\t"
                 << get<1>(tuple_records[it->second]) << "\t"
                 << get<2>(tuple_records[it->second]) << std::endl; 
    }

    output_fp.write(str.str().c_str(), str.str().length());
    output_fp.close();

    return tempfile;
}


