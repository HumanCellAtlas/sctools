#include "sort_write.h"
#include <sstream>

#define STRING_LEN  40

extern sem_t semaphore;
#define SEM_POST(X)                                      \
 ({                                                      \
     if (sem_post(&X) == -1)                             \
        error("sem_post: semaphore");                   \
 })

inline bool sortbyfirst(const std::pair<TRIPLET *, int>& a, 
                  const std::pair<TRIPLET *, int>& b) {
    if ((*get<0>(*a.first)).compare(*get<0>(*b.first)) !=0) {
       return ((*get<0>(*a.first)).compare(*get<0>(*b.first)) < 0);
    } else { 
       if ((*get<1>(*a.first)).compare(*get<1>(*b.first)) !=0) {
           return ((*get<1>(*a.first)).compare(*get<1>(*b.first)) < 0);
       }
    } 
    return ((*get<2>(*a.first)).compare(*get<2>(*b.first)) < 0);
}


using namespace std;
/** @copydoc write_out_partial_txt_file */
void  write_out_partial_txt_file(const vector<TAGTUPLE> &tuple_records, \
              std::string const & tmp_folder,  std::vector<string> &partial_files) {;

    std::string tempfile = tmp_folder + string("/") + random_string(STRING_LEN) + std::string(".txt.gz");

#ifdef GZSTREAM
    ogzstream  output_fp;
#else
    ofstream output_fp;
#endif

    output_fp.open(tempfile.c_str());

    std::vector<std::pair<TRIPLET *, int>> index_pairs;
    int k  = 0;
     
    for(auto it=tuple_records.begin(); it!=tuple_records.end(); it++, k++) { 
        index_pairs.push_back(std::make_pair(get<0>(*it) , k));
    }

    // sort using the lambda function
    std::sort(index_pairs.begin(), index_pairs.end(), sortbyfirst);

    stringstream str(stringstream::out|stringstream::binary);

    for (auto it=index_pairs.begin(); it!=index_pairs.end(); it++, k++) {
          
          // what is you ran out of disk space ???? NEED TO add logic
          if (k%10000==0) {
             output_fp.write(str.str().c_str(), str.str().length());
             str.str("");
             str.clear();
          }

          str << *get<0>(*it->first) /*  frist tag */ << "\t"
              << *get<1>(*it->first) /*  second tag */ << "\t"
              << *get<2>(*it->first) /*  third tag */ << "\t"
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
    // what if you ran out of disk space how do you inform the user? 
    output_fp.write(str.str().c_str(), str.str().length());

    str.str("");
    str.clear();
    output_fp.close();
    partial_files.push_back(tempfile);

    SEM_POST(semaphore);
}


