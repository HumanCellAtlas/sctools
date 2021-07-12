#ifndef __TAG_SORT__
#define __TAG_SORT__

#include "utilities.h"
#include "inputoptions.h"

#include <functional>
#include <string>
#include <vector>
#include <fstream>
#include <queue>
#include <tuple>
#include <cstdint>
#include <iostream>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


struct CustomCompare {
   bool operator()(const QUEUETUPLE& a, const QUEUETUPLE& b) {
     return get<0>(a) < get<0>(b);
   }
};

struct Context {
   vector<vector<std::string>> data;
   vector<ifstream *> file_handles;
   vector<int> data_size;
   vector<int> ptrs ;
   vector<bool> isempty;
   int i = 0 ;
   int num_active_files = 0;
   int BUF_SIZE = 100000;
   int NUM_PARTS = 0;
   std::priority_queue<QUEUETUPLE, std::vector<QUEUETUPLE>,  CustomCompare> heap;
};

#endif
