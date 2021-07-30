#ifndef __TAG_SORT__
#define __TAG_SORT__

#include "utilities.h"
#include "datatypes.h"
#include "inputoptions.h"
#include "gzstream.h"

#include <functional>
#include <queue>

#include <string>
#include <vector>
#include <fstream>
#include <tuple>
#include <cstdint>
#include <iostream>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


struct Context {
   vector<vector<std::string>> data;
#ifdef GZSTREAM
   vector<igzstream *> file_handles;
#else
   vector<ifstream *> file_handles;
#endif

   vector<int> data_size;
   vector<int> ptrs;
   vector<bool> isempty;
   int i = -1;
   int num_active_files;
   int BUF_SIZE;
   int NUM_PARTS;

   Context(unsigned int num_parts, int data_buffer_size): 
        NUM_PARTS(num_parts), BUF_SIZE(data_buffer_size) {

       num_active_files = 0;

       // set the isempty for each file to false
       for (auto i=0; i < this->NUM_PARTS; i++) {
          this->isempty.push_back(false);
       }
   
       // set a vector of vectors of data for each file
       for (auto i=0; i < this->NUM_PARTS; i++) {
          this->data.push_back(std::vector<std::string>());
       }
   
       // set the data_size of the buffer for each file to 0
       for (auto i=0; i < this->NUM_PARTS; i++) {
          this->data_size.push_back(0);
       }
   
       // set the pointer to f each buffer to thist.BUF_SIZE
       for (auto i=0; i < this->NUM_PARTS; i++) {
          this->ptrs.push_back(this->BUF_SIZE);
       }
   }

   void print_status() {
      std::cout << "Contx status " << std::endl;
      for (auto i=0; i < this->NUM_PARTS; i++) {
         this->i = i;
         std::cout << "\t" << this->i << "\t" << this->data[this->i].size() << "\t" << this->data_size[this->i] << "\t"  << this->ptrs[this->i] << std::endl;
      }
   }


   void clear() {
       for(auto it = file_handles.begin(); it != file_handles.end(); ++it) {
          delete *it; 
       }
       data_size.clear();
       ptrs.clear();
       isempty.clear();
   }
};

#endif
