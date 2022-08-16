/**
 *  @file   globals.cpp
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/
#include "globals.h"

sem_t semaphore;
std::mutex mtx;

// TODO figure out exactly how this is being used, and clean it up
std::vector<std::string> partial_files;

std::set<unsigned int> busy_buffers, idle_buffers, threads_to_join;
