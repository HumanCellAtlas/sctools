#ifndef __SCTOOL_GLOBALS__
#define __SCTOOL_GLOBALS__

/**
 *  @file   globals.h
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <semaphore.h>
#include <mutex>
#include <vector>
#include <string>
#include <set>

#define THRESHOLD 30.0
#define NUM_THREADS 10
#define MAX_THREADS 30

extern sem_t semaphore;
extern std::mutex mtx;
extern std::vector<std::string> partial_files;

extern std::set<unsigned int> busy_buffers, idle_buffers, threads_to_join;

#endif
