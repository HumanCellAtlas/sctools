#include "globals.h"
sem_t semaphore;
std::mutex mtx;
std::vector<std::string> partial_files;

std::set<unsigned int> busy_buffers, idle_buffers, threads_to_join;


