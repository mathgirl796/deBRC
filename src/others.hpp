#pragma once
#include <string>

int check_core(const std::string &inputFileName, uint32_t nThreads, bool check_dup);
int view_core(const std::string &inputFileName, bool head, bool tail) ;
int search_core(const std::string &smerFileName, std::string kmerStr) ;