#pragma once
#include <vector>
#include <string>

int split_core(const std::string &inputFileName, const std::string &outputFileName, const bool ksmer, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads);
int split_core_legacy(const std::string &inputFileName, const std::string &outputFileName, const bool ksmer) ;