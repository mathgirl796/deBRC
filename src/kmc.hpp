#pragma once
#include <vector>
#include <string>

int kmc_core(const std::vector<std::string> &inputFiles, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads, uint32_t kmerLen, const std::string &outputFileName, bool canonical);
