#pragma once
#include <vector>
#include <string>

int split_core(const std::string &inputFileName, const std::string &outputFileName, const bool ismer, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads);