#pragma once
#include <vector>
#include <string>

int sort_core(const std::string &inputFileName, const std::string &outputFileName, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads, bool distinct, bool usmer);