#pragma once
#include <vector>
#include <string>

int restore_core(const std::string &smerFileName, const std::string &brcFileName, const std::string &outputFileName, uint32_t nThreads, bool useKmerFormat);