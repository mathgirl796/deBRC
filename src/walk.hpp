#pragma once
#include <vector>
#include <string>

int walk_core(const std::string &smerFileName, const std::string &okFileName, 
    const std::vector<std::string> &faFiles, const std::string &outputFileName, uint32_t nThreads,
    bool passSpecialCharactors, bool useKmerFormat, int maxRamGB);