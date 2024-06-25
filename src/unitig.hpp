#pragma once
#include <vector>
#include <string>

int unitig_core(const std::string &ikFileName, const std::string &okFileName, 
    const std::string &faFileName, const std::string &outputFileName, 
    uint32_t nThreads, int maxRamGB);