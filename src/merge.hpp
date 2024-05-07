#pragma once
#include <vector>
#include <string>

int merge_core(const std::vector<std::string> &inputFiles, const std::string &outputFileName, 
    const bool &distinct, uint32_t maxRamGB) ;

int compare_core(const std::string &refFileName, const std::string &tgtFileName, 
    const std::string &outputFileName);