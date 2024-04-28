#pragma once
#include <vector>
#include <string>

/**
 * extract all distinct kmer from inputFile to outputFile
 * outputFile无需扩展名，会自动添加.mer扩展名
*/
int extract_core(const std::string &inputFileName, const std::string &outputFileName, uint32_t kmerLength);