#include <vector>
#include <string>
#include <map>
#include "utils.hpp"
#include "klib/kthread.hpp"
#include "FastaReader/FastaReader.hpp"

int restore_core(const std::string &kFileName, const std::string &brcFileName, const std::string &outputFileName, uint32_t nThreads) {
    FILE *kFile = xopen(kFileName.c_str(), "rb");
    setvbuf(kFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t k = 0;
    uint64_t kmerNum = 0;
    err_fread_noeof(&k, sizeof(uint32_t), 1, kFile);
    err_fread_noeof(&kmerNum, sizeof(uint64_t), 1, kFile);
    map<uint64_t, uint8_t> kmerMap; // uint8_t的高四位无意义，低四位由低到高表示ACGT
    for (uint64_t i = 0; i < kmerNum; ++i) {
        uint64_t kmer;
        err_fread_noeof(&kmer, sizeof(uint64_t), 1, kFile);
        kmerMap[kmer >> 2] |= (0x1 << ((kmer & 0x3) + 1));
    }

    FastaReader brcFile(brcFileName);

    //TODO:everything

    return 0;
}