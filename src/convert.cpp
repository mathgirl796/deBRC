#include "kmc_api/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "utils.hpp"

int convert_core(const std::string &inputFileName, const std::string &outputFileName)
{
    CKMCFile kmer_database; // it holds 'uint64 KmerCount(void)' and 'uint32 KmerLength(void);'
    bool result = kmer_database.OpenForListing(inputFileName);
    if (!result) {
        err_func_printf(__func__, "read kmer_database failed");
        return 1;
    }
    uint32_t kmerLength = kmer_database.KmerLength();
    uint64_t kmerCount = kmer_database.KmerCount();
    if (kmer_database.KmerLength() > 32) {
        err_func_printf(__func__, "cannot convert kmer whose k > 32");
        return 1;
    }
    CKmerAPI kmer_object(kmerLength); 
    uint32_t counter; 
    std::vector<uint64> kmer_long;
    std::string fullOutputFileName = outputFileName + ".mer";
    FILE* outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    err_fwrite(&kmerLength, sizeof(uint32_t), 1, outputFile);
    err_fwrite(&kmerCount, sizeof(uint64_t), 1, outputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", kmerLength, kmerCount);
    progressbar bar(kmerCount);
    while(kmer_database.ReadNextKmer(kmer_object, counter)) {
        kmer_object.to_long(kmer_long);
        err_fwrite(&kmer_long[0], sizeof(uint64_t), 1, outputFile);
        bar.update();
    }
    bar.end();
    kmer_database.Close();
    err_fclose(outputFile);
    err_func_printf(__func__, "done saving to %s\n", fullOutputFileName.c_str());
    return 0;
}