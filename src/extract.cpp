#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include "utils.hpp"
#include "FastaReader/FastaReader.hpp"
extern unsigned char nst_nt4_table[256];

int extract_core(const std::string &inputFileName, const std::string &outputFileName, uint32_t kmerLength) {
    FastaReader fr(inputFileName);
    unordered_set<uint64_t> kmerSet;
    uint32_t k = kmerLength;
    while (fr.hasNextSequence()) {
        fr.readSequence();
        string seq = fr.Sequence();
        err_func_printf(__func__, "processing %s\n", fr.Id().c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
        progressbar bar(seq.length() - k + 1);
        uint64_t kmer = 0;
        for (size_t i = 0; i < seq.length() - k + 1; ++i) {
            // 循环处理输入fa中的一条sequence的每个kmer
            if (i == 0) {
                for (uint32_t j = 0; j < k - 1; ++j) {
                    kmer = (kmer << 2) | nst_nt4_table[(int)seq[j]];
                }
            }
            // 第i个kmer
            kmer = ((kmer << 2) | nst_nt4_table[(int)seq[i + k - 1]]) & (~(((uint64_t)-1) << (k << 1)));
            kmerSet.insert(kmer);
            bar.update();
        }
        bar.end();
    }

    err_func_printf(__func__, "writing to .smer file...\n");
    FILE *outputFile = xopen((outputFileName + ".mer").c_str(), "wb"); 
    uint64_t kmerCount = kmerSet.size();
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    err_fwrite(&k, sizeof(uint32_t), 1, outputFile);
    err_fwrite(&kmerCount, sizeof(uint64_t), 1, outputFile);
    for(unordered_set<uint64_t>::iterator kmer = kmerSet.begin(); kmer != kmerSet.end(); ++kmer) {
        err_fwrite(&(*kmer), sizeof(uint64_t), 1, outputFile);
    }
    err_fclose(outputFile);
    return 0;
}
