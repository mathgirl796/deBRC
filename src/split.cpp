#include <vector>
#include <string>
#include <set>
#include <iostream>
#include "utils.hpp"
using namespace std;

/**
 * 给定一个文件，文件格式如下：
 * uint32_t kmerLength | uint64_t kmerCount | uint64_t kmer (kmerCount个，从小到大有序排列)
 * 
 * 假定输入文件的kmerLength为k+1
 * 输出两个文件outputFileName.smer和outputFileName.o.smer
 * outputFileName.o.smer中包含所有多出kmer的子kmer（不允许重复值）
 * 
 * ksmer: 需要.k.smer
*/
int split_core(const std::string &inputFileName, const std::string &outputFileName, const bool ksmer) {
    
    FILE* inputFile = xopen(inputFileName.c_str(), "r");
    setvbuf(inputFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t kp1;
    uint64_t kp1merCount;
    err_fread_noeof(&kp1, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kp1merCount, sizeof(uint64_t), 1, inputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", kp1, kp1merCount);

    // 分出kmerSet和okmerSet
    uint32_t k = kp1 - 1;
    set<uint64_t> kmerSet;
    set<uint64_t> okmerSet;
    uint64_t lastLkmer, lastRkmer;
    uint64_t kp1mer;
    progressbar bar(kp1merCount);
    for (uint64_t i = 0; i < kp1merCount; ++i) {
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, inputFile);
        uint64_t lkmer = kp1mer >> 2;
        uint64_t rkmer = kp1mer & (~(0x3 << ((kp1 << 1) - 2)));
        if (ksmer) {
            kmerSet.insert(lkmer);
            kmerSet.insert(rkmer);
        }
        if(i > 0 && lkmer == lastLkmer) {
            okmerSet.insert(lastRkmer);
            okmerSet.insert(rkmer);
        }
        lastLkmer = lkmer;
        lastRkmer = rkmer;
        bar.update();
    }
    bar.end();
    uint64_t kmerCount = kmerSet.size();
    uint64_t okmerCount = okmerSet.size();

    // 分别写入到文件
    string smerFileName = outputFileName + ".k.smer";
    string osmerFileName = outputFileName + ".o.smer";
    FILE* outputFile = NULL;

    if (ksmer){ // 写入.k.smer
        err_func_printf(__func__, "writing to .k.smer file...\n");
        outputFile = xopen(smerFileName.c_str(), "wb");
        setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
        err_fwrite(&k, sizeof(uint32_t), 1, outputFile);
        err_fwrite(&kmerCount, sizeof(uint64_t), 1, outputFile);
        for(set<uint64_t>::iterator kmer = kmerSet.begin(); kmer != kmerSet.end(); ++kmer) {
            err_fwrite(&(*kmer), sizeof(uint64_t), 1, outputFile);
        }
        err_fclose(outputFile);
    }

    err_func_printf(__func__, "writing to .o.smer file...");
    outputFile = xopen(osmerFileName.c_str(), "wb"); // 写入.o.smer
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    err_fwrite(&k, sizeof(uint32_t), 1, outputFile);
    err_fwrite(&okmerCount, sizeof(uint64_t), 1, outputFile);
    for(set<uint64_t>::iterator kmer = okmerSet.begin(); kmer != okmerSet.end(); ++kmer) {
        err_fwrite(&(*kmer), sizeof(uint64_t), 1, outputFile);
    }
    err_fclose(outputFile);

    return 0;
}


// TODO: 先分成mer，再排序的方式（省内存）
// 难点：mer文件的去重
// 解决方案：写一个函数，读入smer文件，对它进行去重，复写到原来的文件上。
// 紧急程度：比较低，应该先写生成原序列brc的代码