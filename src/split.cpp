#include <vector>
#include <string>
#include <set>
#include <iostream>
#include "utils.hpp"
#include "sort.hpp"
#include "klib/kthread.hpp"
#include "merge.hpp"
using namespace std;

static int compare_uint64 (const void *ptr_a, const void *ptr_b) {
    uint64_t a = *(uint64_t *)ptr_a, b = *(uint64_t *)ptr_b;
    if (a < b) return -1;
    else if (a == b) return 0;
    else return 1;
} 

struct SplitData {
    FILE* inputFile;
    uint64_t inputBatchNum;
    pthread_mutex_t inputFileLock;
    uint64_t kp1;
    uint64_t kmerCount;
    uint64_t singleBufferKmerCount;
    uint64_t totalTaskNum;
    vector<string> tmpMerFileNames;
    vector<string> tmpOmerFileNames;
    bool ksmer;
};

void ktf_split(void* data, long i, int tid) {
    // 提取要处理的数据范围
    SplitData *d = (SplitData *)data;
    pthread_mutex_lock(&d->inputFileLock);
    i = (long)d->inputBatchNum;
    uint64_t startCount = i * d->singleBufferKmerCount;
    startCount = (startCount >= 4) ? startCount - 4 : startCount; // 分任务的起始位置应该适当提前（提前4个应该纯纯够了）
    uint64_t endCount = ((i + 1) * d->singleBufferKmerCount < d->kmerCount) ? ((i + 1) * d->singleBufferKmerCount) : (d->kmerCount);
    uint64_t bufferCount = endCount - startCount;
    uint64_t *kp1mers = (uint64_t *)err_calloc(__func__, bufferCount, sizeof(uint64_t));
    // 从文件中读取相应数量的数据到内存中
    err_fread_noeof(kp1mers, sizeof(uint64_t), bufferCount, d->inputFile);
    err_fseek(d->inputFile, -sizeof(uint64_t) * 4, SEEK_CUR); // 给下一个任务提前的四个kp1mer准备好前移过的文件指针
    d->inputBatchNum ++;
    pthread_mutex_unlock(&d->inputFileLock);
    // 把kp1mers分成kmer和okmer存入临时文件
    err_func_printf(__func__, "worker:%d, i:%ld, bufferCount:%lu, from %lu to %lu, start.\n", tid, i, bufferCount, startCount, endCount);
    FILE *merFile = xopen(d->tmpMerFileNames[i].c_str(), "wb+");
    FILE *omerFile = xopen(d->tmpOmerFileNames[i].c_str(), "wb+");
    setvbuf(merFile, NULL, _IOFBF, CommonFileBufSize);
    setvbuf(omerFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t k  = d->kp1 - 1, kp1 = d->kp1;
    uint64_t merCount = 0, omerCount = 0;
    err_fwrite(&k, sizeof(uint32_t), 1, merFile);
    err_fwrite(&k, sizeof(uint32_t), 1, omerFile);
    err_fwrite(&merCount, sizeof(uint64_t), 1, merFile); // 暂时写入mer个数，占位
    err_fwrite(&omerCount, sizeof(uint64_t), 1, omerFile);
    bool ksmer = d->ksmer;
    uint64_t lastLkmer = 0, lastRkmer = 0;
    for (uint64_t i = 0; i < bufferCount; ++i) {
        uint64_t kp1mer = kp1mers[i];
        uint64_t lkmer = kp1mer >> 2;
        uint64_t rkmer = kp1mer & (~(0x3 << ((kp1 << 1) - 2)));
        if (ksmer) {
            err_fwrite(&lkmer, sizeof(uint64_t), 1, merFile);
            err_fwrite(&rkmer, sizeof(uint64_t), 1, merFile);
            merCount += 2;
        }
        if(i > 0 && lkmer == lastLkmer) {
            err_fwrite(&lastRkmer, sizeof(uint64_t), 1, omerFile);
            err_fwrite(&rkmer, sizeof(uint64_t), 1, omerFile);
            omerCount += 2;
        }
        lastLkmer = lkmer;
        lastRkmer = rkmer;
    }

    free(kp1mers);

    err_fseek(merFile, sizeof(uint32_t), SEEK_SET); // 写入实际mer个数
    err_fwrite(&merCount, sizeof(uint64_t), 1, merFile);
    err_fseek(omerFile, sizeof(uint32_t), SEEK_SET);
    err_fwrite(&omerCount, sizeof(uint64_t), 1, omerFile);

    // 读出刚刚提取的mer数据进行排序
    if (ksmer) {
        err_fseek(merFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); 
        uint64_t *kmers = (uint64_t *)err_calloc(__func__, merCount, sizeof(uint64_t));
        err_fread_noeof(kmers, sizeof(uint64_t), merCount, merFile);
        qsort(kmers, merCount, sizeof(uint64_t), compare_uint64);
        err_fseek(merFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); // 写入实际mer个数
        err_fwrite(kmers, sizeof(uint64_t), merCount, merFile);
        free(kmers);
    }
    err_fseek(omerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); 
    uint64_t *kmers = (uint64_t *)err_calloc(__func__, omerCount, sizeof(uint64_t));
    err_fread_noeof(kmers, sizeof(uint64_t), omerCount, omerFile);
    qsort(kmers, omerCount, sizeof(uint64_t), compare_uint64);
    err_fseek(omerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); // 写入实际mer个数
    err_fwrite(kmers, sizeof(uint64_t), omerCount, omerFile);
    free(kmers);

    err_fclose(merFile); // 关闭文件
    err_fclose(omerFile);

    err_func_printf(__func__, "worker:%d done.\n", tid);
}

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
int split_core(const std::string &inputFileName, const std::string &outputFileName, const bool ksmer, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads) {
    
    FILE* inputFile = xopen(inputFileName.c_str(), "r");
    setvbuf(inputFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t kp1;
    uint64_t kp1merCount;
    err_fread_noeof(&kp1, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kp1merCount, sizeof(uint64_t), 1, inputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", kp1, kp1merCount);

    SplitData data;
    data.inputFile = inputFile;
    data.inputBatchNum = 0;
    data.inputFileLock = PTHREAD_MUTEX_INITIALIZER;
    data.kp1 = kp1;
    data.kmerCount = kp1merCount;
    data.singleBufferKmerCount = (uint64_t)maxRamGB * OneGiga / nThreads / sizeof(uint64_t) / 2; // 最后一个除2是调整内存用量
    // data.singleBufferKmerCount = 1000000;
    xassert(data.singleBufferKmerCount > 0, "maxRamGB too small or nThreads too large!");
    data.totalTaskNum = (data.kmerCount - 1) / data.singleBufferKmerCount + 1;
    for (uint64_t i = 0; i < data.totalTaskNum; ++i) {
        data.tmpMerFileNames.push_back(tmpPath + string_format("/%lu.k.mer.tmp", i));
        data.tmpOmerFileNames.push_back(tmpPath + string_format("/%lu.o.mer.tmp", i));
    }
    data.ksmer = ksmer;

    // 执行多线程分类+排序，结果写入到临时文件中
    err_func_printf(__func__, "spliting and sorting... (total %lu tasks)\n", data.totalTaskNum);
    kt_for(nThreads, ktf_split, &data, data.totalTaskNum);
    err_func_printf(__func__, "done spliting and sorting\n");
    // 多路归并
    err_func_printf(__func__, "merging...\n");
    string fullOutputFileName;
    if (ksmer) {
        fullOutputFileName = outputFileName + ".k";
        merge_core(data.tmpMerFileNames, fullOutputFileName, true, maxRamGB); 
        err_func_printf(__func__, "done merging to %s\n", fullOutputFileName.c_str());
    }
    fullOutputFileName = outputFileName + ".o";
    merge_core(data.tmpOmerFileNames, fullOutputFileName, true, maxRamGB); 
    err_func_printf(__func__, "done merging to %s\n", fullOutputFileName.c_str());

    // 删除临时文件
    err_func_printf(__func__, "deleting...\n");
    for (auto tmpFileName : data.tmpMerFileNames) {
        if (remove(tmpFileName.c_str()) != 0) {
            err_func_printf(__func__, "error delete temp file [%s]\n", tmpFileName.c_str());
        }
    }
    for (auto tmpFileName : data.tmpOmerFileNames) {
        if (remove(tmpFileName.c_str()) != 0) {
            err_func_printf(__func__, "error delete temp file [%s]\n", tmpFileName.c_str());
        }
    }
    err_func_printf(__func__, "done deleting temp files\n");

    return 0;
}


int split_core_legacy(const std::string &inputFileName, const std::string &outputFileName, const bool ksmer) {
    
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