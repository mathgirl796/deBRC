#include <string>
#include <vector>
#include <iostream>
#include <pthread.h>
#include <algorithm>
#include "utils.hpp"
#include "klib/kthread.hpp"
#include "merge.hpp"

using namespace std;

struct SortData {
    FILE* inputFile;
    pthread_mutex_t inputFileLock;
    uint64_t kmerLength;
    uint64_t kmerCount;
    uint64_t singleBufferKmerCount;
    string tmpPath;
    uint64_t totalTaskNum;
    pthread_mutex_t tmpFileLock;
    vector<string> tmpFileNames;
};

static int compare_uint64 (const void *ptr_a, const void *ptr_b) {
    uint64_t a = *(uint64_t *)ptr_a, b = *(uint64_t *)ptr_b;
    if (a < b) return -1;
    else if (a == b) return 0;
    else return 1;
} 

void ktf_sort(void* data, long i, int tid) {
    // 提取要排序的数据范围
    SortData *d = (SortData *)data;
    uint64_t startCount = i * d->singleBufferKmerCount;
    uint64_t endCount = ((i + 1) * d->singleBufferKmerCount < d->kmerCount) ? ((i + 1) * d->singleBufferKmerCount) : (d->kmerCount);
    uint64_t bufferCount = endCount - startCount;
    uint64_t *numbers = (uint64_t *)err_calloc(__func__, bufferCount, sizeof(uint64_t));
    // 从文件中读取相应数量的数据到内存中
    pthread_mutex_lock(&d->inputFileLock);
    err_fread_noeof(numbers, sizeof(uint64_t), bufferCount, d->inputFile);
    pthread_mutex_unlock(&d->inputFileLock);
    // 开始快速排序
    err_func_printf(__func__, "tid:%d, i:%ld, bufferCount:%lu, from %lu to %lu, start.\n", tid, i, bufferCount, startCount, endCount);
    qsort(numbers, bufferCount, sizeof(uint64_t), compare_uint64);
    // 输出结果到临时文件中
    char fn[CommonFileNameBufSize];
    sprintf(fn, "%s/%ld.smer", d->tmpPath.c_str(), i);
    pthread_mutex_lock(&d->tmpFileLock);
    d->tmpFileNames.push_back(string(fn));
    pthread_mutex_unlock(&d->tmpFileLock);
    FILE *tmpWrite = xopen(fn, "wb");
    err_fwrite(&d->kmerLength, sizeof(uint32_t), 1, tmpWrite);
    err_fwrite(&bufferCount, sizeof(uint64_t), 1, tmpWrite);
    err_fwrite(numbers, sizeof(uint64_t), bufferCount, tmpWrite);
    free(numbers);
    err_fclose(tmpWrite);
    err_func_printf(__func__, "tid:%d done.\n", tid);
}

/**
 * 首先打开输入文件，读取头部获取里面整数的个数；
 * 然后根据用户规定的maxRamGB，将文件划分成占用存储不大于maxRamGB/nThreads的块
 * 然后用kthread库中的kt_for创建线程池，对这些块进行排序
 * 最后创建输出文件，写入头部，进而用堆排序将这些块归并，写入数据部分
*/
int sort_core(const std::string &inputFileName, const std::string &outputFileName, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads, bool distinct) {
    // 打开待排序文件
    FILE* inputFile = xopen(inputFileName.c_str(), "r");
    uint32_t kmerLength;
    uint64_t kmerCount;
    err_fread_noeof(&kmerLength, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kmerCount, sizeof(uint64_t), 1, inputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", kmerLength, kmerCount);
    // 准备多线程排序所需数据
    SortData data;
    data.inputFile = inputFile;
    data.inputFileLock = PTHREAD_MUTEX_INITIALIZER;
    data.tmpFileLock = PTHREAD_MUTEX_INITIALIZER;
    data.kmerLength = kmerLength;
    data.kmerCount = kmerCount;
    data.singleBufferKmerCount = (uint64_t)maxRamGB * OneGiga / nThreads / sizeof(uint64_t);
    xassert(data.singleBufferKmerCount > 0, "maxRamGB too small or nThreads too large!");
    data.tmpPath = tmpPath;
    err_func_printf(__func__, "singleBufferKmerCount:[%lu], tmpPath:[%s]\n", data.singleBufferKmerCount, tmpPath.c_str());
    data.totalTaskNum = (data.kmerCount - 1) / data.singleBufferKmerCount + 1;
    err_func_printf(__func__, "totalTaskNum:%lu\n", data.totalTaskNum);
    // 执行多线程排序，结果写入到临时文件中
    err_func_printf(__func__, "sorting... (total %lu tasks)\n", data.totalTaskNum);
    kt_for(nThreads, ktf_sort, &data, data.totalTaskNum);
    err_func_printf(__func__, "done sorting\n");
    // 多路归并
    err_func_printf(__func__, "merging...\n");
    merge_core(data.tmpFileNames, outputFileName, false, maxRamGB); 
    err_func_printf(__func__, "done merging to %s\n", outputFileName.c_str());
    // 删除临时文件
    err_func_printf(__func__, "deleting...\n");
    for (auto tmpFileName : data.tmpFileNames) {
        if (remove(tmpFileName.c_str()) != 0) {
            err_func_printf(__func__, "error delete temp file [%s]\n", tmpFileName.c_str());
        }
    }
    err_func_printf(__func__, "done deleting temp files\n");
    return 0;
}

