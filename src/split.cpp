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
    uint64_t k;
    uint64_t kmerCount;
    uint64_t singleBufferKmerCount;
    uint64_t totalTaskNum;
    vector<string> tmpOmerFileNames;

    bool usmer;
    vector<string> tmpImerFileNames;
    
};

void ktf_split(void* data, long i, int tid) {
    // 提取要处理的数据范围
    SplitData *d = (SplitData *)data;
    pthread_mutex_lock(&d->inputFileLock);
    i = (long)d->inputBatchNum;
    uint64_t startCount = i * d->singleBufferKmerCount;
    uint64_t overlapCount = (d->usmer)?16:4; // 分任务的起始位置应该适当提前（smer形式需要提前4个，usmer形式要提前16个）
    startCount = (startCount >= overlapCount) ? (startCount - overlapCount) : startCount; 
    uint64_t endCount = ((i + 1) * d->singleBufferKmerCount < d->kmerCount) ? ((i + 1) * d->singleBufferKmerCount) : (d->kmerCount);
    uint64_t bufferCount = endCount - startCount;
    uint64_t *kpmers = (uint64_t *)err_calloc(__func__, bufferCount, sizeof(uint64_t));
    
    err_func_printf(__func__, "worker:%d, i:%ld, bufferCount:%lu, from %lu to %lu, start.\n", tid, i, bufferCount, startCount, endCount);

    // 从文件中读取相应数量的数据到内存中
    err_fread_noeof(kpmers, sizeof(uint64_t), bufferCount, d->inputFile);
    err_fseek(d->inputFile, -sizeof(uint64_t) * overlapCount, SEEK_CUR); // 给下一个任务提前的四个kp1mer准备好前移过的文件指针
    d->inputBatchNum ++;
    pthread_mutex_unlock(&d->inputFileLock);

    // 把kp1mers中的omer存入临时文件
    uint32_t kp1;
    uint32_t k;
    FILE *omerFile;
    FILE *imerFile;
    uint64_t omerCount = 0;
    uint64_t imerCount = 0;

    kp1 = d->k; k = kp1 - 1;
    if (d->usmer) {
        kp1 = d->k - 1; k = kp1 - 1;
    }

    omerFile = xopen(d->tmpOmerFileNames[i].c_str(), "wb+");
    setvbuf(omerFile, NULL, _IOFBF, CommonFileBufSize);
    err_fwrite(&kp1, sizeof(uint32_t), 1, omerFile);
    err_fwrite(&omerCount, sizeof(uint64_t), 1, omerFile);

    if (d->usmer) {
        imerFile = xopen(d->tmpImerFileNames[i].c_str(), "wb+");
        setvbuf(imerFile, NULL, _IOFBF, CommonFileBufSize);
        err_fwrite(&kp1, sizeof(uint32_t), 1, imerFile);
        err_fwrite(&imerCount, sizeof(uint64_t), 1, imerFile);
    }

    // 分析读入的kmer，写入到临时文件
    if (!d->usmer) {
        uint64_t lastLkmer = 0;
        for (uint64_t i = 0; i < bufferCount; ++i) {
            uint64_t kp1mer = kpmers[i];
            uint64_t lkmer = kp1mer >> 2;
            if(i > 0 && lkmer == lastLkmer) {
                err_fwrite(&(kpmers[i-1]), sizeof(uint64_t), 1, omerFile);
                err_fwrite(&(kpmers[i]), sizeof(uint64_t), 1, omerFile);
                omerCount += 2;
            }
            lastLkmer = lkmer;
        }
    }
    else {
        // TODO: 暂时通过test k4肉眼检测，等待walk和restore的实验验证
        uint64_t kmer = 0;
        uint64_t i = 0;
        uint64_t j = 0;
        while (1) {
            // 解析当前kmer的出入base有哪些
            kmer = kpmers[i] >> 4;
            bool ibase[4] = {false, false, false, false};
            bool obase[4] = {false, false, false, false};
            while (1) {
                uint64_t tmpKp2mer = kpmers[j];
                if (j >= bufferCount || kmer != (tmpKp2mer >> 4)) {
                    break;
                }
                ibase[tmpKp2mer & ((uint64_t)0x3)] = true;
                obase[(tmpKp2mer >> 2) & ((uint64_t)0x3)] = true;
                j += 1;
            }
            i = j;
            // 把当前kmer的出入base情况以kp1mer形式存入临时文件
            bool multiIn = (ibase[0] + ibase[1] + ibase[2] + ibase[3]) > 1;
            bool multiOut = (obase[0] + obase[1] + obase[2] + obase[3]) > 1;
            for (uint64_t base = 0; base < 4; ++base) {
                if (multiIn && ibase[base]) {
                    uint64_t tmpKp1mer = kmer | (base << (k<<1));
                    err_fwrite(&tmpKp1mer, sizeof(uint64_t), 1, imerFile);
                    imerCount += 1;
                }
                if (multiOut && obase[base]) {
                    uint64_t tmpKp1mer = (kmer << 2) | base;
                    err_fwrite(&tmpKp1mer, sizeof(uint64_t), 1, omerFile);
                    omerCount += 1;
                }
            }
            // 判断是否结束
            if (i >= bufferCount) {
                break;
            }
        }
    }

    free(kpmers);

    err_fseek(omerFile, sizeof(uint32_t), SEEK_SET);
    err_fwrite(&omerCount, sizeof(uint64_t), 1, omerFile); // 写入实际mer个数
    if (d->usmer) {
        err_fseek(imerFile, sizeof(uint32_t), SEEK_SET);
        err_fwrite(&imerCount, sizeof(uint64_t), 1, imerFile); // 写入实际mer个数
    }

    // 读出刚刚提取的mer数据进行排序
    err_fseek(omerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); // 跳过头部读出数据
    uint64_t *kmers = (uint64_t *)err_calloc(__func__, omerCount, sizeof(uint64_t));
    err_fread_noeof(kmers, sizeof(uint64_t), omerCount, omerFile);
    qsort(kmers, omerCount, sizeof(uint64_t), compare_uint64);
    err_fseek(omerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); // 跳过头部写入数据
    err_fwrite(kmers, sizeof(uint64_t), omerCount, omerFile);
    free(kmers);

    if (d->usmer) {
        err_fseek(imerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); 
        kmers = (uint64_t *)err_calloc(__func__, imerCount, sizeof(uint64_t));
        err_fread_noeof(kmers, sizeof(uint64_t), imerCount, imerFile);
        qsort(kmers, imerCount, sizeof(uint64_t), compare_uint64);
        err_fseek(imerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET);
        err_fwrite(kmers, sizeof(uint64_t), imerCount, imerFile);
        free(kmers);
    }

    // 关闭文件
    err_fclose(omerFile);
    if (d->usmer) {err_fclose(imerFile);}

    err_func_printf(__func__, "worker:%d done.\n", tid);
}

/**
 * 给定一个文件，文件格式如下：
 * uint32_t kmerLength | uint64_t kmerCount | uint64_t kmer (kmerCount个，从小到大有序排列)
 * 
 * 假定输入文件的kmerLength为k+1
 * 输出outputFileName.o.smer
 * outputFileName.o.smer中包含所有多出kmer的子kmer和该多出kmer组成的kp1mer（不允许重复值）
*/
int split_core(const std::string &inputFileName, const std::string &outputFileName, const bool usmer, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads) {
    
    FILE* inputFile = xopen(inputFileName.c_str(), "r");
    setvbuf(inputFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t k;
    uint64_t kmerCount;
    err_fread_noeof(&k, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kmerCount, sizeof(uint64_t), 1, inputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", k, kmerCount);

    SplitData data;
    data.inputFile = inputFile;
    data.inputBatchNum = 0;
    data.inputFileLock = PTHREAD_MUTEX_INITIALIZER;
    data.k = k;
    data.kmerCount = kmerCount;
    data.singleBufferKmerCount = (uint64_t)maxRamGB * OneGiga / nThreads / sizeof(uint64_t);
    data.usmer = usmer;
    xassert(data.singleBufferKmerCount > 0, "maxRamGB too small or nThreads too large!");
    data.totalTaskNum = (data.kmerCount - 1) / data.singleBufferKmerCount + 1;
    for (uint64_t i = 0; i < data.totalTaskNum; ++i) {
        data.tmpOmerFileNames.push_back(tmpPath + string_format("/%lu.o.mer.tmp", i));
        if (usmer) {data.tmpImerFileNames.push_back(tmpPath + string_format("/%lu.i.mer.tmp", i));}
    }

    // 执行多线程分类+排序，结果写入到临时文件中
    err_func_printf(__func__, "spliting and sorting... (total %lu tasks)\n", data.totalTaskNum);
    kt_for(nThreads, ktf_split, &data, data.totalTaskNum);
    err_func_printf(__func__, "done spliting and sorting\n");

    // 多路归并
    err_func_printf(__func__, "merging...\n");
    string fullOutputFileName;
    fullOutputFileName = outputFileName + ".o.smer";
    merge_core(data.tmpOmerFileNames, fullOutputFileName, true, maxRamGB); 
    err_func_printf(__func__, "done merging to %s\n", fullOutputFileName.c_str());
    if (usmer) {
        fullOutputFileName = outputFileName + ".i.smer";
        merge_core(data.tmpImerFileNames, fullOutputFileName, true, maxRamGB); 
        err_func_printf(__func__, "done merging to %s\n", fullOutputFileName.c_str());
    }

    // 删除临时文件
    err_func_printf(__func__, "deleting...\n");
    for (auto tmpFileName : data.tmpOmerFileNames) {
        if (remove(tmpFileName.c_str()) != 0) {
            err_func_printf(__func__, "error delete temp file [%s]\n", tmpFileName.c_str());
        }
    }
    if (usmer) {
        for (auto tmpFileName : data.tmpImerFileNames) {
            if (remove(tmpFileName.c_str()) != 0) {
                err_func_printf(__func__, "error delete temp file [%s]\n", tmpFileName.c_str());
            }
        }
    }
    err_func_printf(__func__, "done deleting temp files\n");

    return 0;
}

