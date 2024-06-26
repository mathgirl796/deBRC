#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <parallel/losertree.h>
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
    uint64_t kmerLength;
    uint64_t kmerCount;
    uint64_t singleBufferKmerCount;
    uint64_t totalTaskNum;
    vector<string> tmpOmerFileNames;

    uint64_t L2[5] = {UINT64_MAX};
    pthread_mutex_t L2Lock;
};

void ktf_split(void* data, long i, int tid) {
    // 提取要处理的数据范围
    SplitData *d = (SplitData *)data;
    pthread_mutex_lock(&d->inputFileLock);
    i = (long)d->inputBatchNum;
    uint64_t startCount = i * d->singleBufferKmerCount;
    uint64_t overlapCount = 4; // 分任务的起始位置应该适当提前（smer形式需要提前4个，usmer形式要提前16个）
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
    FILE *omerFile;
    uint64_t omerCount = 0;

    kp1 = d->kmerLength;

    omerFile = xopen(d->tmpOmerFileNames[i].c_str(), "wb+");
    setvbuf(omerFile, NULL, _IOFBF, CommonFileBufSize);
    err_fwrite(&kp1, sizeof(uint32_t), 1, omerFile);
    err_fwrite(&omerCount, sizeof(uint64_t), 1, omerFile);

    // 分析读入的kmer，写入到临时文件
    uint64_t lastLkmer = 0;
    uint64_t lastLbase = 0;
    for (uint64_t i = 0; i < bufferCount; ++i) {
        uint64_t kp1mer = kpmers[i];
        uint64_t lkmer = kp1mer >> 2;
        uint64_t lbase = kp1mer >> ((kp1 - 1) << 1);
        if (i > 0 && lkmer == lastLkmer) {
            err_fwrite(&(kpmers[i-1]), sizeof(uint64_t), 1, omerFile);
            err_fwrite(&(kpmers[i]), sizeof(uint64_t), 1, omerFile);
            omerCount += 2;
        }
        if (i > 0 && lbase != lastLbase) {
            pthread_mutex_lock(&d->L2Lock);
            d->L2[lbase] = startCount + i;
            pthread_mutex_unlock(&d->L2Lock);
        }
        lastLkmer = lkmer;
        lastLbase = lbase;
    }

    free(kpmers);

    err_fseek(omerFile, sizeof(uint32_t), SEEK_SET);
    err_fwrite(&omerCount, sizeof(uint64_t), 1, omerFile); // 写入实际mer个数

    // 读出刚刚提取的mer数据进行排序
    err_fseek(omerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); // 跳过头部读出数据
    uint64_t *kmers = (uint64_t *)err_calloc(__func__, omerCount, sizeof(uint64_t));
    err_fread_noeof(kmers, sizeof(uint64_t), omerCount, omerFile);
    qsort(kmers, omerCount, sizeof(uint64_t), compare_uint64);
    err_fseek(omerFile, sizeof(uint32_t) + sizeof(uint64_t), SEEK_SET); // 跳过头部写入数据
    err_fwrite(kmers, sizeof(uint64_t), omerCount, omerFile);
    free(kmers);

    // 关闭文件
    err_fclose(omerFile);

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
int split_core(const std::string &inputFileName, const std::string &outputFileName, const bool ismer, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads) {
    
    FILE* inputFile = xopen(inputFileName.c_str(), "r");
    setvbuf(inputFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t kmerLength;
    uint64_t kmerCount;
    err_fread_noeof(&kmerLength, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kmerCount, sizeof(uint64_t), 1, inputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", kmerLength, kmerCount);

    SplitData data;
    data.inputFile = inputFile;
    data.inputBatchNum = 0;
    data.inputFileLock = PTHREAD_MUTEX_INITIALIZER;
    data.L2Lock = PTHREAD_MUTEX_INITIALIZER;
    data.kmerLength = kmerLength;
    data.kmerCount = kmerCount;
    data.singleBufferKmerCount = (uint64_t)maxRamGB * OneGiga / nThreads / sizeof(uint64_t);
    xassert(data.singleBufferKmerCount > 0, "maxRamGB too small or nThreads too large!");
    data.totalTaskNum = (data.kmerCount - 1) / data.singleBufferKmerCount + 1;
    for (uint64_t i = 0; i < data.totalTaskNum; ++i) {
        data.tmpOmerFileNames.push_back(tmpPath + string_format("/%lu.o.mer.tmp", i));
    }

    // 执行多线程分类+排序，结果写入到临时文件中
    err_func_printf(__func__, "spliting and sorting... (total %lu tasks)\n", data.totalTaskNum);
    kt_for(nThreads, ktf_split, &data, data.totalTaskNum);
    fclose(inputFile);
    err_func_printf(__func__, "done spliting and sorting\n");

    // 多路归并
    err_func_printf(__func__, "merging...\n");
    string fullOutputFileName;
    fullOutputFileName = outputFileName + ".o.smer";
    merge_core(data.tmpOmerFileNames, fullOutputFileName, true, maxRamGB); 
    err_func_printf(__func__, "done merging to %s\n", fullOutputFileName.c_str());

    // 删除临时文件
    err_func_printf(__func__, "deleting...\n");
    for (auto tmpFileName : data.tmpOmerFileNames) {
        if (remove(tmpFileName.c_str()) != 0) {
            err_func_printf(__func__, "error delete temp file [%s]\n", tmpFileName.c_str());
        }
    }
    err_func_printf(__func__, "done deleting temp files\n");

    // 多入kmer
    if (ismer) {
        uint64_t *L2 = data.L2;
        L2[0] = 0, L2[4] = kmerCount;
        for (int i = 3; i >= 1; --i) {if (L2[i] == UINT64_MAX) {L2[i] = L2[i +1];}} // 未经测试的，出了bug再说，欸嘿
        err_func_printf(__func__, "%lu\t%lu\t%lu\t%lu\t%lu\n", L2[0], L2[1], L2[2], L2[3], L2[4]);
        FILE* f[4];
        for (int i = 0; i < 4; ++i) {
            f[i] = xopen(inputFileName.c_str(), "rb");
            setvbuf(f[i], NULL, _IOFBF, CommonFileBufSize);
            err_fseek(f[i], sizeof(uint32_t) + sizeof(uint64_t) * (L2[i] + 1), SEEK_CUR);
        }


        // 四路归并
        // 初始化优先级队列
        __gnu_parallel::_LoserTree<false, uint64_t, std::less<>> minHeap(4, std::less<>());
        uint64_t minHeapBuf[4];
        // 初始化堆
        uint64_t LP[4] = {L2[0], L2[1], L2[2], L2[3]}; // 四个分区游标位置
        for (int i = 0; i < 4; ++i) {
            if (L2[i+1] - L2[i] > 0) {
                uint64_t kmer;
                err_fread_noeof(&kmer, sizeof(uint64_t), 1, f[i]);
                LP[i] += 1;
                kmer = ror_kmer(kmer, kmerLength, -1); // 把首个base甩到最后面
                minHeap.__insert_start(kmer, i, false);
                minHeapBuf[i] = kmer;
                // err_func_printf(__func__, "%s\n", uint64_to_str(kmer, kmerLength).c_str());
            }
        }
        minHeap.__init();
        // 打开输出文件
        fullOutputFileName = outputFileName + ".i.smer";
        FILE* outputFile = xopen(fullOutputFileName.c_str(), "wb");
        setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
        err_fwrite(&kmerLength, sizeof(uint32_t), 1, outputFile); // 写入kmer长度
        err_fseek(outputFile, sizeof(uint64_t), SEEK_CUR); // kmer数量留空

        // 使用败者树进行归并
        uint64_t closedFileNum = 0;
        bool isFirstKmer = true;
        uint64_t lastKmer = 0;
        uint64_t lastWriteKmer = UINT64_MAX;
        uint64_t imerCount = 0;
        while (closedFileNum < 4) {
            // 弹出堆顶元素（最小值）
            int idx = minHeap.__get_min_source();
            uint64_t minValue = minHeapBuf[idx];
            // err_func_printf(__func__, "%s\t%s\t%s\t%s\n" 
            //     , uint64_to_str(minHeapBuf[0], kmerLength).c_str()
            //     , uint64_to_str(minHeapBuf[1], kmerLength).c_str()
            //     , uint64_to_str(minHeapBuf[2], kmerLength).c_str()
            //     , uint64_to_str(minHeapBuf[3], kmerLength).c_str()
            //     );
            

            // 从对应文件中读取下一个元素到败者树中
            if (LP[idx] >= L2[idx+1]) {
                minHeap.__delete_min_insert(UINT64_MAX, false);
                err_func_printf(__func__, "lbase part [%c] run out\n", "ACGT"[idx]);
                err_fclose(f[idx]);
                closedFileNum += 1;
            }
            else {
                uint64_t kmer;
                err_fread_noeof(&kmer, sizeof(uint64_t), 1, f[idx]);
                LP[idx] += 1;
                kmer = ror_kmer(kmer, kmerLength, -1);
                minHeap.__delete_min_insert(kmer, false);
                minHeapBuf[idx] = kmer;
            }

            // 处理split逻辑
            if (isFirstKmer) {
                isFirstKmer = false;
            }
            else if ((lastKmer >> 2) == (minValue >> 2)){
                uint64_t tmpKmer;
                if (lastKmer != lastWriteKmer) {
                    tmpKmer = ror_kmer(lastKmer, kmerLength, 1);
                    err_fwrite(&tmpKmer, sizeof(uint64_t), 1, outputFile);
                    imerCount += 1;
                }
                tmpKmer = ror_kmer(minValue, kmerLength, 1);
                err_fwrite(&tmpKmer, sizeof(uint64_t), 1, outputFile);
                imerCount += 1;
                lastWriteKmer = minValue;
            }

            // 写入最小值到输出文件
            lastKmer = minValue;
        }
        // 写入kmer数量
        err_fseek(outputFile, sizeof(uint32_t), SEEK_SET); 
        err_fwrite(&imerCount, sizeof(uint64_t), 1, outputFile); // 写入kmer数量
        // 关闭输出文件
        err_fclose(outputFile);
    }

    return 0;
}

