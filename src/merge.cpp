#include <vector>
#include <string>
#include <queue>
#include <parallel/losertree.h>
#include "utils.hpp"

using namespace std;
/**
 * 给定一系列文件，每个文件格式如下：
 * uint32_t kmerLength | uint64_t kmerCount | uint64_t kmer (kmerCount个，从小到大有序排列)
 * 
 * 首先判断所有文件的kmerLength是否一致，如果不一致，返回1；
 * 一致，则写入kmerLength和一个占位的kmerCount
 * 然后，使用堆排序方法，将输入文件归并到输出文件中，归并时要去掉重复值（必要时要修改kmerCount）；
 * 最后，写入正确的kmerLength、kmerCount
 * 
 * 注：
 * 打开、关闭文件，读写文件，可以使用utils.hpp中的辅助函数
 * 打印调试信息可以用err_func_printf
 * 打开文件后应用setvbuf调整文件缓冲区大小，CommonFileBufSize即可(在utils.hpp中定义)
*/
int merge_core(const std::vector<std::string> &inputFiles, const std::string &outputFileName, const bool &distinct, uint32_t maxRamGB) {
    // 初始化优先级队列
    __gnu_parallel::_LoserTree<false, uint64_t, std::less<>> minHeap(inputFiles.size(), std::less<>());
    uint64_t minHeapBuf[inputFiles.size()];

    // 打开所有输入文件
    FILE* files[inputFiles.size()];
    uint32_t kmerLength = 0;
    uint64_t kmerCount = 0;
    for (size_t i = 0; i < inputFiles.size(); ++i) {
        uint32_t tmpLength;
        uint64_t tmpCount;
        files[i] = xopen(inputFiles[i].c_str(), "rb");
        setvbuf(files[i], NULL, _IOFBF, maxRamGB * OneGiga / inputFiles.size());
        // 读取头部
        err_fread_noeof(&tmpLength, sizeof(uint32_t), 1, files[i]);
        err_fread_noeof(&tmpCount, sizeof(uint64_t), 1, files[i]);
        if (kmerLength == 0) {
            kmerLength = tmpLength;
        }
        else if (kmerLength != tmpLength) {
            err_func_printf(__func__, "error: kmerLengths are not consistent!\n");
            return 1;
        }
        kmerCount += tmpCount;
        // 初始化堆
        uint64_t kmer;
        err_fread_noeof(&kmer, sizeof(uint64_t), 1, files[i]);
        minHeap.__insert_start(kmer, i, false);
        minHeapBuf[i] = kmer;
    }

    minHeap.__init();

    // 打开输出文件
    string fullOutputFileName = outputFileName + ".smer";
    FILE* outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    // 暂时写入头部
    err_fwrite(&kmerLength, sizeof(uint32_t), 1, outputFile);
    err_fwrite(&kmerCount, sizeof(uint64_t), 1, outputFile);

    // 使用败者树进行归并
    uint64_t closedFileNum = 0;
    while (closedFileNum < inputFiles.size()) {
        // 弹出堆顶元素（最小值）
        int fileIndex = minHeap.__get_min_source();
        uint64_t minValue = minHeapBuf[fileIndex];

        // 写入最小值到输出文件
        err_fwrite(&minValue, sizeof(uint64_t), 1, outputFile);

        // 从对应文件中读取下一个元素
        uint64_t kmer;
        size_t ret = fread(&kmer, sizeof(uint64_t), 1, files[fileIndex]);
        if (ret == 1) {
            minHeap.__delete_min_insert(kmer, false);
            minHeapBuf[fileIndex] = kmer;
        } else {
            minHeap.__delete_min_insert(UINT64_MAX, false);
            err_func_printf(__func__, "input file [%s] run out\n", inputFiles[fileIndex].c_str());
            err_fclose(files[fileIndex]);
            closedFileNum += 1;
        }
    }

    // 关闭输出文件
    err_fclose(outputFile);

    return 0;
}