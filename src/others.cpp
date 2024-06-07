#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <string>
#include "utils.hpp"
using namespace std;

struct CheckData{
    string inputFileName;
    int taskID;
    int taskNum;
    uint64_t taskQuan;
    uint64_t lastTaskQuan;
    bool check_dup;
};

void *multi_check(void *data) {
    CheckData *d = (CheckData *)data;
    // 打开文件
    FILE* inputFile = xopen(d->inputFileName.c_str(), "rb");
    setvbuf(inputFile, NULL, _IOFBF, CommonFileBufSize);
    // 跳转到任务开始处
    err_fseek(inputFile, sizeof(uint32_t) + sizeof(uint64_t) + d->taskID * d->taskQuan * sizeof(uint64_t), SEEK_SET);
    uint64_t kmerCount = d->taskID != d->taskNum - 1 ? d->taskQuan + 1 : d->lastTaskQuan;
    err_func_printf(__func__, "tid:%d, check from %lu to %lu\n", d->taskID, d->taskID * d->taskQuan, d->taskID * d->taskQuan + kmerCount - 1);
    // 读出第一个kmer
    uint64_t last_kmer;
    uint64_t kmer;
    fread(&last_kmer, sizeof(uint64_t), 1, inputFile);
    // 检查后续kmer
    for (uint64_t i = 1; i < kmerCount; ++i) {
        int ret = fread(&kmer, sizeof(uint64_t), 1, inputFile);
        if (ret != 1) {
            err_func_printf(__func__, "tid:%d, error: unexpected end of file at [%lu]\n", d->taskID, i + d->taskID * d->taskQuan);
            exit(1);
        }
        if (kmer < last_kmer) {
            err_func_printf(__func__, "tid:%d, error: decrease at [%lu]\n", d->taskID, i + d->taskID * d->taskQuan);
            exit(1);
        }
        if (d->check_dup && kmer == last_kmer) {
            err_func_printf(__func__, "tid:%d, error: duplicate at [%lu]\n", d->taskID, i + d->taskID * d->taskQuan);
            exit(1);
        }
        last_kmer = kmer;
    }
    err_func_printf(__func__, "tid:%d, passed check from %lu to %lu\n", d->taskID, d->taskID * d->taskQuan, d->taskID * d->taskQuan + kmerCount);
    err_fclose(inputFile);
    return NULL;
}

int check_core(const std::string &inputFileName, uint32_t nThreads, bool check_dup) {
    FILE* inputFile = xopen(inputFileName.c_str(), "rb");
    uint32_t kmerLength;
    uint64_t kmerCount;
    // 读出头部
    err_fread_noeof(&kmerLength, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kmerCount, sizeof(uint64_t), 1, inputFile);
    err_func_printf(__func__, "kmerLength:%u, kmerCount:%lu\n", kmerLength, kmerCount);
    err_fclose(inputFile);
    // 分任务
    int taskNum = kmerCount > nThreads ? nThreads : 1;
    uint64_t taskQuan = kmerCount / taskNum;
    err_func_printf(__func__, "check by %d threads\n", taskNum);
    // 多线程检查
    CheckData data_blocks[taskNum];
    pthread_t tid[taskNum];
    for (uint32_t i = 0; i < nThreads; ++i) {
        data_blocks[i].inputFileName = inputFileName;
        data_blocks[i].taskID = i;
        data_blocks[i].taskNum = taskNum;
        data_blocks[i].taskQuan = taskQuan;
        data_blocks[i].lastTaskQuan = kmerCount - (taskNum - 1) * taskQuan;
        data_blocks[i].check_dup = check_dup;
        pthread_create(&tid[i], NULL, multi_check, &data_blocks[i]);
    }
    for (uint32_t i = 0; i < nThreads; ++i) {
        pthread_join(tid[i], NULL);
    }
    err_func_printf(__func__, "all pass (check_dup: %d)\n", (int)check_dup);
    return 0;
}



int view_core(const std::string &inputFileName, bool head) {
    FILE* inputFile = xopen(inputFileName.c_str(), "r");
    setvbuf(inputFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t kmerLength;
    uint64_t kmerCount;
    err_fread_noeof(&kmerLength, sizeof(uint32_t), 1, inputFile);
    err_fread_noeof(&kmerCount, sizeof(uint64_t), 1, inputFile);
    printf("kmerLength\tkmerCount\n");
    printf("%u\t%lu\n", kmerLength, kmerCount);
    if (head) return 0;
    uint64_t kmer;
    for (uint64_t i = 0; i < kmerCount; ++i) {
        if (i % 5 == 0 && i != 0) {std::cerr << std::endl; }
        if (i % 20 == 0 && i != 0) {getchar();}
        err_fread_noeof(&kmer, sizeof(uint64_t), 1, inputFile);
        std::cerr << uint64_to_str(kmer, kmerLength) << "\t";
    }
    std::cerr << std::endl;
    return 0;
}

static int compare_uint64 (const void *ptr_a, const void *ptr_b) {
    uint64_t a = *(uint64_t *)ptr_a, b = *(uint64_t *)ptr_b;
    if (a < b) return -1;
    else if (a == b) return 0;
    else return 1;
} 

int search_core(const std::string &smerFileName, std::string kmerStr) {
    FILE* smerFile = xopen(smerFileName.c_str(), "rb");
    uint32_t smerK;
    uint64_t kmerCount;
    // 读出头部
    err_fread_noeof(&smerK, sizeof(uint32_t), 1, smerFile);
    err_fread_noeof(&kmerCount, sizeof(uint64_t), 1, smerFile);

    // 检测输入数据是否合法（只能是kmer或km1mer或0）
    uint32_t inputK = (uint32_t)kmerStr.length();
    bool interactive = false;
    if (inputK == 0) interactive = true;
    else if (inputK != smerK && inputK != smerK - 1) {
        err_fclose(smerFile);
        err_fprintf(stderr, "inputK must equal smerK or smerKm1!\n");
        return 1;
    }

    // 读出smerList数据
    uint64_t* kmerList = (uint64_t *)err_calloc(__func__, kmerCount, sizeof(uint64_t));
    err_fread_noeof(kmerList, sizeof(uint64_t), kmerCount, smerFile);

    for(;;) {
        vector<string> workList;
        if (interactive) {
            err_fprintf(stderr, "smerK=%u, input a kmer or km1mer (enter to exit): ", smerK);
            getline(cin, kmerStr);
        }

        inputK = kmerStr.length();
        if (inputK == 0 || kmerStr == "\r") {
            return 0;
        }
        else if (inputK == smerK) {
            workList.push_back(kmerStr);
        }
        else if (inputK == smerK - 1) {
            workList = {"A"+kmerStr,"C"+kmerStr,"G"+kmerStr,"T"+kmerStr,kmerStr+"A",kmerStr+"C",kmerStr+"G",kmerStr+"T"};
        }
        else {
            err_fprintf(stderr, "wrong inputK, inputK must equal smerK or smerKm1!\n");
            if (interactive) continue;
            else break;
        }
        for (auto tmpKmerStr: workList) {
            inputK = tmpKmerStr.length();
            uint64_t kmer = str_to_uint64(tmpKmerStr, inputK);
            void* retPtr = bsearch(&kmer, kmerList, kmerCount, sizeof(uint64_t), compare_uint64);
            err_printf("%s\t%s\n", uint64_to_str(kmer, inputK).c_str(), retPtr ? string_format("%lu", ((uint64_t)retPtr - (uint64_t)kmerList) / sizeof(uint64_t)).c_str() : "not exist");
        }
    }

    free(kmerList);
    err_fclose(smerFile);
    return 0;
}