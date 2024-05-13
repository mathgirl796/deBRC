#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "utils.hpp"
#include "klib/kthread.hpp"
#define is_single(state) (((state) == 1)||((state) == 2)||((state) == 4)||((state) == 8))

int restore_core(const std::string &kFileName, const std::string &brcFileName, const std::string &outputFileName, uint32_t nThreads) {
    FILE *kFile = xopen(kFileName.c_str(), "rb");
    setvbuf(kFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t k = 0;
    uint64_t kmerNum = 0;
    err_fread_noeof(&k, sizeof(uint32_t), 1, kFile);
    err_fread_noeof(&kmerNum, sizeof(uint64_t), 1, kFile);
    map<uint64_t, uint8_t> km1merMap; // uint8_t的高四位无意义，低四位由低到高表示ACGT
    for (uint64_t i = 0; i < kmerNum; ++i) {
        uint64_t kmer;
        err_fread_noeof(&kmer, sizeof(uint64_t), 1, kFile);
        km1merMap[kmer >> 2] |= (0x1 << (kmer & 0x3));
    }
    err_fclose(kFile);

    for (auto it : km1merMap) {
        fprintf(stderr, "%s\t%d\t", uint64_to_kmer(it.first, 32).c_str(), it.second);
        if ((it.second & (uint8_t)0x1) != 0) fprintf(stderr, "%c\t", 'A');
        if ((it.second & (uint8_t)0x2) != 0) fprintf(stderr, "%c\t", 'C');
        if ((it.second & (uint8_t)0x4) != 0) fprintf(stderr, "%c\t", 'G');
        if ((it.second & (uint8_t)0x8) != 0) fprintf(stderr, "%c\t", 'T');
        fprintf(stderr, "\n");
    }

    gzFile brcFile = xzopen(brcFileName.c_str(), "r");
    kseq_t *kseq = kseq_init(brcFile);

    string fullOutputFileName = outputFileName + ".restore.fa";
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "w");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);

    extern unsigned char nst_nt4_table[256];
    char single2char[9] = {0, 'A', 'C', 0, 'G', 0, 0, 0, 'T'};

    uint64_t brcSeqNum = 0;
    while (kseq_read(kseq) > 0) {
        brcSeqNum ++;
        string brc(kseq->seq.s);
        // uint64_t brcLength = kseq->seq.l;
        string name(kseq->name.s);
        vector<string> info = split(name, '|');
        if (info.size() < 5) {
            err_func_printf(__func__, "brc file format error\n");
            err_fclose(outputFile);
            err_gzclose(brcFile);
            return 1;
        }
        for (size_t i = 0; i < info.size(); ++i) err_fprintf(stderr, "%s\t", info[i].c_str());
        err_fprintf(stderr, "\n");
        uint64_t brcId = std::stoull(info[0]);
        uint64_t seqLength = std::stoull(info[1]);
        // uint64_t startPosOnSeq = std::stoull(info[2]);
        string type = info[3];
        string id = info[4];

        if (brcId == 0) { // brcId为0，写入seq名字和头k-1个base
            if (brcSeqNum > 1) err_fprintf(outputFile, "\n\n");
            err_fprintf(outputFile, ">%s\n%s", id.c_str(), brc.substr(0, k - 1).c_str());
        }

        if (type == "bad") {
            err_fprintf(outputFile, "%s", brc.substr(k-1).c_str());
        }
        else if (type == "good") {
            uint64_t tmpSeqLength = k - 1;
            uint64_t tmpBrcPos = k - 1;
            uint64_t km1mer = 0;
            for (uint32_t j = 0; j < k - 1; ++j) { // 头一个kmer的前k-1个base
                km1mer = (km1mer << 2) | nst_nt4_table[(int)brc[j]];
            }
            err_fprintf(stderr, "%s\n", uint64_to_kmer(km1mer, 32).c_str());
            while (tmpSeqLength < seqLength) {
                uint8_t sucState = km1merMap[km1mer];
                unsigned char newBase; // 0, 1, 2, 3
                if (is_single(sucState)) {
                    err_fprintf(outputFile, "%c", single2char[sucState]);
                    newBase = nst_nt4_table[(int)single2char[sucState]];
                }
                else {
                    err_fprintf(outputFile, "%c", brc[tmpBrcPos]);
                    newBase = nst_nt4_table[(int)brc[tmpBrcPos]];
                    tmpBrcPos ++;
                }
                km1mer = ((km1mer << 2) | (newBase & 0x3)) & (~(((uint64_t)-1) << ((k - 1) << 1)));
                err_fprintf(stderr, "%s\n", uint64_to_kmer(km1mer, 32).c_str());
                tmpSeqLength ++;
            }
        }
    }

    err_fprintf(outputFile, "\n\n");
    err_fclose(outputFile);
    err_gzclose(brcFile);
    return 0;
}