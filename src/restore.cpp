#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "utils.hpp"
#include "klib/kthread.hpp"

static bool is_single(uint8_t state) {return (((state) == 1)||((state) == 2)||((state) == 4)||((state) == 8));}

int restore_core(const std::string &smerFileName, const std::string &brcFileName, const std::string &outputFileName, uint32_t nThreads) {
    FILE *smerFile = xopen(smerFileName.c_str(), "rb");
    setvbuf(smerFile, NULL, _IOFBF, CommonFileBufSize);
    uint32_t kp1 = 0;
    uint64_t kp1merNum = 0;
    err_fread_noeof(&kp1, sizeof(uint32_t), 1, smerFile);
    err_fread_noeof(&kp1merNum, sizeof(uint64_t), 1, smerFile);
    uint32_t k = kp1 - 1;
    map<uint64_t, uint8_t> kmerMap; // uint8_t的高四位无意义，低四位由低到高表示ACGT
    for (uint64_t i = 0; i < kp1merNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, smerFile);
        kmerMap[kp1mer >> 2] |= (((uint8_t)0x1) << (kp1mer & 0x3));
    }
    err_fclose(smerFile);

    // for (auto it : kmerMap) {
    //     if (!is_single(it.second)) {
    //         fprintf(stderr, "%s\t%d\t", uint64_to_kmer(it.first, k).c_str(), it.second);
    //         if ((it.second & (uint8_t)0x1) != 0) fprintf(stderr, "%c\t", 'A');
    //         if ((it.second & (uint8_t)0x2) != 0) fprintf(stderr, "%c\t", 'C');
    //         if ((it.second & (uint8_t)0x4) != 0) fprintf(stderr, "%c\t", 'G');
    //         if ((it.second & (uint8_t)0x8) != 0) fprintf(stderr, "%c\t", 'T');
    //         fprintf(stderr, "\n");
    //     }
    // }

    gzFile brcFile = xzopen(brcFileName.c_str(), "r");
    kseq_t *seqs = kseq_init(brcFile);

    string fullOutputFileName = outputFileName + ".restore.fa";
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "w");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);

    string output = "";
    extern unsigned char nst_nt4_table[256];
    char single2char[9] = {0, 'A', 'C', 0, 'G', 0, 0, 0, 'T'};

    uint64_t brcSeqNum = 0;
    while (kseq_read(seqs) > 0) {
        // 解析brc中一条seq的头部
        brcSeqNum ++;
        string brc(seqs->seq.s);
        uint64_t brcLength = seqs->seq.l;
        string name(seqs->name.s);
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
        string type = info[3];
        string fullSeqName = info[4];
        if (seqs->comment.l > 0) fullSeqName += " " + string(seqs->comment.s);

        if (brcId == 0) { // brcId为0，写入seq名字和头一个kmer
            if (brcSeqNum > 1) output += "\n\n"; // brcId为0，但不是第一个brcId为零的，说明不是第一条fasta序列，要加换行
            output += string_format(">%s\n", fullSeqName.c_str()); // fasta id
            if (type == "good") output += string_format("%s", brc.substr(0, k).c_str());
        }

        // 读入seq开始恢复fasta seq
        uint64_t tmpSeqLength = 0;
        uint64_t tmpBrcPos = 0;
        if (type == "bad") {
            output += string_format("%s", brc.substr(0, brc.length()).c_str());
        }
        else if (type == "good") {
            uint64_t kmer = 0;
            for (uint32_t i = 0; i < k; ++i) { // 头一个kmer
                kmer = (kmer << 2) | nst_nt4_table[(int)brc[i]];
            }
            tmpSeqLength = k;
            tmpBrcPos = k;
            uint8_t sucState = 0;
            while (tmpSeqLength < seqLength) {
                sucState = kmerMap[kmer];
                if (sucState == 0) {
                    err_func_printf(__func__, "sucState == 0, there must be some bugs, error key: %s, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
                        uint64_to_kmer(kmer, k).c_str(), tmpSeqLength, seqLength, tmpBrcPos, brcLength);
                }
                unsigned char newBase; // 0, 1, 2, 3
                if (is_single(sucState)) {
                    output += string_format("%c", single2char[sucState]);
                    newBase = nst_nt4_table[(int)single2char[sucState]];
                    if (newBase >= 4) {
                        err_func_printf(__func__, "newBase >= 4, there must be some bugs, error key: %s, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
                            uint64_to_kmer(kmer, k).c_str(), tmpSeqLength, seqLength, tmpBrcPos, brcLength);
                    }
                }
                else {
                    if (tmpBrcPos >= brcLength) {
                        err_func_printf(__func__, "tmpBrcPos >= brcLength, there must be some bugs, error key: %s, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
                            uint64_to_kmer(kmer, k).c_str(), tmpSeqLength, seqLength, tmpBrcPos, brcLength);
                    }
                    output += string_format("%c", brc[tmpBrcPos]);
                    newBase = nst_nt4_table[(int)brc[tmpBrcPos]];
                    tmpBrcPos ++;
                }
                kmer = ((kmer << 2) | (newBase & (uint8_t)0x3)) & (~(((uint64_t)-1) << (k << 1)));
                // err_fprintf(stderr, "%s\n", uint64_to_kmer(kmer, 32).c_str());
                tmpSeqLength ++;
            }
        }

        err_func_printf(__func__, "%s brcId:%lu done, type:%s, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
            fullSeqName.c_str(), brcId, type.c_str(), tmpSeqLength, seqLength, tmpBrcPos, brcLength);
    }

    err_fprintf(outputFile, "%s\n\n", output.data());
    err_fclose(outputFile);
    err_gzclose(brcFile);
    return 0;
}