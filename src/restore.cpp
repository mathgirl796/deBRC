#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <unordered_map>
#include "utils.hpp"
#include "klib/kthread.hpp"

static bool is_single(uint8_t state) {return (((state) == 1)||((state) == 2)||((state) == 4)||((state) == 8));}
extern unsigned char nst_nt4_table[256];
static char single2char[9] = {0, 'A', 'C', 0, 'G', 0, 0, 0, 'T'};

struct RestoreData {
    vector<string> infoList;
    vector<string> brcList;
    vector<string> outputList;
    unordered_map<uint64_t, uint8_t> kmerMap;
    uint32_t kp1;
    bool useKmerFormat;
};

void ktf_restore(void* data, long i, int tid) {
    RestoreData *d = (RestoreData *)data;
    // 解析brc中一条seq的头部
    string &brc = d->brcList[i];
    uint64_t brcLength = brc.length();
    bool useKmerFormat = d->useKmerFormat;
    vector<string> info = split(d->infoList[i], '|');
    if (info.size() < 5) {
        err_func_printf(__func__, "brc file format error\n");
    }
    uint64_t brcId = std::stoull(info[0]);
    uint64_t seqLength = std::stoull(info[1]);
    string type = info[3];
    string fullSeqName = info[4];
    uint32_t kp1 = d->kp1;
    uint32_t k = kp1 - 1;
    unordered_map<uint64_t, uint8_t> &kmerMap = d->kmerMap;
    err_func_printf(__func__, "worker:%d, task %ld start, processing %s, brcId:%lu, type:%s\n", tid, i, fullSeqName.c_str(), brcId, type.c_str());

    string output = "";

    if (brcId == 0 && type == "good") {output += string_format("%s", brc.substr(0, k).c_str());}

    // 读入seq开始恢复fasta seq
    uint64_t tmpSeqLength = 0;
    uint64_t tmpBrcPos = 0;
    if (type == "bad") {
        output += string_format("%s", brc.substr(0, brc.length()).c_str());
        tmpSeqLength = seqLength;
        tmpBrcPos = brcLength;
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
                exit(1);
            }
            unsigned char newBase; // 0, 1, 2, 3
            if (is_single(sucState)) {
                output += string_format("%c", single2char[sucState]);
                newBase = nst_nt4_table[(int)single2char[sucState]];
                if (newBase >= 4) {
                    err_func_printf(__func__, "newBase >= 4, there must be some bugs, error key: %s, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
                        uint64_to_kmer(kmer, k).c_str(), tmpSeqLength, seqLength, tmpBrcPos, brcLength);
                    exit(1);
                }
            }
            else {
                if (tmpBrcPos >= brcLength) {
                    err_func_printf(__func__, "tmpBrcPos >= brcLength, there must be some bugs, error key: %s, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
                        uint64_to_kmer(kmer, k).c_str(), tmpSeqLength, seqLength, tmpBrcPos, brcLength);
                    exit(1);
                }
                if (useKmerFormat) {tmpBrcPos += k - 1;}
                output += string_format("%c", brc[tmpBrcPos]);
                newBase = nst_nt4_table[(int)brc[tmpBrcPos]];
                tmpBrcPos ++;
            }
            kmer = ((kmer << 2) | (newBase & (uint8_t)0x3)) & (~(((uint64_t)-1) << (k << 1)));
            tmpSeqLength ++;
        }
    }

    d->outputList[i] += output;

    err_func_printf(__func__, "worker:%d done, tmpSeqLength: %lu/%lu, tmpBrcPos:%lu/%lu\n", 
        tid, tmpSeqLength, seqLength, tmpBrcPos, brcLength);
}

int restore_core(const std::string &smerFileName, const std::string &brcFileName, const std::string &outputFileName, uint32_t nThreads, bool useKmerFormat) {

    // 创建ktf worker所需数据结构（读入brcFile数据）
    err_func_printf(__func__, "loading brc file...\n");
    gzFile brcFile = xzopen(brcFileName.c_str(), "r");
    kseq_t *seqs = kseq_init(brcFile);
    struct RestoreData data;
    data.useKmerFormat = useKmerFormat;
    while (kseq_read(seqs) > 0) {
        string info(seqs->name.s);
        if (seqs->comment.l > 0) {info += " " + string(seqs->comment.s);}
        string brc(seqs->seq.s);
        data.infoList.push_back(info);
        data.brcList.push_back(brc);
        data.outputList.push_back("");
    }
    err_gzclose(brcFile);

    // 读入smer文件
    err_func_printf(__func__, "loading smer file...\n");
    FILE *smerFile = xopen(smerFileName.c_str(), "rb");
    setvbuf(smerFile, NULL, _IOFBF, CommonFileBufSize);
    uint64_t kp1merNum = 0;
    err_fread_noeof(&(data.kp1), sizeof(uint32_t), 1, smerFile);
    err_fread_noeof(&kp1merNum, sizeof(uint64_t), 1, smerFile);
    for (uint64_t i = 0; i < kp1merNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, smerFile);
        data.kmerMap[kp1mer >> 2] |= (((uint8_t)0x1) << (kp1mer & 0x3));
    }
    err_fclose(smerFile);
    err_func_printf(__func__, "done, total %lu kp1mers, kmerMap total %lu keys\n", kp1merNum, data.kmerMap.size());

    // ktf执行restore任务
    err_func_printf(__func__, "total %lu tasks\n", data.infoList.size());
    kt_for(nThreads, ktf_restore, &data, data.infoList.size());

    // 输出到文件
    err_func_printf(__func__, "output to file, total %lu tasks\n", data.infoList.size());
    string fullOutputFileName = outputFileName;
    if (useKmerFormat) {fullOutputFileName += "kmer";}
    fullOutputFileName += ".restore.fa";
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "w");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    for (size_t i = 0; i < data.infoList.size(); ++i) {
        vector<string> info = split(data.infoList[i], '|');
        uint64_t brcId = std::stoull(info[0]);
        // uint64_t seqLength = std::stoull(info[1]);
        string type = info[3];
        string fullSeqName = info[4];
        if (brcId == 0) { // fasta seq 头部
            if (i != 0) {err_fprintf(outputFile, "\n\n");}
            err_fprintf(outputFile, ">%s\n", fullSeqName.c_str());
        }
        err_fprintf(outputFile, "%s", data.outputList[i].c_str());
    }
    err_fclose(outputFile);

    err_func_printf(__func__, "done\n");

    return 0;
}