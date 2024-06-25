#include <vector>
#include <string>
#include <sys/resource.h>
#include <unordered_map>
#include "unitig.hpp"
#include "utils.hpp"
#include "klib/kthread.hpp"

using namespace std;

extern unsigned char nst_nt4_table[256];

struct UnitigData {
    vector<string> idList;
    vector<string> seqList;
    uint32_t kp1;
    vector<string> outputList;

    unordered_map<uint64_t, bool> *imap;
    unordered_map<uint64_t, bool> *omap;
};

void unitig_one_seq(string seq, uint32_t kmerLength, unordered_map<uint64_t, bool> &imap, unordered_map<uint64_t, bool> &omap, vector<string> &ret_string, vector<uint64_t> &ret_start_pos) {
    uint32_t kp1 = kmerLength + 1;
    string unitig = "";
    uint64_t kp1mer = 0;
    unsigned char newBase;
    bool multiIn, multiOut;
    bool badMod;
    
    uint64_t kp1merUpdateMask = (~(((uint64_t)-1) << (kp1 << 1)));
    if (kp1 == 32) {kp1merUpdateMask = ~kp1merUpdateMask;}

    size_t i = 0;
    for (; i < seq.length(); ++i) {

        newBase = nst_nt4_table[(int)seq[i]];
        kp1mer = ((kp1mer << 2) | (newBase & 0x3)) & kp1merUpdateMask;

        if (i == 0) {badMod = (newBase >= 4) ? (true) : (false);}

        if (newBase >= 4) {
            if (badMod == false) {
                if (unitig.length() > 0) {
                    ret_string.push_back(unitig);
                    ret_start_pos.push_back(i - unitig.length());
                    unitig = seq[i];
                    badMod = true;
                }
            }
            else if (badMod == true) {
                unitig += seq[i];
            }
        }
        else if (newBase < 4) {
            if (badMod == true) {
                if (unitig.length() > 0) {
                    ret_string.push_back(unitig);
                    ret_start_pos.push_back(i - unitig.length());
                    unitig = "";
                    badMod = false;
                }
            }
            if (unitig.length() < kmerLength) {
                unitig += seq[i];
            }
            else if (unitig.length() >= kmerLength) {
                multiIn = imap.find(kp1mer) != imap.end();
                multiOut = omap.find(kp1mer) != omap.end();

                if (!multiIn && !multiOut) {
                    unitig += seq[i];
                }
                else{
                    ret_string.push_back(unitig);
                    ret_start_pos.push_back(i - unitig.length());
                    unitig = seq.substr(i - kmerLength + 1, kmerLength);
                }
            }
        }
    }

    if (unitig.length() > 0) {
        ret_string.push_back(unitig);
        ret_start_pos.push_back(i - unitig.length());
    }
}

void ktf_unitig(void* shared, long i, int tid) {
    UnitigData *data = (UnitigData *)shared;
    uint32_t kp1 = data->kp1;

    string fullSeqName = data->idList[i];

    vector<string> ret_string;
    vector<uint64_t> ret_start_pos;
    unitig_one_seq(data->seqList[i], kp1 - 1, *(data->imap), *(data->omap), ret_string, ret_start_pos);
    data->outputList[i] = "";
    for (size_t idx = 0; idx < ret_string.size(); ++idx) {
        data->outputList[i] += string_format(">%lu|%lu|%s\n%s\n\n", idx, ret_start_pos[idx], fullSeqName.c_str(), ret_string[idx].c_str());
    }
}

int unitig_core(const std::string &ikFileName, const std::string &okFileName, 
    const std::string &faFileName, const std::string &outputFileName, 
    uint32_t nThreads, int maxRamGB) {

    uint32_t kp1, tmpKp1;
    uint64_t ikmerNum, okmerNum;

    FILE *ikFile = xopen(ikFileName.c_str(), "rb");
    setvbuf(ikFile, NULL, _IOFBF, CommonFileBufSize);
    err_fread_noeof(&kp1, sizeof(uint32_t), 1, ikFile);
    err_fread_noeof(&ikmerNum, sizeof(uint64_t), 1, ikFile);

    FILE *okFile = xopen(okFileName.c_str(), "rb");
    setvbuf(okFile, NULL, _IOFBF, CommonFileBufSize);
    err_fread_noeof(&tmpKp1, sizeof(uint32_t), 1, okFile);
    err_fread_noeof(&okmerNum, sizeof(uint64_t), 1, okFile);

    if (kp1 != tmpKp1) {
        err_func_printf(__func__, "ikFile[%s](k=%u) not consistent with okFile[%s](k=%u)", ikFileName.c_str(), kp1, okFileName.c_str(), tmpKp1);
        return 1;
    }
    
    unordered_map<uint64_t, bool> imap;
    unordered_map<uint64_t, bool> omap;

    imap.reserve(ikmerNum);
    omap.reserve(okmerNum);

    for (uint64_t i = 0; i < ikmerNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, ikFile);
        imap[kp1mer] = false;
    }

    for (uint64_t i = 0; i < okmerNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, okFile);
        omap[kp1mer] = false;
    }

    // 打开输入输出文件
    gzFile inputFile = xzopen(faFileName.c_str(), "r");
    gzbuffer(inputFile, CommonFileBufSize);
    kseq_t *seqs = kseq_init(inputFile);

    string fullOutputFileName = outputFileName + ".unitig";
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);

    // 创建线程池
    void *pool = kt_forpool_init(nThreads);
    struct UnitigData data;
    data.kp1 = kp1;
    data.imap = &imap;
    data.omap = &omap;

    // 估算剩余内存
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    double maxRam = (double)(maxRamGB - 1) * OneGiga - (double)usage.ru_maxrss * OneKilo; // 最后减1G是因为人最长染色体0.25G，加上输出brc是0.5G，多给1G应该就不会超标了
    err_func_printf(__func__, "total ram: %lf GB, used %lf GB, io buffer remaining %lf GB\n", (double)maxRamGB, (double)usage.ru_maxrss / OneMega, maxRam / OneGiga);

    // 分批读入数据处理
    double tmpRam;
    long taskNum;
    while (1) {
        tmpRam = 0;
        taskNum = 0;
        data.idList.clear();
        data.seqList.clear();
        data.outputList.clear();
        while (kseq_read(seqs) > 0) {
            tmpRam += seqs->seq.l * 2;
            taskNum += 1;
            string fullSeqName = string(seqs->name.s);
            if (seqs->comment.l > 0) fullSeqName += " " + string(seqs->comment.s);
            data.idList.push_back(fullSeqName);
            data.seqList.push_back(string(seqs->seq.s));
            data.outputList.push_back("");
            if (tmpRam >= maxRam) {
                break;
            }
        }

        if (taskNum == 0) break;

        getrusage(RUSAGE_SELF, &usage);
        err_func_printf(__func__, "Peak memory usage: %lfGB\n", double(usage.ru_maxrss) / OneMega); 
        err_func_printf(__func__, "kt_forpool summit %ld tasks\n", taskNum); 
        kt_forpool(pool, ktf_unitig, &data, taskNum);

        for (size_t i = 0; i < data.outputList.size(); ++i) {
            err_fprintf(outputFile, data.outputList[i].c_str());
        }

    }
    while (kseq_read(seqs) > 0) {
        string fullSeqName = string(seqs->name.s);
        if (seqs->comment.l > 0) fullSeqName += " " + string(seqs->comment.s);

        vector<string> ret_string;
        vector<uint64_t> ret_start_pos;
        unitig_one_seq(string(seqs->seq.s), kp1 - 1, imap, omap, ret_string, ret_start_pos);
        for (size_t i = 0; i < ret_string.size(); ++i) {
            err_fprintf(outputFile, ">%lu|%lu|%s\n%s\n\n", i, ret_start_pos[i], fullSeqName.c_str(), ret_string[i].c_str());
        }
        ret_string.clear();
        ret_start_pos.clear();
    }

    kseq_destroy(seqs);
    err_gzclose(inputFile);
    err_fclose(outputFile);

    return 0;
}