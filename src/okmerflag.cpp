#include <vector>
#include <string>
#include <sys/resource.h>
#include <unordered_map>
#include <unordered_set>
#include "utils.hpp"
#include "klib/kthread.hpp"
#include "okmerflag.hpp"
#include <string>

extern unsigned char nst_nt4_table[256];

struct OkmerFlagData {
    vector<string> seqList;
    vector<string> idList;
    vector<string> outputList;
    unordered_set<uint64_t> okmerSet;
    uint32_t kp1;
};

void ktf_okerflag(void *shared, long i, int tid) {
    OkmerFlagData *data = (OkmerFlagData *)shared;
    err_func_printf(__func__, "worker %d start work, processing %s (job %ld)\n", tid, data->idList[i].c_str(), i);
    string &seq = data->seqList[i];
    string &output = data->outputList[i];
    uint32_t &kp1 = data->kp1;

    if (seq.length() < kp1) {
        data->outputList[i] += seq;
    }
    else {
        unsigned char newBase;
        uint64_t kp1merUpdateMask = (kp1 == 32)? (~(~(((uint64_t)-1) << (kp1 << 1)))) : (~(((uint64_t)-1) << (kp1 << 1)));
        int badBasePos = -1;
        uint64_t kp1mer;
        for (size_t i = 0; i < seq.length(); ++i) {
            newBase = nst_nt4_table[(int)seq[i]];
            if (newBase >= 4) {badBasePos = kp1 - 1;}
            else {badBasePos -= 1;};
            kp1mer = ((kp1mer << 2) | (newBase & 0x3)) & kp1merUpdateMask;
            // for (auto x : data->okmerSet) {
            //     err_printf("%s\t", uint64_to_str(x, kp1).c_str());
            // }
            // err_printf("%s\n", uint64_to_str(kp1mer, kp1).c_str());

            if (i >= kp1-1 && badBasePos < 0 && data->okmerSet.find(kp1mer) != data->okmerSet.end()){
                output += (string[]){"0","1","2","3"}[newBase];
                // err_printf("%s\n", output.c_str());
            }
            else {
                output += seq[i];
            }
        }
    }
    output += "\n";
    err_func_printf(__func__, "worker %d finish work\n", tid);
}

int okmerflag_core(std::string refFileName, std::string okFileName, std::string outputFileName, uint32_t nThreads){

    OkmerFlagData data;

    // 加载参考基因组到内存
    gzFile fp = xzopen(refFileName.c_str(), "r");
    gzbuffer(fp, CommonFileBufSize);
    kseq_t *seqs = kseq_init(fp);
    data.idList.clear();
    data.seqList.clear();
    while (kseq_read(seqs) > 0) {
        string fullSeqName = string(seqs->name.s);
        if (seqs->comment.l > 0) fullSeqName += " " + string(seqs->comment.s);
        data.idList.push_back(fullSeqName);
        data.seqList.push_back(string(seqs->seq.s));
        data.outputList.push_back(">"+fullSeqName+"\n");
    }

    // 加载okmer到内存
    FILE *omer = xopen(okFileName.c_str(), "rb");
    setvbuf(omer, NULL, _IOFBF, CommonFileBufSize);
    uint32_t kp1 = 0;
    uint64_t omerNum = 0;
    err_fread_noeof(&kp1, sizeof(uint32_t), 1, omer);
    err_fread_noeof(&omerNum, sizeof(uint64_t), 1, omer);
    err_func_printf(__func__, "loading %s\n", okFileName.c_str());
    data.kp1 = kp1;
    data.okmerSet.reserve(omerNum);
    for (uint64_t i = 0; i < omerNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, omer);
        data.okmerSet.insert(kp1mer);
    }
    err_fclose(omer);
    err_func_printf(__func__, "done, total %lu omers\n", data.okmerSet.size());

    // 进行标注
    kt_for(nThreads, ktf_okerflag, &data, data.idList.size());

    // 输出
    FILE *outputFile = xopen((outputFileName + ".okmerflag").c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    for (size_t i = 0; i < data.outputList.size(); ++i) {
        err_fprintf(outputFile, data.outputList[i].c_str());
    }
    err_fclose(outputFile);

    return 0;
}