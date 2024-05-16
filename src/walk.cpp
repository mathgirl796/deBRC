#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include "utils.hpp"
#include "klib/kthread.hpp"

extern unsigned char nst_nt4_table[256];

enum class MerType {x, i, o};
enum class BrcType {good, bad};

struct WalkData {
    vector<string> seqList;
    vector<string> idList;
    unordered_set<uint64_t> okmerSet;
    uint32_t kp1;
    bool passSpecialCharactors;

    // fastaFormat output
    vector<string> outputList;

    // kmerFormat output
    bool useKmerFormat;
};

void ktf_walk(void* data, long i, int tid) {
    WalkData *walkData = (WalkData *)data;
    // 处理输入中的一个seq
    err_func_printf(__func__, "worker %d start work, processing %s (job %ld)\n", tid, walkData->idList[i].c_str(), i); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    string id = walkData->idList[i];
    string seq = walkData->seqList[i];
    string output = "";
    uint32_t kp1 = walkData->kp1;
    uint32_t k = kp1 - 1;
    uint64_t brc_id = 0;
    uint64_t original_seq_length = 0;
    uint64_t kp1mer = 0;
    uint64_t startPosOnSeq = 0;
    BrcType brcType;
    MerType merType;
    int badBasePos = -1;
    string brc;
    unsigned char newBase;
    uint64_t kp1merUpdateMask = (~(((uint64_t)-1) << (kp1 << 1)));
    if (kp1 == 32) {kp1merUpdateMask = ~kp1merUpdateMask;}
    // 遍历所有kp1mer
    for (size_t i = 0; i < seq.length() - kp1 + 1; ++i) {
        /* 一条seq头一个kp1mer的处理！（初始化读入kmer） */
        if (i == 0) {
            for (uint32_t j = 0; j < kp1 - 1; ++j) { // 读入头一个kmer的前k-1个base
                newBase = nst_nt4_table[(int)seq[j]];
                if (newBase >= 4) { // 处理特殊base
                    badBasePos = kp1 - 1;
                }
                else {
                    badBasePos --;
                }
                kp1mer = ((kp1mer << 2) | (newBase & 0x3)) & kp1merUpdateMask;
            }
            brc = seq.substr(0, k); // 一条sequence的头一个kmer一定要写入brc
            original_seq_length = k;
            startPosOnSeq = 0;
        }

        // 读入一个base
        newBase = nst_nt4_table[(int)seq[i + kp1 - 1]];
        // 处理特殊base
        if (newBase >= 4) {
            badBasePos = kp1 - 1;
        }
        else {
            badBasePos --;
        }
        // 把读入的base加入kp1mer
        kp1mer = ((kp1mer << 2) | (newBase & 0x3)) & kp1merUpdateMask;
        if (badBasePos >= 0) merType = MerType::x;
        else if (walkData->okmerSet.find(kp1mer) == walkData->okmerSet.end()) merType = MerType::i;
        else merType = MerType::o;

        /* 一条seq头一个kp1mer的处理！（初始化BrcType） */
        if (i == 0) {
            brcType = (merType == MerType::x) ? (BrcType::bad) : (BrcType::good);
        }

        // err_printf("%d\t %d\n", brcType, merType);

        if (brcType == BrcType::good && merType == MerType::o) {
            if (walkData->useKmerFormat) brc += seq.substr(i, k + 1);
            else brc += seq[i + kp1 - 1];
            original_seq_length ++;
        }
        else if (brcType == BrcType::good && merType == MerType::i) {
            original_seq_length ++;
        }
        else if (brcType == BrcType::good && merType == MerType::x) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            output += string_format(">%lu|%lu|%lu|%s|%s\n%s\n\n", brc_id, original_seq_length, startPosOnSeq, (brcType == BrcType::good) ? ("good"):("bad"), id.c_str(), brc.c_str());
            // 确定下一个brc的类型
            brcType = BrcType::bad;
            // 初始化下一个brc到内存中
            brc_id ++;
            brc = seq[i + kp1 - 1];
            original_seq_length = 1;
            startPosOnSeq = i + kp1 - 1;
        }
        else if (brcType == BrcType::bad && merType == MerType::i) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            if (walkData->passSpecialCharactors != true)
                output += string_format(">%lu|%lu|%lu|%s|%s\n%s\n\n", brc_id, original_seq_length - k, startPosOnSeq, (brcType == BrcType::good) ? ("good"):("bad"), id.c_str(), brc.substr(0, brc.length() - k).c_str());
            // 确定下一个brc的类型
            brcType = BrcType::good;
            // 初始化下一个brc到内存中
            brc_id ++;
            brc = seq.substr(i, k);
            original_seq_length = kp1;
            startPosOnSeq = i;
        }
        else if (brcType == BrcType::bad && merType == MerType::o) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            if (walkData->passSpecialCharactors != true)
                output += string_format(">%lu|%lu|%lu|%s|%s\n%s\n\n", brc_id, original_seq_length - k, startPosOnSeq, (brcType == BrcType::good) ? ("good"):("bad"), id.c_str(), brc.substr(0, brc.length() - k).c_str());
            // 确定下一个brc的类型
            brcType = BrcType::good;
            // 初始化下一个brc到内存中
            brc_id ++;
            brc = seq.substr(i, k);
            if (walkData->useKmerFormat) brc += seq.substr(i, k + 1);
            else brc += seq[i + kp1 - 1];
            original_seq_length = kp1;
            startPosOnSeq = i;
        }
        else if (brcType == BrcType::bad && merType == MerType::x) {
            brc += seq[i + kp1 - 1];
            original_seq_length ++;
        }
    }
    // 输出最后一个brc到文件
    if (brcType == BrcType::good || (brcType == BrcType::bad && walkData->passSpecialCharactors != true)) {
        output += string_format(">%lu|%lu|%lu|%s|%s\n%s\n\n", brc_id, original_seq_length, startPosOnSeq, (brcType == BrcType::good) ? ("good"):("bad"), id.c_str(), brc.c_str());
        walkData->outputList[i] = output;
    }
    err_func_printf(__func__, "worker %d finish work\n", tid);
}

int walk_core(const std::string &smerFileName, const std::string &okFileName, 
    const std::string &faFileName, const std::string &outputFileName, uint32_t nThreads,
    bool passSpecialCharactors, bool useKmerFormat) {

    uint32_t kp1 = 0;
    uint64_t omerNum = 0;

    FILE *omer = xopen(okFileName.c_str(), "rb");
    setvbuf(omer, NULL, _IOFBF, CommonFileBufSize);
    err_fread_noeof(&kp1, sizeof(uint32_t), 1, omer);
    err_fread_noeof(&omerNum, sizeof(uint64_t), 1, omer);

    string fullOutputFileName = outputFileName;
    if (passSpecialCharactors) fullOutputFileName += ".passN"; // 若不保留未知字符导致的xBrc，则加后缀注明
    if (useKmerFormat) fullOutputFileName += ".kmer"; 
    fullOutputFileName += ".brc";

    // 创建ktf worker所需数据结构
    struct WalkData walkData;
    walkData.kp1 = kp1;
    walkData.outputList = vector<string>();
    walkData.passSpecialCharactors = passSpecialCharactors;
    walkData.useKmerFormat = useKmerFormat;
    gzFile fp = xzopen(faFileName.c_str(), "r");
    kseq_t *seqs = kseq_init(fp);
    err_func_printf(__func__, "loading %s\n", faFileName.c_str()); 
    while (kseq_read(seqs) > 0) {
        // fprintf(stderr, "%s|%s\n", seqs->name.s, seqs->comment.s);
        walkData.seqList.push_back(string(seqs->seq.s));
        string fullSeqName = string(seqs->name.s);
        if (seqs->comment.l > 0) fullSeqName += " " + string(seqs->comment.s);
        walkData.idList.push_back(fullSeqName);
        walkData.outputList.push_back("");
    }
    gzclose(fp);

    // 加载okmer到内存
    err_func_printf(__func__, "loading %s\n", okFileName.c_str());
    for (uint64_t i = 0; i < omerNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, omer);
        walkData.okmerSet.insert(kp1mer);
    }
    err_fclose(omer);
    err_func_printf(__func__, "done, total %lu omers\n", walkData.okmerSet.size());


    // 执行walk任务
    err_func_printf(__func__, "total %lu tasks\n", walkData.idList.size());
    kt_for(nThreads, ktf_walk, &walkData, walkData.idList.size());

    // fasta格式的输出
    err_func_printf(__func__, "writing results to %s\n", fullOutputFileName.c_str());
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    for (size_t i = 0; i < walkData.outputList.size(); ++i) {
        err_fprintf(outputFile, walkData.outputList[i].c_str());
    }
    err_fclose(outputFile);
    err_func_printf(__func__, "done\n");

    return 0;
}
