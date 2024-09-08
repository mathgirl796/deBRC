#include <vector>
#include <string>
#include <set>
#include <sys/resource.h>
#include <unordered_set>
#include <algorithm>
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

void ktf_walk(void* shared, long i, int tid) {
    WalkData *data = (WalkData *)shared;
    // 处理输入中的一个seq
    err_func_printf(__func__, "worker %d start work, processing %s (job %ld)\n", tid, data->idList[i].c_str(), i); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    string id = data->idList[i];
    string seq = data->seqList[i];
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    string output = "";
    uint32_t kp1 = data->kp1;
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

    // seq长度过短，直接视为一条bad brc，结束循环
    if (seq.length() < kp1) {
        brc_id = 0;
        original_seq_length = seq.length();
        startPosOnSeq = 0;
        brcType = BrcType::bad;
        brc = seq;
        goto END_OF_SINGLE_WALK;
    }

    // 遍历所有kp1mer
    for (size_t i = 0; i < seq.length() - kp1 + 1; ++i) {

        /* 一条seq头一个kp1mer的处理！（初始化读入kmer） */
        if (i == 0) {
            for (uint32_t j = 0; j < kp1 - 1; ++j) { // 读入头一个kmer的前k个base
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
        else if (data->okmerSet.find(kp1mer) == data->okmerSet.end()) merType = MerType::i;
        else merType = MerType::o;

        /* 一条seq头一个kp1mer的处理！（初始化BrcType） */
        if (i == 0) {
            brcType = (merType == MerType::x) ? (BrcType::bad) : (BrcType::good);
        }

        // err_printf("%d\t %d\n", brcType, merType);

        if (brcType == BrcType::good && merType == MerType::o) {
            if (data->useKmerFormat) brc += seq.substr(i, k + 1);
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
            if (data->passSpecialCharactors != true)
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
            if (data->passSpecialCharactors != true)
                output += string_format(">%lu|%lu|%lu|%s|%s\n%s\n\n", brc_id, original_seq_length - k, startPosOnSeq, (brcType == BrcType::good) ? ("good"):("bad"), id.c_str(), brc.substr(0, brc.length() - k).c_str());
            // 确定下一个brc的类型
            brcType = BrcType::good;
            // 初始化下一个brc到内存中
            brc_id ++;
            brc = seq.substr(i, k);
            if (data->useKmerFormat) brc += seq.substr(i, k + 1);
            else brc += seq[i + kp1 - 1];
            original_seq_length = kp1;
            startPosOnSeq = i;
        }
        else if (brcType == BrcType::bad && merType == MerType::x) {
            brc += seq[i + kp1 - 1];
            original_seq_length ++;
        }
    }

END_OF_SINGLE_WALK:
    // 输出最后一个brc到文件
    if (brcType == BrcType::good || (brcType == BrcType::bad && data->passSpecialCharactors != true)) {
        output += string_format(">%lu|%lu|%lu|%s|%s\n%s\n\n", brc_id, original_seq_length, startPosOnSeq, (brcType == BrcType::good) ? ("good"):("bad"), id.c_str(), brc.c_str());
    }
    data->outputList[i] = output;
    err_func_printf(__func__, "worker %d finish work\n", tid);
}

int walk_core(const std::string &smerFileName, const std::string &okFileName, 
    const std::vector<std::string> &faFiles, const std::string &outputFileName, uint32_t nThreads,
    bool passSpecialCharactors, bool useKmerFormat, int maxRamGB) {

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

    // 创建多线程共享数据结构
    struct WalkData data;
    data.kp1 = kp1;
    data.passSpecialCharactors = passSpecialCharactors;
    data.useKmerFormat = useKmerFormat;

    // 加载okmer到内存
    err_func_printf(__func__, "loading %s\n", okFileName.c_str());
    data.okmerSet.reserve(omerNum);
    for (uint64_t i = 0; i < omerNum; ++i) {
        uint64_t kp1mer;
        err_fread_noeof(&kp1mer, sizeof(uint64_t), 1, omer);
        data.okmerSet.insert(kp1mer);
    }
    err_fclose(omer);
    err_func_printf(__func__, "done, total %lu omers\n", data.okmerSet.size());

    // 创建线程池
    void *pool = kt_forpool_init(nThreads);

    // 估算剩余内存
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    double maxRam = (double)(maxRamGB - 1) * OneGiga - (double)usage.ru_maxrss * OneKilo; // 最后减1G是因为人最长染色体0.25G，加上输出brc是0.5G，多给1G应该就不会超标了
    err_func_printf(__func__, "total ram: %lf GB, used %lf GB, io buffer remaining %lf GB\n", (double)maxRamGB, (double)usage.ru_maxrss / OneMega, maxRam / OneGiga);

    // 打开输出文件
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);

    // 处理每条输入文件
    for (auto faFileName : faFiles) {
        gzFile fp = xzopen(faFileName.c_str(), "r");
        gzbuffer(fp, CommonFileBufSize);
        kseq_t *seqs = kseq_init(fp);

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
            kt_forpool(pool, ktf_walk, &data, taskNum);

            for (size_t i = 0; i < data.outputList.size(); ++i) {
                err_fprintf(outputFile, data.outputList[i].c_str());
            }
        }

        // 关闭输入文件
        kseq_destroy(seqs);
        err_gzclose(fp);
    }

    // 关闭输出文件
    err_fclose(outputFile);

    return 0;
}
