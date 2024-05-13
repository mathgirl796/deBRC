#include <vector>
#include <string>
#include <set>
#include "utils.hpp"
#include "FastaReader/FastaReader.hpp"
#include "klib/kthread.hpp"

// TODO: 用kt_pipeline改写成读算写三阶段，这样就可以处理大输入文件了

extern unsigned char nst_nt4_table[256];

enum class MerType {x, i, o};
enum class BrcType {good, bad};

struct WalkData {
    vector<string> seqList;
    vector<string> idList;
    set<uint64_t> okmerSet;
    uint32_t k;
    bool passSpecialCharactors;

    // fastaFormat output
    vector<string> outputList;

    // kmerFormat output
    bool useKmerFormat;
};

void ktf_walk(void* data, long i, int tid) {
    WalkData *walkData = (WalkData *)data;
    // 循环处理输入fa中的每个sequence
    err_func_printf(__func__, "worker %d start work, processing %s (job %ld)\n", tid, walkData->idList[i].c_str(), i); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    string id = walkData->idList[i];
    string seq = walkData->seqList[i];
    string output = "";
    uint32_t k = walkData->k;
    uint64_t brc_id = 0;
    uint64_t brc_length = 0;
    uint64_t kmer = 0;
    uint64_t startPosOnSeq = 0;
    BrcType brcType;
    MerType merType;
    int badBasePos = -1;
    string brc;
    for (size_t i = 0; i < seq.length() - k + 1; ++i) {

        if (i == 0) { /* 一条seq头一个kmer的处理！ */
            for (uint32_t j = 0; j < k - 1; ++j) { // 头一个kmer的前k-1个base
                kmer = (kmer << 2) | nst_nt4_table[(int)seq[j]];
                brc += seq[j]; // 一条sequence的头k-1个base一定要写进去
                brc_length ++;
            }
            startPosOnSeq = 0;
        }
        // 读入一个base
        unsigned char newBase = nst_nt4_table[(int)seq[i + k - 1]];
        // 处理特殊base
        if (newBase >= 4) {
            badBasePos = k - 1;
        }
        else {
            badBasePos --;
        }

        kmer = ((kmer << 2) | (newBase & 0x3)) & (~(((uint64_t)-1) << (k << 1)));
        if (badBasePos >= 0) merType = MerType::x;
        else if (walkData->okmerSet.find(kmer) == walkData->okmerSet.end()) merType = MerType::i;
        else merType = MerType::o;

        if (i == 0) { /* 一条seq头一个kmer的处理！ */
            brcType = (merType == MerType::x) ? (BrcType::bad) : (BrcType::good);
        }

        // err_printf("%d\t %d\n", brcType, merType);

        if (brcType == BrcType::good && merType == MerType::o) {
            if (walkData->useKmerFormat) brc += seq.substr(i, k);
            else brc += seq[i + k - 1];
            brc_length ++;
        }
        else if (brcType == BrcType::good && merType == MerType::i) {
            brc_length ++;
        }
        else if (brcType == BrcType::good && merType == MerType::x) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            output += string_format(">%lu|%lu|%lu|%s\n%s\n\n", brc_id, brc_length, startPosOnSeq, id.c_str(), brc.c_str());
            // 确定下一个brc的类型
            brcType = BrcType::bad;
            // 初始化下一个brc到内存中
            brc_length = k;
            brc_id ++;
            brc = seq.substr(i, k);
            startPosOnSeq = i;
        }
        else if (brcType == BrcType::bad && merType == MerType::o) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            output += string_format(">%lu|%lu|%lu|%s\n%s\n\n", brc_id, brc_length, startPosOnSeq, id.c_str(), brc.c_str());
            // 确定下一个brc的类型
            brcType = BrcType::good;
            // 初始化下一个brc到内存中
            brc_length = k;
            brc_id ++;
            brc = seq.substr(i, k);
            startPosOnSeq = i;
            // 直接保存当前的omer
            if (walkData->useKmerFormat) brc += seq.substr(i, k);
            else brc += seq[i + k - 1];
        }
        else if (brcType == BrcType::bad && merType == MerType::i) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            output += string_format(">%lu|%lu|%lu|%s\n%s\n\n", brc_id, brc_length, startPosOnSeq, id.c_str(), brc.c_str());
            // 确定下一个brc的类型
            brcType = BrcType::good;
            // 初始化下一个brc到内存中
            brc_length = k;
            brc_id ++;
            brc = seq.substr(i, k);
            startPosOnSeq = i;
        }
        else if (brcType == BrcType::bad && merType == MerType::x) {
            brc += seq[i + k - 1];
            brc_length ++;
        }
    }
    // 输出最后一个brc到文件
    if (brcType == BrcType::good || walkData->passSpecialCharactors != true) {
        output += string_format(">%lu|%lu|%lu|%s\n%s\n\n", brc_id, brc_length, startPosOnSeq, id.c_str(), brc.c_str());
        walkData->outputList[i] = output;
    }
    err_func_printf(__func__, "worker %d finish work\n", tid);
}

int walk_core(const std::string &kFileName, const std::string &okFileName, 
    const std::string &faFileName, const std::string &outputFileName, uint32_t nThreads,
    bool passSpecialCharactors, bool useKmerFormat) {

    uint32_t ok = 0;
    uint64_t okmerNum = 0;

    FILE *okFile = xopen(okFileName.c_str(), "rb");
    setvbuf(okFile, NULL, _IOFBF, CommonFileBufSize);
    err_fread_noeof(&ok, sizeof(uint32_t), 1, okFile);
    err_fread_noeof(&okmerNum, sizeof(uint64_t), 1, okFile);

    uint32_t k = ok;
    // uint64_t kmerNum = 0;
    // FILE *kFile = xopen(kFileName.c_str(), "rb");
    // setvbuf(kFile, NULL, _IOFBF, CommonFileBufSize);
    // err_fread_noeof(&k, sizeof(uint32_t), 1, kFile);
    // err_fread_noeof(&kmerNum, sizeof(uint64_t), 1, kFile);

    if (k != ok) {
        perror("inputed kFile is not compatable with okFile");
    }

    string fullOutputFileName = outputFileName;
    if (passSpecialCharactors) fullOutputFileName += ".passN"; // 若不保留未知字符导致的xBrc，则加后缀注明
    fullOutputFileName += ".brc";

    // 创建ktf worker所需数据结构
    struct WalkData walkData;
    walkData.k = k;
    walkData.outputList = vector<string>();
    walkData.passSpecialCharactors = passSpecialCharactors;
    walkData.useKmerFormat = useKmerFormat;
    gzFile fp = xzopen(faFileName.c_str(), "r");
    kseq_t *seqs = kseq_init(fp);
    while (kseq_read(seqs) >= 0) {
        // fprintf(stderr, "%s | %s\n", seqs->name.s, seqs->seq.s);
        walkData.seqList.push_back(string(seqs->seq.s));
        walkData.idList.push_back(string(seqs->name.s));
        walkData.outputList.push_back("");
    }
    gzclose(fp);

    // 加载okmer到内存
    err_func_printf(__func__, "loading %s\n", okFileName.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    set<uint64_t> okmerSet;
    progressbar bar(okmerNum);
    for (uint64_t i = 0; i < okmerNum; ++i) {
        uint64_t kmer;
        err_fread_noeof(&kmer, sizeof(uint64_t), 1, okFile);
        walkData.okmerSet.insert(kmer);
        bar.update();
    }
    bar.end();

    // 执行walk任务
    err_func_printf(__func__, "total %lu tasks\n", walkData.idList.size());
    kt_for(nThreads, ktf_walk, &walkData, walkData.idList.size());

    // fasta格式的输出
    err_func_printf(__func__, "writing results to %s\n", fullOutputFileName.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);
    progressbar bar2(walkData.outputList.size());
    for (size_t i = 0; i < walkData.outputList.size(); ++i) {
        err_fprintf(outputFile, walkData.outputList[i].c_str());
        bar2.update();
    }
    bar2.end();
    err_fclose(outputFile);

    err_fclose(okFile);
    return 0;
}
