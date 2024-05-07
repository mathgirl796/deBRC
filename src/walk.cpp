#include <vector>
#include <string>
#include <set>
#include "utils.hpp"
#include "FastaReader/FastaReader.hpp"
#include "klib/kthread.hpp"

extern unsigned char nst_nt4_table[256];

enum class MerType {x, i, o};
enum class BrcType {good, bad};

struct WalkData {
    vector<string> seqList;
    vector<string> idList;
    vector<string> outputList;
    set<uint64_t> okmerSet;
    uint32_t k;
    bool passSpecialCharactors;
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
            brc += seq[i + k - 1];
            brc_length ++;
        }
        else if (brcType == BrcType::good && merType == MerType::i) {
            brc_length ++;
        }
        else if (brcType == BrcType::good && merType == MerType::x) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            output.append(string_format(">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, id.c_str(), brc.c_str()));
            // 初始化下一个brc到内存中
            brcType = BrcType::bad;
            brc_length = k;
            brc_id ++;
            brc = brc.substr(brc.size() - k + 1);
            brc += seq[i + k - 1];
        }
        else if (brcType == BrcType::bad && merType == MerType::o) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            if (walkData->passSpecialCharactors != true)
                output.append(string_format(">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, id.c_str(), brc.c_str()));
            // 初始化下一个brc到内存中
            brcType = BrcType::good;
            brc_length = k;
            brc_id ++;
            brc = brc.substr(brc.size() - k + 1);
            brc += seq[i + k - 1];
        }
        else if (brcType == BrcType::bad && merType == MerType::i) {
            // 输出brc到文件, 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
            if (walkData->passSpecialCharactors != true)
                output.append(string_format(">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, id.c_str(), brc.c_str()));
            // 初始化下一个brc到内存中
            brcType = BrcType::good;
            brc_length = k;
            brc_id ++;
            brc = brc.substr(brc.size() - k + 1);
        }
        else if (brcType == BrcType::bad && merType == MerType::x) {
            brc += seq[i + k - 1];
            brc_length ++;
        }
    }
    // 输出最后一个brc到文件
    if (brcType == BrcType::good || walkData->passSpecialCharactors != true)
        output += string_format(">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, id.c_str(), brc.c_str());
    walkData->outputList[i] = output;
    err_func_printf(__func__, "worker %d finish work\n", tid);
}

int walk_core(const std::string &kFileName, const std::string &okFileName, 
    const std::string &faFileName, const std::string &outputFileName, uint32_t nThreads,
    bool passSpecialCharactors) {

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

    string fullOutputFileName = outputFileName + ".brc";
    FILE *outputFile = xopen(fullOutputFileName.c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);

    // 创建ktf worker所需数据结构
    struct WalkData walkData;
    walkData.k = k;
    walkData.outputList = vector<string>();
    walkData.passSpecialCharactors = passSpecialCharactors;
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

    kt_for(nThreads, ktf_walk, &walkData, walkData.idList.size());

    err_func_printf(__func__, "writing results to %s\n", fullOutputFileName.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    progressbar bar2(walkData.outputList.size());
    for (size_t i = 0; i < walkData.outputList.size(); ++i) {
        err_fprintf(outputFile, walkData.outputList[i].c_str());
        bar2.update();
    }
    bar2.end();

    err_fclose(okFile);
    err_fclose(outputFile);
    return 0;
}
