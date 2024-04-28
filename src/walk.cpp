#include <vector>
#include <string>
#include <set>
#include "utils.hpp"
#include "FastaReader/FastaReader.hpp"

extern unsigned char nst_nt4_table[256];

enum class MerType {x, i, o};
enum class BrcType {good, bad};

int walk_core(const std::string &kFileName, const std::string &okFileName, const std::string &faFileName, const std::string &outputFileName) {

    FastaReader fr(faFileName);
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

    FILE *outputFile = xopen((outputFileName + ".brc").c_str(), "wb");
    setvbuf(outputFile, NULL, _IOFBF, CommonFileBufSize);

    // 加载okmer到内存
    err_func_printf(__func__, "loading %s\n", okFileName.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
    set<uint64_t> okmerSet;
    progressbar bar(okmerNum);
    for (uint64_t i = 0; i < okmerNum; ++i) {
        uint64_t kmer;
        err_fread_noeof(&kmer, sizeof(uint64_t), 1, okFile);
        okmerSet.insert(kmer);
        bar.update();
    }
    bar.end();

    while (fr.hasNextSequence()) {
        // 循环处理输入fa中的每个sequence
        fr.readSequence();
        std::string seq = fr.Sequence();
        uint64_t brc_id = 0;
        uint64_t brc_length = 0;
        err_func_printf(__func__, "processing %s\n", fr.Id().c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
        progressbar bar(seq.length() - k + 1);
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
            else if (okmerSet.find(kmer) == okmerSet.end()) merType = MerType::i;
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
                // 输出brc到文件
                err_fprintf(outputFile, ">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, fr.Id().c_str(), brc.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
                // 初始化下一个brc到内存中
                brcType = BrcType::bad;
                brc_length = k;
                brc_id ++;
                brc = brc.substr(brc.size() - k + 1);
                brc += seq[i + k - 1];
            }
            else if (brcType == BrcType::bad && merType == MerType::o) {
                // 输出brc到文件
                err_fprintf(outputFile, ">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, fr.Id().c_str(), brc.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
                // 初始化下一个brc到内存中
                brcType = BrcType::good;
                brc_length = k;
                brc_id ++;
                brc = brc.substr(brc.size() - k + 1);
                brc += seq[i + k - 1];
            }
            else if (brcType == BrcType::bad && merType == MerType::i) {
                // 输出brc到文件
                err_fprintf(outputFile, ">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, fr.Id().c_str(), brc.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
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
            bar.update();
        }
        // 输出最后一个brc到文件
        err_fprintf(outputFile, ">%ld|%ld|%s\n%s\n\n", brc_length, brc_id, fr.Id().c_str(), brc.c_str()); // 输出的fasta头包含该条brc恢复后应该的长度，以及该brc是所属seq的第几条brc，以及该brc原来所属seq的id
        bar.end();
    }

    // err_fclose(kFile);
    err_fclose(okFile);
    err_fclose(outputFile);
    return 0;
}
