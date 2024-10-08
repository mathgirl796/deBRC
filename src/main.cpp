#include <iostream>
#include <cstring>
#include <string>
#include <bitset>
#include <sys/resource.h>
#include "cmdline.hpp"
#include "utils.hpp"

#include "kmc.hpp"
#include "convert.hpp"
#include "sort.hpp"
#include "merge.hpp"
#include "others.hpp"
#include "split.hpp"
#include "walk.hpp"
#include "restore.hpp"
#include "unitig.hpp"
#include "okmerflag.hpp"

using namespace std;

int kmc_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<int>("maxRamGB", 'r', "max amount of RAM in GB", false, 8, cmdline::range(1, 65535));
    a.add<int>("nThreads", 't', "total number of threads", false, 4, cmdline::range(1, 65535));
    a.add<int>("kmerLen", 'k', "k-mer length", false, 25, cmdline::range(2, 32));
    a.add<string>("tmpPath", 'w', "place to hold temp files", false, "./");
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add("canonical", '\0', "if count canonical-form kmers");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() < 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return kmc_core(
            a.rest(), a.get<string>("tmpPath"),
            a.get<int>("maxRamGB"), a.get<int>("nThreads"), a.get<int>("kmerLen"),
            a.get<string>("outputFileName"), a.exist("canonical")
        );
    }
}

int convert_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return convert_core(
            a.rest()[0],
            a.get<string>("outputFileName")
        );
    }
}

int merge_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("outputFileName", 'o', "output file path", true);
    a.add("distinct", 'd', "if make result distinct");
    a.footer("inputFileNames");
    a.parse_check(argc, argv);
    if (a.rest().size() < 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return merge_core(
            a.rest(), a.get<string>("outputFileName"),
            a.exist("distinct"), 1
        );
    }
}

int compare_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("refFileName", 'r', "reference .smer file", true);
    a.add<string>("tgtFileName", 't', "target .smer file", true);
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.parse_check(argc, argv);
    return compare_core(
        a.get<string>("refFileName"), a.get<string>("tgtFileName"),
        a.get<string>("outputFileName")
    );
}

int sort_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add<int>("maxRamGB", 'r', "max amount of RAM in GB", false, 8, cmdline::range(1, 65535));
    a.add<int>("nThreads", 't', "total number of threads", false, 4, cmdline::range(1, 65535));
    a.add<string>("tmpPath", 'w', "place to hold temp files", false, "./");
    a.add("distinct", '\0', "delete duplicate kmers");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return sort_core(
            a.rest()[0], a.get<string>("outputFileName"),
            a.get<string>("tmpPath"), a.get<int>("maxRamGB"), a.get<int>("nThreads"), 
            a.exist("distinct")
        );
    }
}

int view_main(int argc, char *argv[]) {
    cmdline::parser a;
    a.add("head", 'h', "only view kmerLength & kmerNum");
    a.add("tail", 't', "view from tail");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return view_core(a.rest()[0], a.exist("head"), a.exist("tail"));
    }
}


int check_main(int argc, char *argv[]) {
    cmdline::parser a;
    a.add<int>("nThreads", 't', "total number of threads", false, 4, cmdline::range(1, 65535));
    a.add("check_dup", 'd', "check duplicate rather than only decrease");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return check_core(a.rest()[0], a.get<int>("nThreads"), a.exist("check_dup"));
    }
}

int split_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add<int>("maxRamGB", 'r', "max amount of RAM in GB", false, 8, cmdline::range(1, 65535));
    a.add<int>("nThreads", 't', "total number of threads", false, 4, cmdline::range(1, 65535));
    a.add<string>("tmpPath", 'w', "place to hold temp files", false, "./");
    a.add("ismer", '\0', "also output .i.smer");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return split_core(
            a.rest()[0],
            a.get<string>("outputFileName"),
            a.exist("ismer"),
            a.get<string>("tmpPath"), a.get<int>("maxRamGB"), a.get<int>("nThreads")
        );
    }
}

int walk_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("smerFileName", 's', ".smer file path(kp1)", false, "");
    a.add<string>("okFileName", 'l', ".o.mer file path", true);
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add<int>("nThreads", 't', "number of workers dealing with seqs", false, 8);
    a.add<int>("maxRamGB", 'r', "max amount of RAM in GB (not very smart, at least be okmerFileSize+biggestChr*2+secondBigChr*2+severalRedundantGB)", false, 32, cmdline::range(1, 65535));
    a.add("passSpecialCharactors", '\0', "don't store bad brc which has unknown charactors");
    a.add("useKmerFormat", '\0', "output kmer format");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() < 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return walk_core(
            a.get<string>("smerFileName"),
            a.get<string>("okFileName"),
            a.rest(),
            a.get<string>("outputFileName"),
            a.get<int>("nThreads"),
            a.exist("passSpecialCharactors"),
            a.exist("useKmerFormat"),
            a.get<int>("maxRamGB")
        );
    }
}

int restore_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("smerFileName", 's', ".smer(kp1) file path", false, "");
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add<int>("nThreads", 't', "number of workers dealing with seqs", false, 8);
    a.add("useKmerFormat", '\0', "output kmer format");
    a.footer("brcFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return restore_core(
            a.get<string>("smerFileName"),
            a.rest()[0],
            a.get<string>("outputFileName"),
            a.get<int>("nThreads"),
            a.exist("useKmerFormat")
        );
    }
}

int search_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("smerFileName", 's', ".smer file path", true);
    a.add<string>("kmer", 'm', "kmer", false, "");
    a.parse_check(argc, argv);
    return search_core(
        a.get<string>("smerFileName"),
        a.get<string>("kmer")
    );
}

int unitig_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("ikFileName", 'i', "multi-in kp1mer set", true);
    a.add<string>("okFileName", 'k', "multi-out kp1mer set", true);
    a.add<string>("outputFileName", 'o', "output file path (without extention)", true);
    a.add<int>("nThreads", 't', "number of workers dealing with seqs", false, 4);
    a.add<int>("maxRamGB", 'r', "max amount of RAM in GB (not very smart, at least be ikmerFileSize+okmerFileSize+biggestChr*2+secondBigChr*2+severalRedundantGB)", false, 8);
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else
        return unitig_core(
            a.get<string>("ikFileName"), a.get<string>("okFileName"), a.rest()[0], a.get<string>("outputFileName"), 
            a.get<int>("nThreads"), a.get<int>("maxRamGB")
        );
}

int okmerflag_main(int argc, char *argv[]) {
    cmdline::parser a;
    a.add<string>("outputFileName", 'o', "output file path (without extention)", true);
    a.add<string>("okFileName", 'k', "multi-out kp1mer set", true);
    a.add<int>("nThreads", 't', "number of workers dealing with seqs", false, 4);
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    return okmerflag_core(
        a.rest()[0], a.get<string>("okFileName"), a.get<string>("outputFileName"), 
        a.get<int>("nThreads")
    );
}

int main(int argc, char *argv[])
{
    err_func_printf(__func__, "program start, argv:\n");
    for (int i = 0 ; i < argc; ++i) {
        fprintf(stderr, "%s ", argv[i]);
    }fprintf(stderr, "\n");

    COMMAND_HANDLER h;
    h.add_function("kmc", "count kmers in one/several fasta files", kmc_main);
    h.add_function("convert", "convert file format from .kmc_* to .mer", convert_main);
    h.add_function("merge", "merge consistent .smer files, duplicate values will be removed", merge_main);
    h.add_function("compare", "compare two .smer files", compare_main);
    h.add_function("sort", "sort file format from .mer to .smer", sort_main);
    h.add_function("view", "view .mer or .smer file", view_main);
    h.add_function("check", "check .mer or .smer file is no-decrease or not", check_main);
    h.add_function("split", "analysis kp1.smer to generate o.smer", split_main);
    h.add_function("walk", "generate brc using fasta & o.smer & .smer(optional)", walk_main);
    h.add_function("restore", "restore fasta using k.mer & brc", restore_main);
    h.add_function("search", "search a kmer in an smer file", search_main);
    h.add_function("unitig", "represent fasta in consecutive unitigs format", unitig_main);
    h.add_function("okmerflag", "flag a fasta file indicate which pos is brc", okmerflag_main);
    int ret = h.run(argc, argv);

    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0)
    {
        err_func_printf(__func__, "Memory usage: %lfGB\n", double(usage.ru_maxrss) / OneMega);
    }
    else
    {
        err_func_printf(__func__, "Failed to get memory usage\n");
    }
    exit(ret);
}