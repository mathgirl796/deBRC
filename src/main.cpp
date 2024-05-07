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

using namespace std;

int kmc_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<int>("maxRamGB", 'r', "max amount of RAM in GB", false, 8, cmdline::range(1, 65535));
    a.add<int>("nThreads", 't', "total number of threads", false, 4, cmdline::range(1, 65535));
    a.add<int>("kmerLen", 'k', "k-mer length", false, 25, cmdline::range(2, 32));
    a.add<string>("tmpPath", 'w', "place to hold temp files", false, ".");
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
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
            a.get<string>("outputFileName")
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
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.footer("inputFileNames");
    a.parse_check(argc, argv);
    if (a.rest().size() < 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return merge_core(
            a.rest(), a.get<string>("outputFileName"),
            true, 1
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
    a.add<string>("tmpPath", 'w', "place to hold temp files", false, ".");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return sort_core(
            a.rest()[0], a.get<string>("outputFileName"),
            a.get<string>("tmpPath"), a.get<int>("maxRamGB"), a.get<int>("nThreads")
        );
    }
}

int view_main(int argc, char *argv[]) {
    cmdline::parser a;
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return view_core(a.rest()[0]);
    }
}


int check_main(int argc, char *argv[]) {
    cmdline::parser a;
    a.add<int>("nThreads", 't', "total number of threads", false, 4, cmdline::range(1, 65535));
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return check_core(a.rest()[0], a.get<int>("nThreads"));
    }
}

int split_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add("ksmer", '\0', "need .k.smer");
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
            a.exist("ksmer")
        );
    }
}

int walk_main(int argc, char *argv[])
{
    cmdline::parser a;
    a.add<string>("kFileName", 'k', "output file path (without extension)", false, "");
    a.add<string>("okFileName", 'l', "output file path (without extension)", true);
    a.add<string>("outputFileName", 'o', "output file path (without extension)", true);
    a.add<int>("nThreads", 't', "number of workers dealing with seqs", false, 8);
    a.add("passSpecialCharactors", '\0', "don't store bad brc which has unknown charactors");
    a.footer("inputFileName");
    a.parse_check(argc, argv);
    if (a.rest().size() != 1) {
        cerr << a.usage();
        return 1;
    }
    else {
        return walk_core(
            a.get<string>("kFileName"),
            a.get<string>("okFileName"),
            a.rest()[0],
            a.get<string>("outputFileName"),
            a.get<int>("nThreads"),
            a.exist("passSpecialCharactors")
        );
    }
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
    // h.add_function("extract", "proper version of kmc&convert", extract_main);
    h.add_function("view", "view .mer or .smer file", view_main);
    h.add_function("check", "check .mer or .smer file is no-decrease or not", check_main);
    h.add_function("split", "split kp1.smer to k.smer & o.smer", split_main);
    h.add_function("walk", "generate brc using k.smer(optional) & o.smer", walk_main);
    int ret = h.run(argc, argv);

    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0)
    {
        err_func_printf(__func__, "Memory usage: %lfGB\n", double(usage.ru_maxrss) / 1024 / 1024);
    }
    else
    {
        err_func_printf(__func__, "Failed to get memory usage\n");
    }
    return 0;

    return ret;
}