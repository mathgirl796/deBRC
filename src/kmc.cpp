#include <string>
#include "kmc_api/kmc_runner.h"
#include "kmc_api/kmc_file.h"
#include "utils.hpp"

int kmc_core(const std::vector<std::string> &inputFiles, const std::string &tmpPath, uint32_t maxRamGB, uint32_t nThreads, uint32_t kmerLen, const std::string &outputFileName) {
    try
    {       
        KMC::Runner runner;

        KMC::Stage1Params stage1Params;
        stage1Params
            .SetInputFileType(KMC::InputFileType::MULTILINE_FASTA) // must invoke this function to tell kmc the input type
            .SetCanonicalKmers(false)
            .SetInputFiles(inputFiles)
            .SetTmpPath(tmpPath)
            .SetMaxRamGB(maxRamGB)
            .SetNThreads(nThreads)
            .SetKmerLen(kmerLen)
            ;
        
        auto stage1Result = runner.RunStage1(stage1Params);
        
        KMC::Stage2Params stage2Params;

        stage2Params
            .SetCutoffMin(1)
            .SetOutputFileType(KMC::OutputFileType::KMC)
            .SetOutputFileName(outputFileName)
            .SetMaxRamGB(maxRamGB)
            .SetNThreads(nThreads)
            ;

        auto stage2Result = runner.RunStage2(stage2Params);

        //print some stats
        std::cerr << "total k-mers: " << stage2Result.nTotalKmers << "\n";
        std::cerr << "total unique k-mers: " << stage2Result.nUniqueKmers << "\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
