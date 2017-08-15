// ============================================================================
//                                 CNVetti
// ============================================================================
// Copyright (C) 2016-2017 Berlin Institute for Health
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 3 of the License, at (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
//
// ============================================================================
// Author:  Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>
// ============================================================================

#include "cnvetti/program_options.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/synced_bcf_reader.h>

#include "cnvetti/program_options.h"
#include "cnvetti/version.h"
#include "cnvetti/utils.h"

namespace {  // anonymous namespace

const double MAD_RES = 10000.0;

// ---------------------------------------------------------------------------
// Struct MedianStats
// ---------------------------------------------------------------------------

struct MedianStats
{
    int gcContent;
    int windowCount;
    double cov;
    double cov0;
    double rc;
    double rc0;
};

// ---------------------------------------------------------------------------
// Struct MadHistoEntry
// ---------------------------------------------------------------------------

struct MadHistoEntry
{
    MadHistoEntry() : cov(0), cov0(0), rc(0), rc0(0)
    {}

    int cov;
    int cov0;
    int rc;
    int rc0;
};

// ---------------------------------------------------------------------------
// Class CnvettiNormalizeApp
// ---------------------------------------------------------------------------

class CnvettiNormalizeApp
{
public:
    CnvettiNormalizeApp(CnvettiNormalizeOptions const & options) :
            options(options), vcfHeaderIn(nullptr), vcfHeaderOut(nullptr)
    {}

    void run();

private:

    // Open the input files and store them in bamFilesIn and open vcfFileOutPtr;
    void openFiles();
    // Check that the references FASTA file fits to the BAM file.
    void checkFiles();
    // Load the per-GC medians from "__cnvetti_stats" contig
    void loadMedians();
    // Perform chromosome-wise processing
    void processChromosomes();
    // Process a region
    void processRegion(std::string const & contig, int beginPos, int endPos);
    // Print statistics such as MAD.
    void printStatistics(std::ostream & out);

    CnvettiNormalizeOptions const & options;

    // medianStats[sampleID][gcContent]
    std::vector<
        std::map<int, MedianStats>
    > medianStats;

    // Median absolute deviation of normalized value.
    //
    // madHisto[sampleID][(int)(MAD * 100)]
    std::vector<
        std::map<int, MadHistoEntry>
    > madHistos;

    std::shared_ptr<bcf_srs_t> vcfReaderPtr;
    std::shared_ptr<vcfFile> vcfFileOutPtr;
    std::shared_ptr<bcf_hdr_t> vcfHeaderOutPtr;
    bcf_hdr_t * vcfHeaderIn;
    bcf_hdr_t * vcfHeaderOut;
};


void CnvettiNormalizeApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti normalize\n"
                  << "=================\n\n";
        options.print(std::cerr);
    }

    if (options.verbosity >= 1)
        std::cerr << "\nOPENING FILES\n\n";
    openFiles();

    if (options.verbosity >= 1)
        std::cerr << "\nLOADING SUMMARIES\n\n";
    loadMedians();

    if (options.verbosity >= 1)
        std::cerr << "\nPROCESSING\n\n";
    processChromosomes();

    if (options.verbosity >= 1)
        printStatistics(std::cerr);

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiNormalizeApp::openFiles()
{
    // Open input VCF file and read header
    if (options.verbosity >= 1)
        std::cerr << "Opening " << options.inputFileName << " ...\n";
    vcfReaderPtr = std::shared_ptr<bcf_srs_t>(bcf_sr_init(), bcf_sr_destroy);
    bcf_sr_set_opt(vcfReaderPtr.get(), BCF_SR_REQUIRE_IDX);
    if (bcf_sr_set_threads(vcfReaderPtr.get(), options.numIOThreads) < 0)
        throw std::runtime_error("Failed to create threads");
    if (!bcf_sr_add_reader(vcfReaderPtr.get(), options.inputFileName.c_str()))
        throw std::runtime_error(std::string("Failed to open ") + options.inputFileName);
    vcfHeaderIn = vcfReaderPtr.get()->readers[0].header;

    // Check input VCF file
    checkFiles();

    // Prepare file opening, flags etc.
    if (options.verbosity >= 1)
        std::cerr << "Opening " << options.outputFileName << "\n";
    enum htsExactFormat openFormat = vcf;
    enum htsCompression openCompression = no_compression;
    std::string openMode = "w";
    if (hasSuffix(options.outputFileName, ".bcf")) {
        openMode += "b";
        openFormat = bcf;
        openCompression = bgzf;
    } else if (hasSuffix(options.outputFileName, ".vcf.gz")) {
        openMode += "b";
        openFormat = vcf;
        openCompression = bgzf;
    }
    // Actually open file, apply flags
    vcfFileOutPtr = std::shared_ptr<vcfFile>(
        hts_open(options.outputFileName.c_str(), openMode.c_str()),
        [&](vcfFile * ptr) {
            hts_close(ptr);
            if (hasSuffix(options.outputFileName, ".vcf.gz"))
            {
                if (tbx_index_build(options.outputFileName.c_str(), 0, &tbx_conf_vcf))
                    throw std::runtime_error("Problem writing TBI index!");
            }
            else if (hasSuffix(options.outputFileName, ".bcf"))
            {
                if (bcf_index_build(options.outputFileName.c_str(), 14))
                    throw std::runtime_error("Problem writing CSI index!");
            }
        });
    hts_set_threads(vcfFileOutPtr.get(), options.numIOThreads);
    vcfFileOutPtr->format.format = openFormat;
    vcfFileOutPtr->format.compression = openCompression;

    // Copy input to output header
    vcfHeaderOutPtr = std::shared_ptr<bcf_hdr_t>(bcf_hdr_dup(vcfHeaderIn), bcf_hdr_destroy);
    vcfHeaderOut = vcfHeaderOutPtr.get();
    // Remove pseudo-contig with statistics and FILTER
    bcf_hdr_remove(vcfHeaderOut, BCF_HL_CTG, "__cnvetti_stats");
    bcf_hdr_remove(vcfHeaderOut, BCF_HL_FLT, "CNVETTI_STATS");

    bcf_hdr_printf(vcfHeaderOut,
        ("##FILTER=<ID=FEW_GCWINDOWS,Description="
         "\"Masked because of few windows with this GC content\">"));

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeaderOut, "##cnvetti_normalizeCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOutPtr.get(), vcfHeaderOut);
}

void CnvettiNormalizeApp::checkFiles()
{
    // Check that the file looks like created by "cnvetti coverage"
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_FLT, "ID", "CNVETTI_STATS", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, FILTER/CNVETTI_STATS header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_STR, "ID", "STATS", "ALT"))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, ALT/STATS header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_INFO, "ID", "END", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, INFO/END header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_INFO, "ID", "GC", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, INFO/GC header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_INFO, "ID", "MAPABILITY", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, INFO/MAPABILITY header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_INFO, "ID", "GAP", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, INFO/GAP header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_INFO, "ID", "GCWINDOWS", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, INFO/GCWINDOWS header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_FMT, "ID", "COV", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, FORMAT/COV header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_FMT, "ID", "COV0", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, FORMAT/COV0 header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_FMT, "ID", "RC", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, FORMAT/RC header missing");
    if (!bcf_hdr_get_hrec(vcfHeaderIn, BCF_HL_FMT, "ID", "RC0", NULL))
        throw std::runtime_error(
            "Input file does not look like cnvetti coverage file, FORMAT/RC0 header missing");
}

void CnvettiNormalizeApp::loadMedians()
{
    if (bcf_sr_seek(vcfReaderPtr.get(), "__cnvetti_stats", 0) != 0)
        throw std::runtime_error("Could not seek to statistics pseudo-contig");

    // Allocate median stats for all samples
    int const numSamples = bcf_hdr_nsamples(vcfHeaderIn);
    medianStats.resize(numSamples);
    madHistos.resize(numSamples);

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);

        float * ptrGC = nullptr;
        int32_t valGC;
        int32_t sizeGC = 0;
        if (bcf_get_info_float(vcfHeaderIn, line, "GC", &ptrGC, &sizeGC) == 1)
        {
            valGC = *ptrGC;
            free(ptrGC);
        }
        else
        {
            free(ptrGC);
            throw std::runtime_error("Could not get INFO/GC value");
        }

        int32_t * ptrGCWindows = nullptr;
        int32_t valGCWindows;
        int32_t sizeGCWindows = 0;
        if (bcf_get_info_int32(vcfHeaderIn, line, "GCWINDOWS", &ptrGCWindows, &sizeGCWindows) == 1)
        {
            valGCWindows = *ptrGCWindows;
            free(ptrGCWindows);
        }
        else
        {
            free(ptrGCWindows);
            throw std::runtime_error("Could not get INFO/GCWINDOWS value");
        }

        typedef std::unique_ptr<float, std::function<void(float*)>> TCleanedPtr;

        TCleanedPtr covs, cov0s, rcs, rc0s;
        int res;

        float * tmpCovs = nullptr;
        int32_t sizeCovs = 0;
        res = bcf_get_format_float(vcfHeaderIn, line, "COV", &tmpCovs, &sizeCovs);
        covs = TCleanedPtr(tmpCovs, [](float * data) { free(data); });
        if (!res)
            throw std::runtime_error("Could not get FORMAT/COV value");

        float * tmpCov0s = nullptr;
        int32_t sizeCov0s = 0;
        res = bcf_get_format_float(vcfHeaderIn, line, "COV0", &tmpCov0s, &sizeCov0s);
        cov0s = TCleanedPtr(tmpCov0s, [](float * data) { free(data); });
        if (!res)
            throw std::runtime_error("Could not get FORMAT/COV0 value");

        float * tmpRCs = nullptr;
        int32_t sizeRCs = 0;
        res = bcf_get_format_float(vcfHeaderIn, line, "RC", &tmpRCs, &sizeRCs);
        rcs = TCleanedPtr(tmpRCs, [](float * data) { free(data); });
        if (!res)
            throw std::runtime_error("Could not get FORMAT/RC value");

        float * tmpRC0s = nullptr;
        int32_t sizeRC0s = 0;
        res = bcf_get_format_float(vcfHeaderIn, line, "RC0", &tmpRC0s, &sizeRC0s);
        rc0s = TCleanedPtr(tmpRC0s, [](float * data) { free(data); });
        if (!res)
            throw std::runtime_error("Could not get FORMAT/RC0 value");

        for (int sampleID = 0; sampleID < numSamples; ++sampleID)
            medianStats[sampleID][valGC] = MedianStats {
                valGC, valGCWindows, covs.get()[sampleID], cov0s.get()[sampleID],
                rcs.get()[sampleID], rc0s.get()[sampleID]
            };
    }

    if (options.verbosity >= 1)
    {
        std::cerr << "\nCOVERAGE MEDIANS\n\n";

        std::cerr << "GC_CONTENT\tWINDOWS";
        for (int sampleNo = 0; sampleNo < bcf_hdr_nsamples(vcfHeaderIn); ++sampleNo)
            std::cerr
                << "\t" << vcfHeaderIn->samples[sampleNo] << ":COV"
                << "\t" << vcfHeaderIn->samples[sampleNo] << ":COV0"
                << "\t" << vcfHeaderIn->samples[sampleNo] << ":RC"
                << "\t" << vcfHeaderIn->samples[sampleNo] << ":RC0";
        std::cerr << "\n";

        for (int gcContent = 0; gcContent <= 100; ++gcContent)
        {
            std::cerr << gcContent;
            if (medianStats[0].count(gcContent) > 0u)
                std::cerr << "\t" << medianStats[0][gcContent].windowCount;
            else
                std::cerr << "\t0";
            for (int sampleNo = 0; sampleNo < bcf_hdr_nsamples(vcfHeaderIn); ++sampleNo)
            {
                std::cerr
                    << "\t" << medianStats[sampleNo][gcContent].cov
                    << "\t" << medianStats[sampleNo][gcContent].cov0
                    << "\t" << medianStats[sampleNo][gcContent].rc
                    << "\t" << medianStats[sampleNo][gcContent].rc0;
            }
            std::cerr << "\n";
        }
    }
}

void CnvettiNormalizeApp::processChromosomes()
{
    for (int i = 0; i < vcfHeaderIn->nhrec; ++i)  // skip pseudo contig
    {
        bcf_hrec_t * hrec = vcfHeaderIn->hrec[i];
        if (hrec->type == BCF_HL_CTG)
        {
            for (int j = 0; j < hrec->nkeys; ++j)
            {
                if (strcmp("ID", hrec->keys[j]) == 0 &&
                        strcmp("__cnvetti_stats", hrec->vals[j]) != 0)
                {
                    processRegion(hrec->vals[j], 0, 0);
                }
            }
        }
    }
}

void CnvettiNormalizeApp::processRegion(std::string const & contig, int beginPos, int endPos)
{
    if (options.verbosity >= 1)
    {
        if (beginPos != 0 && endPos != 0)
            std::cerr << "Processing " << contig << ":" << (beginPos + 1) << "-" << endPos << "\n";
        else
            std::cerr << "Processing " << contig << "\n";
    }

    if (bcf_sr_seek(vcfReaderPtr.get(), contig.c_str(), 0) != 0)
    {
        std::cerr << "WARNING: Could not seek to beginning of " << contig << "\n";
        return;
    }

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (contig != vcfReaderPtr.get()->regions->seq_names[line->rid])
            break;  // on different contig
        if (bcf_translate(vcfHeaderOut, vcfHeaderIn, line))
            throw std::runtime_error("Could not translate from old to new header");

        float * ptrGC = nullptr;
        int32_t gcContent;
        int32_t sizeGC = 0;
        if (bcf_get_info_float(vcfHeaderOut, line, "GC", &ptrGC, &sizeGC) == 1)
        {
            gcContent = *ptrGC;
            free(ptrGC);
        }
        else
        {
            free(ptrGC);
            throw std::runtime_error("Could not get INFO/GC value");
        }

        int32_t windowCount = medianStats[0][gcContent].windowCount;
        bcf_update_info_int32(vcfHeaderOut, line, "GCWINDOWS", &windowCount, 1);
        if (windowCount < options.minGCWindowCount)
        {
            int * filters = (int *)calloc(1, sizeof(int));
            *filters = bcf_hdr_id2int(vcfHeaderOut, BCF_DT_ID, "FEW_GCWINDOWS");
            bcf_update_filter(vcfHeaderOut, line, filters, 1);
            free(filters);
        }

        // Normalize COV
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderOut, line, "COV", &values, &sizeVals);
            if (!res)
            {
                free(values);
                throw std::runtime_error("Could not get FORMAT/COV value");
            }
            for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
            {
                if (medianStats[sampleID][gcContent].cov)
                    values[sampleID] /= medianStats[sampleID][gcContent].cov;
                else
                    values[sampleID] = 0;

                madHistos[sampleID][abs(MAD_RES - round(MAD_RES * values[sampleID]))].cov += 1;
            }
            bcf_update_format_float(vcfHeaderOut, line, "COV", values, sizeVals);
            free(values);
        }

        // Normalize COV0
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderOut, line, "COV0", &values, &sizeVals);
            if (!res)
            {
                free(values);
                throw std::runtime_error("Could not get FORMAT/COV0 value");
            }
            for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
            {
                if (medianStats[sampleID][gcContent].cov0)
                    values[sampleID] /= medianStats[sampleID][gcContent].cov0;
                else
                    values[sampleID] = 0;

                madHistos[sampleID][abs(MAD_RES - round(MAD_RES * values[sampleID]))].cov0 += 1;
            }
            bcf_update_format_float(vcfHeaderOut, line, "COV0", values, sizeVals);
            free(values);
        }

        // Normalize WINSD
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderOut, line, "WINSD", &values, &sizeVals);
            if (!res)
            {
                free(values);
                throw std::runtime_error("Could not get FORMAT/WINSD value");
            }
            for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
                if (medianStats[sampleID][gcContent].cov)
                    values[sampleID] /= medianStats[sampleID][gcContent].cov;
                else
                    values[sampleID] = 0;
            bcf_update_format_float(vcfHeaderOut, line, "WINSD", values, sizeVals);
            free(values);
        }

        // Normalize WINSD0
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderOut, line, "WINSD0", &values, &sizeVals);
            if (!res)
            {
                free(values);
                throw std::runtime_error("Could not get FORMAT/WINSD0 value");
            }
            for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
                if (medianStats[sampleID][gcContent].cov0)
                    values[sampleID] /= medianStats[sampleID][gcContent].cov0;
                else
                    values[sampleID] = 0;
            bcf_update_format_float(vcfHeaderOut, line, "WINSD0", values, sizeVals);
            free(values);
        }

        // Normalize RC
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderOut, line, "RC", &values, &sizeVals);
            if (!res)
            {
                free(values);
                throw std::runtime_error("Could not get FORMAT/RC value");
            }
            for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
            {
                if (medianStats[sampleID][gcContent].rc)
                    values[sampleID] /= medianStats[sampleID][gcContent].rc;
                else
                    values[sampleID] = 0;

                madHistos[sampleID][abs(MAD_RES - round(MAD_RES * values[sampleID]))].rc += 1;
            }
            bcf_update_format_float(vcfHeaderOut, line, "RC", values, sizeVals);
            free(values);
        }

        // Normalize RC0
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderOut, line, "RC0", &values, &sizeVals);
            if (!res)
            {
                free(values);
                throw std::runtime_error("Could not get FORMAT/RC0 value");
            }
            for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
            {
                if (medianStats[sampleID][gcContent].rc0)
                    values[sampleID] /= medianStats[sampleID][gcContent].rc0;
                else
                    values[sampleID] = 0;

                madHistos[sampleID][abs(MAD_RES - round(MAD_RES * values[sampleID]))].rc0 += 1;
            }
            bcf_update_format_float(vcfHeaderOut, line, "RC0", values, sizeVals);
            free(values);
        }

        bcf_write(vcfFileOutPtr.get(), vcfHeaderOut, line);
    }
}

double medianValue(std::map<int, int> const & histo)
{
    // Determine total number of value registered in histogram
    int totalCount = 0;
    for (auto const & pair : histo)
        totalCount += pair.second;
    if (totalCount == 0u)
        return 0.0;  // Short-circuit if empty

    // Compute median
    int low;
    int up;
    if (totalCount % 2 == 0)
    {
        low = (totalCount / 2) - 1;
        up = low + 1;
    }
    else
    {
        low = totalCount / 2;
        up = low;
    }

    int a = -1, b = -1;
    for (auto const & pair : histo)
    {
        int const val = pair.first;
        int const num = pair.second;
        if (a == -1 && num > low)
            a = val;
        else if (a == -1 && num <= low)
            low -= num;
        if (b == -1 && num > up)
            b = val;
        else if (b == -1 && num <= up)
            up -= num;
        if (a != -1 && b != -1)
            break;
    }

    if (a == -1 && b == -1)
        return 0.0;
    else
        return static_cast<double>(a + b) / 2.0;
}

void CnvettiNormalizeApp::printStatistics(std::ostream & out)
{
    std::cerr << "\nSTATISTICS\n\n";

    for (int sampleID = 0; sampleID < bcf_hdr_nsamples(vcfHeaderOut); ++sampleID)
    {
        std::cerr << "SAMPLE: " << vcfHeaderOut->samples[sampleID] << "\n\n";

        std::map<int, MadHistoEntry> const & madHisto = madHistos[sampleID];

        // std::cerr << "\nAbsolute deviations coverage\n\n";
        std::map<int, int> covHisto;
        for (auto const & entry : madHisto)
        {
            // std::cerr << entry.first / MAD_RES << "\t" << entry.second.cov << "\n";
            covHisto[entry.first] = entry.second.cov;
        }
        std::cerr << "\nMAD coverage: " << medianValue(covHisto) / MAD_RES * 100 << "%\n";

        // std::cerr << "\nAbsolute deviations coverage (incl. q0)\n\n";
        std::map<int, int> cov0Histo;
        for (auto const & entry : madHisto)
        {
            // std::cerr << entry.first / MAD_RES << "\t" << entry.second.cov0 << "\n";
            cov0Histo[entry.first] = entry.second.cov0;
        }
        std::cerr << "\nMAD coverage (incl. q0): " << medianValue(cov0Histo) / MAD_RES * 100 << "%\n";

        // std::cerr << "\nAbsolute deviations read count\n\n";
        std::map<int, int> rcHisto;
        for (auto const & entry : madHisto)
        {
            // std::cerr << entry.first / MAD_RES << "\t" << entry.second.rc << "\n";
            rcHisto[entry.first] = entry.second.rc;
        }
        std::cerr << "\nMAD read count: " << medianValue(rcHisto) / MAD_RES * 100 << "%\n";

        // std::cerr << "\nAbsolute deviations read count (incl. q0)\n\n";
        std::map<int, int> rc0Histo;
        for (auto const & entry : madHisto)
        {
            // std::cerr << entry.first / MAD_RES << "\t" << entry.second.rc0 << "\n";
            rc0Histo[entry.first] = entry.second.rc0;
        }
        std::cerr << "\nMAD read count (incl. q0): " << medianValue(rc0Histo) / MAD_RES * 100 << "%\n";
    }
}

}  // anonymous namespace

int mainNormalize(CnvettiNormalizeOptions const & options)
{
    CnvettiNormalizeApp app(options);
    app.run();

    return 0;
}
