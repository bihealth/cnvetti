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

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>

#include "cnvetti/program_options.h"
#include "cnvetti/version.h"
#include "cnvetti/utils.h"

namespace {  // anonymous namespace

// ---------------------------------------------------------------------------
// Cosntant STATS_NAMES
// ---------------------------------------------------------------------------

// Names of the statistics fields to compute
std::vector<std::string> const STATS_NAMES { "COV", "COV0", "RC", "RC0" };

// ---------------------------------------------------------------------------
// Class CnvettiBackgroundApp
// ---------------------------------------------------------------------------

class CnvettiBackgroundApp
{
public:
    CnvettiBackgroundApp(CnvettiBackgroundOptions const & options) :
            options(options), vcfHeaderIn(nullptr), vcfHeaderOut(nullptr)
    {}

    void run();

private:

    // Open the input files and store them in bamFilesIn and open vcfFileOutPtr;
    void openFiles();
    // Check that the opened files look correct
    void checkFiles();
    // Perform chromosome-wise processing
    void processChromosomes();
    // Process a region
    void processRegion(std::string const & contig, int beginPos, int endPos);

    CnvettiBackgroundOptions const & options;

    std::shared_ptr<bcf_srs_t> vcfReaderPtr;
    std::shared_ptr<vcfFile> vcfFileOutPtr;
    std::shared_ptr<bcf_hdr_t> vcfHeaderOutPtr;
    bcf_hdr_t * vcfHeaderIn;
    bcf_hdr_t * vcfHeaderOut;
};


void CnvettiBackgroundApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti background\n"
                  << "=================\n\n";
        options.print(std::cerr);
    }

    if (options.verbosity >= 1)
        std::cerr << "\nOPENING FILES\n\n";
    openFiles();

    if (options.verbosity >= 1)
        std::cerr << "\nPROCESSING\n\n";
    processChromosomes();

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiBackgroundApp::openFiles()
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

    // Copy input to output header, maybe removing samples
    if (options.writeSamples)
        vcfHeaderOutPtr = std::shared_ptr<bcf_hdr_t>(bcf_hdr_dup(vcfHeaderIn), bcf_hdr_destroy);
    else
        vcfHeaderOutPtr = std::shared_ptr<bcf_hdr_t>(
            bcf_hdr_subset(vcfHeaderIn, 0, NULL, NULL), bcf_hdr_destroy);
    vcfHeaderOut = vcfHeaderOutPtr.get();

    for (std::string const & stat : STATS_NAMES)
    {
        bcf_hdr_printf(vcfHeaderOut,
            ("##INFO=<ID=SUMMARY_%s,Type=Float,Number=5,Description="
             "\"Five-number summary of the %s data (Q0, Q1, Q2, Q3, Q4)\">"),
            stat.c_str(), stat.c_str());
        bcf_hdr_printf(vcfHeaderOut,
            ("##INFO=<ID=IQR_%s,Type=Float,Number=1,Description="
             "\"Inter-quartile range of the %s data\">"),
            stat.c_str(), stat.c_str());
    }

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeaderOut, "##cnvetti_backgroundCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOutPtr.get(), vcfHeaderOut);
}

void CnvettiBackgroundApp::checkFiles()
{
    // Check that the file looks like created by "cnvetti coverage"
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

void CnvettiBackgroundApp::processChromosomes()
{
    for (int i = 0; i < vcfHeaderIn->nhrec; ++i)  // skip pseudo contig
    {
        bcf_hrec_t * hrec = vcfHeaderIn->hrec[i];
        if (hrec->type == BCF_HL_CTG)
        {
            for (int j = 0; j < hrec->nkeys; ++j)
            {
                if (strcmp("ID", hrec->keys[j]) == 0)
                    processRegion(hrec->vals[j], 0, 0);
            }
        }
    }
}

void CnvettiBackgroundApp::processRegion(std::string const & contig, int beginPos, int endPos)
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
        std::cerr << "WARNING: could not seek to beginning of region " << contig << "\n";
        return;
    }

    // Five-number summaries for each statistics value
    std::map<std::string, std::vector<float>> summaries;
    for (std::string const & name : STATS_NAMES)
        summaries[name].resize(5, 0.0);

    std::vector<float> vals;
    float * valsPtr = nullptr;
    int32_t sizeVals = 0;
    int res;

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);

        for (std::string const & name : STATS_NAMES)
        {
            res = bcf_get_format_float(vcfHeaderIn, line, name.c_str(), &valsPtr, &sizeVals);
            if (!res)
            {
                free(valsPtr);
                throw std::runtime_error(std::string("Could not get value of FORMAT/") + name);
            }

            double index, h, lo, hi;

            std::sort(valsPtr, valsPtr + sizeVals);
            vals.clear();
            std::copy(valsPtr, valsPtr + sizeVals, std::back_inserter(vals));

            // Q0
            summaries[name][0] = valsPtr[0];

            // Q1
            index = (sizeVals - 1) * 0.25;
            lo = floor(index);
            hi = ceil(index);
            h = index - lo;
            summaries[name][1] = (1 - h) * vals.at(size_t(lo)) + h * vals.at(size_t(hi));

            // Q2
            if (sizeVals % 2 == 1)
                summaries[name][2] = vals.at(sizeVals / 2);
            else
                summaries[name][2] = (vals.at(sizeVals / 2) + vals.at((sizeVals + 1) / 2)) / 2.0;

            // Q3
            index = (sizeVals - 1) * 0.75;
            lo = floor(index);
            hi = ceil(index);
            h = index - lo;
            summaries[name][3] = (1 - h) * vals.at(size_t(lo)) + h * vals.at(size_t(hi));

            // Q4
            summaries[name][4] = vals.at(sizeVals - 1);
        }

        if (!options.writeSamples)
            if (bcf_subset(vcfHeaderIn, line, 0, NULL))
                throw std::runtime_error("Could not remove samples from BCF record");
        if (bcf_translate(vcfHeaderOut, vcfHeaderIn, line))
           throw std::runtime_error("Could not translate from old to new header");

        for (std::string const & name : STATS_NAMES)
        {
            std::string const keySummary = "SUMMARY_" + name;
            bcf_update_info_float(vcfHeaderOut, line, keySummary.c_str(), &summaries[name][0], 5);

            std::string const keyIQR = "IQR_" + name;
            float const iqr = summaries[name][3] - summaries[name][1];
            bcf_update_info_float(vcfHeaderOut, line, keyIQR.c_str(), &iqr, 1);
        }

        bcf_write(vcfFileOutPtr.get(), vcfHeaderOut, line);
    }

    free(valsPtr);
}

}  // anonymous namespace

int mainBackground(CnvettiBackgroundOptions const & options)
{
    CnvettiBackgroundApp app(options);
    app.run();

    return 0;
}
