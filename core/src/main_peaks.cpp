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
#include <sstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/synced_bcf_reader.h>

#include "cnvetti/program_options.h"
#include "cnvetti/version.h"
#include "cnvetti/utils.h"

namespace {  // anonymous namespace

// ---------------------------------------------------------------------------
// Struct HistogramEntry
// ---------------------------------------------------------------------------

struct HistogramEntry
{
    uint64_t cov0 { 0 };
};

// ---------------------------------------------------------------------------
// Class CnvettiPeaksApp
// ---------------------------------------------------------------------------

class CnvettiPeaksApp
{
public:
    CnvettiPeaksApp(CnvettiPeaksOptions const & options) :
            options(options), vcfHeaderIn(nullptr)
    {}

    void run();

private:

    // Open the input files and store them in bamFilesIn.
    void openFiles();
    // Perform chromosome-wise processing
    void processChromosomes(
            std::function<void(std::string const &, int, int)> func);
    // Augment histogram for the given region.
    void processRegionCreateHistogram(
            std::string const & contig, int beginPos, int endPos);
    // Write out BED file for the given region.
    void processRegionWriteBed(
            std::string const & contig, int beginPos, int endPos, double threshold);

    CnvettiPeaksOptions const & options;

    // histogram[percentile]
    std::map<int, HistogramEntry> histogram;

    std::shared_ptr<bcf_srs_t> vcfReaderPtr;
    bcf_hdr_t * vcfHeaderIn;
    std::unique_ptr<BGZF, std::function<void(BGZF*)>> outBgzf;

    uint64_t peakCount { 0 };
};


void CnvettiPeaksApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti peaks\n"
                  << "=============\n\n";
        options.print(std::cerr);
    }

    if (options.verbosity >= 1)
        std::cerr << "\nOPENING FILES\n\n";
    openFiles();

    if (options.verbosity >= 1)
        std::cerr << "\nCOMPUTING HISTOGRAM\n\n";
    processChromosomes(
            [&](std::string const & chrom, int beginPos, int endPos) {
                processRegionCreateHistogram(chrom, beginPos, endPos);
            });

    if (options.verbosity >= 1)
        std::cerr << "\nDETERMINING THRESHOLD\n\n";
    uint64_t totalCov0 = 0;
    for (auto entry : histogram)
    {
        if (options.verbosity >= 2)
            std::cerr << (entry.first / 100.0) << "\t" << entry.second.cov0 << "\n";
        totalCov0 += entry.second.cov0;
    }
    if (options.verbosity >= 1)
        std::cerr << "\n\n";
    uint64_t val = totalCov0 * (options.percentile / 100.0);
    double threshold = 0.0;
    for (auto entry : histogram)
    {
        threshold = entry.first / 100.0;
        if (val >= entry.second.cov0)
        {
            val -= entry.second.cov0;
        }
        else
        {
            break;
        }
    }
    threshold *= options.threshFactor;
    if (options.verbosity >= 1)
        std::cerr << "\n\n=> threshold == " << threshold << "\n";

    if (options.verbosity >= 1)
        std::cerr << "\nGENERATING REGIONS BED\n\n";
    processChromosomes(
            [&](std::string const & chrom, int beginPos, int endPos) {
                processRegionWriteBed(chrom, beginPos, endPos, threshold);
            });

    if (options.verbosity >= 1)
        std::cerr << "peak count: " << peakCount << "\n";

    if (options.verbosity >= 1)
        std::cerr << "\nCLOSE BGZF, CREATE INDEX\n\n";
    outBgzf.reset();
    tbx_index_build(options.outputFileName.c_str(), -1, &tbx_conf_bed);

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiPeaksApp::openFiles()
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

    // Open output file.
    if (options.verbosity >= 1)
        std::cerr << "Opening " << options.outputFileName << " ...\n";
    outBgzf = std::unique_ptr<BGZF, std::function<void(BGZF*)>>(
            bgzf_open(options.outputFileName.c_str(), "w"), [](BGZF * bgzf) { bgzf_close(bgzf); } );
}

void CnvettiPeaksApp::processChromosomes(
        std::function<void(std::string const &, int, int)> func)
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
                    func(hrec->vals[j], 0, 0);
                }
            }
        }
    }
}

void CnvettiPeaksApp::processRegionCreateHistogram(
        std::string const & contig, int beginPos, int endPos)
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

    // Count samples.
    int const numSamples = bcf_hdr_nsamples(vcfHeaderIn);

    // Compute histogram in bins in resolution of 0.01x.
    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (contig != vcfReaderPtr.get()->regions->seq_names[line->rid])
            break;  // on different contig

        float * tmpCov0s = nullptr;
        int32_t sizeCov0s = 0;
        int res = bcf_get_format_float(vcfHeaderIn, line, "COV0", &tmpCov0s, &sizeCov0s);
        typedef std::unique_ptr<float, std::function<void(float*)>> TCleanedPtr;
        TCleanedPtr cov0s = TCleanedPtr(tmpCov0s, [](float * data) { free(data); });
        if (!res)
            throw std::runtime_error("Could not get FORMAT/COV0 value");

        for (int sampleID = 0; sampleID < numSamples; ++sampleID)
        {
            ++histogram[(int)(std::ceil(cov0s.get()[sampleID] * 100))].cov0;
        }
    }
}

void CnvettiPeaksApp::processRegionWriteBed(
        std::string const & contig, int beginPos, int endPos, double threshold)
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

    // Count samples.
    int const numSamples = bcf_hdr_nsamples(vcfHeaderIn);
    std::stringstream ss;

    int chunkBegin = 0;
    int chunkEnd = 0;

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (contig != vcfReaderPtr.get()->regions->seq_names[line->rid])
            break;  // on different contig

        float * tmpCov0s = nullptr;
        int32_t sizeCov0s = 0;
        int res = bcf_get_format_float(vcfHeaderIn, line, "COV0", &tmpCov0s, &sizeCov0s);
        typedef std::unique_ptr<float, std::function<void(float*)>> TCleanedPtr;
        TCleanedPtr cov0s = TCleanedPtr(tmpCov0s, [](float * data) { free(data); });
        if (!res)
            throw std::runtime_error("Could not get FORMAT/COV0 value");

        int32_t end[] = { 0 };
        int32_t * endPtr = &end[0];
        int32_t numOut = 1;
        if (!bcf_get_info_int32(vcfHeaderIn, line, "END", &endPtr, &numOut))
            throw std::runtime_error("Could not get INFO/END value.");

        bool anyPeak = false;
        for (int sampleID = 0; sampleID < numSamples; ++sampleID)
        {
            if (cov0s.get()[sampleID] >= threshold)
            {
                anyPeak = true;
                break;
            }
        }

        if (anyPeak)
        {
            peakCount += 1;

            if (chunkBegin == chunkEnd)  // first chunk
            {
                chunkBegin = line->pos;
                chunkEnd = end[0];
            }
            else
            {
                if (chunkEnd == line->pos)  // extend
                {
                    chunkEnd = end[0];
                }
                else  // write out and create new chunk
                {
                    ss.str("");
                    ss.clear();
                    ss << contig << "\t" << chunkBegin << "\t" << chunkEnd << "\n";
                    if (bgzf_write(outBgzf.get(), ss.str().c_str(), ss.str().size()) != ss.str().size())
                    {
                        throw std::runtime_error("Problem writing to BGZF BED file.");
                    }

                    chunkBegin = line->pos;
                    chunkEnd = end[0];
                }
            }
        }
    }

    if (chunkBegin != chunkEnd)
    {
        ss.str("");
        ss.clear();
        ss << contig << "\t" << chunkBegin << "\t" << chunkEnd << "\n";
        if (bgzf_write(outBgzf.get(), ss.str().c_str(), ss.str().size()) != ss.str().size())
        {
            throw std::runtime_error("Problem writing to BGZF BED file.");
        }
    }
}

}  // anonymous namespace

int mainPeaks(CnvettiPeaksOptions const & options)
{
    CnvettiPeaksApp app(options);
    app.run();

    return 0;
}
