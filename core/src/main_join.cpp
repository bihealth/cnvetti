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

// ---------------------------------------------------------------------------
// Struct MedianStats
// ---------------------------------------------------------------------------

// TODO: rename to statsCounters

struct MedianStats
{
    double cov { 0 };
    double cov0 { 0 };
    double rc { 0 };
    double rc0 { 0 };
};

// ---------------------------------------------------------------------------
// Class CnvettiJoinApp
// ---------------------------------------------------------------------------

class CnvettiJoinApp
{
public:
    CnvettiJoinApp(CnvettiJoinOptions const & options) :
            options(options), totalLength(0), vcfHeaderIn(nullptr), vcfHeaderOut(nullptr)
    {}

    void run();

private:

    // Open the input files and store them in bamFilesIn and open vcfFileOutPtr;
    void openFiles();
    // Check that the references FASTA file fits to the BAM file.
    void checkFiles();
    // Perform chromosome-wise processing
    void processChromosomes();
    // Process a region
    void processRegion(std::string const & contig, int beginPos, int endPos);
    // Write out medians
    void writeMedians();

    // Handle writing out of the current window.
    void _handleBetweenWindows(
        int & rID,
        int & prevBeginPos,
        int & prevEndPos,
        std::vector<double> & cov,
        std::vector<double> & cov0,
        std::vector<double> & rc,
        std::vector<double> & rc0,
        bcf1_t * line
    );

    CnvettiJoinOptions const & options;

    // medianStats[sampleID]
    std::vector<MedianStats> medianStats;
    // total window length sum
    int32_t totalLength;

    std::shared_ptr<bcf_srs_t> vcfReaderPtr;
    std::shared_ptr<vcfFile> vcfFileOutPtr;
    std::shared_ptr<bcf_hdr_t> vcfHeaderOutPtr;
    bcf_hdr_t * vcfHeaderIn;
    bcf_hdr_t * vcfHeaderOut;
};


void CnvettiJoinApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti join\n"
                  << "============\n\n";
        options.print(std::cerr);
    }

    if (options.verbosity >= 1)
        std::cerr << "\nOPENING FILES\n\n";
    openFiles();

    if (options.verbosity >= 1)
        std::cerr << "\nPROCESSING\n\n";
    processChromosomes();

    if (options.verbosity >= 1)
        std::cerr << "\nWRITING MEDIANS\n\n";
    writeMedians();

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiJoinApp::openFiles()
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

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeaderOut, "##cnvetti_joinCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOutPtr.get(), vcfHeaderOut);
}

void CnvettiJoinApp::checkFiles()
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

void CnvettiJoinApp::processChromosomes()
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

void CnvettiJoinApp::_handleBetweenWindows(
    int & rID,
    int & prevBeginPos,
    int & prevEndPos,
    std::vector<double> & cov,
    std::vector<double> & cov0,
    std::vector<double> & rc,
    std::vector<double> & rc0,
    bcf1_t * line
) {
    // TODO: Fugly, needs to be kept in sync with final window processing at the moment
    if (prevBeginPos == -1)
    {
        if (line)
        {
            prevBeginPos = line->pos;
        }
    }
    else if (prevBeginPos != -1 && (!line || line->pos != prevEndPos))
    {
        std::shared_ptr<bcf1_t> out(bcf_init(), &bcf_destroy);

        // CHROM, POS, REF, ALT
        out->rid = rID;
        out->pos = prevBeginPos;
        bcf_update_alleles_str(vcfHeaderOut, out.get(), "N,<COUNT>");

        int const windowLength = prevEndPos - prevBeginPos;

        // Set INFO/GC
        float gc = 0.0;
        bcf_update_info_float(vcfHeaderOut, out.get(), "GC", &gc, 1);

        // Set INFO/END
        std::vector<int32_t> end;
        end.resize(1, prevEndPos);
        int32_t sizeEnd = 1;
        if (bcf_update_info_int32(vcfHeaderOut, out.get(), "END", &end[0], sizeEnd) != 0)
        {
            throw std::runtime_error("Could not set INFO/END");
        }

        // GT
        std::vector<int32_t> gts(2 * cov.size(), bcf_gt_missing);
        bcf_update_genotypes(vcfHeaderOut, out.get(), &gts[0], gts.size());

        // Set FORMAT/COV
        std::vector<float> floatCov(cov.begin(), cov.end());
        for (float & val : floatCov)
        {
            val /= windowLength;
        }
        int32_t sizeCov = floatCov.size();
        if (bcf_update_format_float(vcfHeaderOut, out.get(), "COV", &floatCov[0], sizeCov) != 0)
        {
            throw std::runtime_error("Could not set FORMAT/COV value!");
        }

        // Set FORMAT/COV0
        std::vector<float> floatCov0(cov0.begin(), cov0.end());
        for (float & val : floatCov0)
        {
            val /= windowLength;
        }
        int32_t sizeCov0 = floatCov0.size();
        if (bcf_update_format_float(vcfHeaderOut, out.get(), "COV0", &floatCov0[0], sizeCov0) != 0)
        {
            throw std::runtime_error("Could not set FORMAT/COV0 value!");
        }

        // Set FORMAT/RC
        std::vector<float> floatRc(rc.begin(), rc.end());
        for (float & val : floatRc)
        {
            val /= windowLength;
        }
        int32_t sizeRc = floatRc.size();
        if (bcf_update_format_float(vcfHeaderOut, out.get(), "RC", &floatRc[0], sizeRc) != 0)
        {
            throw std::runtime_error("Could not set FORMAT/RC value!");
        }

        // Set FORMAT/RC0
        std::vector<float> floatRc0(rc0.begin(), rc0.end());
        for (float & val : floatRc0)
        {
            val /= windowLength;
        }
        int32_t sizeRc0 = floatRc0.size();
        if (bcf_update_format_float(vcfHeaderOut, out.get(), "RC0", &floatRc0[0], sizeRc0) != 0)
        {
            throw std::runtime_error("Could not set FORMAT/RC0 value!");
        }

        if (bcf_write(vcfFileOutPtr.get(), vcfHeaderOut, out.get()) != 0)
        {
            throw std::runtime_error("Could not write to output file");
        }

        // Reset counters and prev*
        std::fill(cov.begin(), cov.end(), 0.0);
        std::fill(cov0.begin(), cov0.end(), 0.0);
        std::fill(rc.begin(), rc.end(), 0.0);
        std::fill(rc0.begin(), rc0.end(), 0.0);

        if (line)
        {
            prevBeginPos = line->pos;
        }
        // prevEndPos will be updated at bottom of loop
    }
}

void CnvettiJoinApp::processRegion(std::string const & contig, int rBeginPos, int rEndPos)
{
    if (options.verbosity >= 1)
    {
        if (rBeginPos != 0 && rEndPos != 0)
            std::cerr << "Processing " << contig << ":" << (rBeginPos + 1) << "-" << rEndPos << "\n";
        else
            std::cerr << "Processing " << contig << "\n";
    }

    if (bcf_sr_seek(vcfReaderPtr.get(), contig.c_str(), 0) != 0)
    {
        std::cerr << "WARNING: Could not seek to beginning of " << contig << "\n";
        return;
    }

    // Increase stats counters if necessary.
    medianStats.resize(bcf_hdr_nsamples(vcfHeaderOut));

    // Values kept for building the summarized statistics.
    int rID = -1;
    int prevBeginPos = -1;
    int prevEndPos = 0;
    // Sum of the values, to be divided by totalLength.
    std::vector<double> cov(bcf_hdr_nsamples(vcfHeaderOut), 0.0);
    std::vector<double> cov0(bcf_hdr_nsamples(vcfHeaderOut), 0.0);
    std::vector<double> rc(bcf_hdr_nsamples(vcfHeaderOut), 0.0);
    std::vector<double> rc0(bcf_hdr_nsamples(vcfHeaderOut), 0.0);

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (contig != vcfReaderPtr.get()->regions->seq_names[line->rid])
            break;  // on different contig
        if (bcf_translate(vcfHeaderOut, vcfHeaderIn, line))
            throw std::runtime_error("Could not translate from old to new header");

        // Get Contig ID.
        rID = line->rid;

        _handleBetweenWindows(
            rID, prevBeginPos, prevEndPos, cov, cov0, rc, rc0, line);

        // Get begin position.
        int beginPos = line->pos;

        // Get end position
        int32_t * ptrEnd = nullptr;
        int32_t endPos;
        int32_t sizeEnd = 0;
        if (bcf_get_info_int32(vcfHeaderOut, line, "END", &ptrEnd, &sizeEnd) == 1)
        {
            endPos = *ptrEnd;
            free(ptrEnd);
        }
        else
        {
            free(ptrEnd);
            throw std::runtime_error("Could not get INFO/END value");
        }


        // Increase window length.
        int const windowLength = (endPos - beginPos);
        totalLength += windowLength;

        // Accumulate COV
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
                cov[sampleID] += values[sampleID] * windowLength;
                medianStats[sampleID].cov += values[sampleID] * windowLength;
            }
            free(values);
        }

        // Accumulate COV0
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
                cov0[sampleID] += values[sampleID] * windowLength;
                medianStats[sampleID].cov0 += values[sampleID] * windowLength;
            }
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
                rc[sampleID] += values[sampleID];
                medianStats[sampleID].rc += values[sampleID];
            }
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
                rc0[sampleID] += values[sampleID];
                medianStats[sampleID].rc0 += values[sampleID];
            }
            free(values);
        }

        prevEndPos = endPos;
    }

    _handleBetweenWindows(
        rID, prevBeginPos, prevEndPos, cov, cov0, rc, rc0, nullptr);
}

void CnvettiJoinApp::writeMedians()
{
    bcf_hdr_t * const hdr = vcfHeaderOutPtr.get();
    int rID = 0;
    for (int i = 0; i < hdr->nhrec; ++i)  // skip pseudo contig
    {
        bcf_hrec_t * hrec = hdr->hrec[i];
        if (hrec->type == BCF_HL_CTG)
        {
            for (int j = 0; j < hrec->nkeys; ++j)
            {
                if (strcmp("ID", hrec->keys[j]) == 0)
                {
                    if (strcmp("__cnvetti_stats", hrec->vals[j]) == 0)
                    {
                        std::cerr << "STATS FOUND " << rID << "\n";
                        break;
                    }
                    else
                    {
                        rID += 1;
                    }
                }
            }
        }
    }

    std::unique_ptr<bcf1_t, void (&)(bcf1_t*)> recordPtr(bcf_init(), bcf_destroy);
    bcf1_t * record = recordPtr.get();

    // CHROM, POS, REF, ALT, FILTER
    record->rid = rID;
    record->pos = 0;
    bcf_update_alleles_str(hdr, record, "N,<STATS>");

    std::unique_ptr<int, std::function<void(int*)>> filters(
        (int*)malloc(sizeof(int)), [](int * ptr) { free(ptr); });
    *filters = bcf_hdr_id2int(hdr, BCF_DT_ID, "CNVETTI_STATS");
    bcf_update_filter(hdr, record, filters.get(), 1);

    // INFO fields
    float valGCContent = 0.0f;
    bcf_update_info_float(hdr, record, "GC", &valGCContent, 1);
    int32_t valGCWindows = 1234;
    bcf_update_info_int32(hdr, record, "GCWINDOWS", &valGCWindows, 1);
    int32_t end = 100;
    bcf_update_info_int32(hdr, record, "END", &end, 1);

    int numSamples = medianStats.size();

    // GT
    std::vector<int32_t> gts(2 * numSamples, bcf_gt_missing);
    bcf_update_genotypes(hdr, record, &gts[0], gts.size());

    // COV
    std::vector<float> cov(numSamples, 0.0);
    for (unsigned i = 0; i < numSamples; ++i)
        cov[i] = medianStats[i].cov / totalLength;
    bcf_update_format_float(hdr, record, "COV", &cov[0], cov.size());

    // COV0
    std::vector<float> cov0(numSamples, 0.0);
    for (unsigned i = 0; i < numSamples; ++i)
        cov0[i] = medianStats[i].cov0 / totalLength;
    bcf_update_format_float(hdr, record, "COV0", &cov0[0], cov0.size());

    // RC
    std::vector<float> rc(numSamples, 0);
    for (unsigned i = 0; i < numSamples; ++i)
        rc[i] = medianStats[i].rc / totalLength * 100.0;
    bcf_update_format_float(hdr, record, "RC", &rc[0], rc.size());

    // RC0
    std::vector<float> rc0(numSamples, 0);
    for (unsigned i = 0; i < numSamples; ++i)
        rc0[i] = medianStats[i].rc0 / totalLength * 100.0;
    bcf_update_format_float(hdr, record, "RC0", &rc0[0], rc.size());

    bcf_write(vcfFileOutPtr.get(), hdr, record);
}

}  // anonymous namespace

int mainJoin(CnvettiJoinOptions const & options)
{
    CnvettiJoinApp app(options);
    app.run();

    return 0;
}
