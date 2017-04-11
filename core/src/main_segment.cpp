// ============================================================================
//                                 CNVetti
// ============================================================================
// Copyright 2016-2017 Berlin Institute for Health
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this softwareand associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// ============================================================================
// Author:  Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>
// ============================================================================

#include <algorithm>
#include <cmath>
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
#include "cnvetti/segmentation.h"

namespace {  // anonymous namespace

// ---------------------------------------------------------------------------
// Class CnvettiSegmentApp
// ---------------------------------------------------------------------------

class CnvettiSegmentApp
{
public:
    CnvettiSegmentApp(CnvettiSegmentOptions const & options) :
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
    // Perform segmentation for one sample
    void doSegmentation(std::vector<int> const & pos, std::vector<float> & values) const;

    CnvettiSegmentOptions const & options;

    std::shared_ptr<bcf_srs_t> vcfReaderPtr;
    std::shared_ptr<vcfFile> vcfFileOutPtr;
    std::shared_ptr<bcf_hdr_t> vcfHeaderOutPtr;
    bcf_hdr_t * vcfHeaderIn;
    bcf_hdr_t * vcfHeaderOut;
};


void CnvettiSegmentApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti segment\n"
                  << "===============\n\n";
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

void CnvettiSegmentApp::openFiles()
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
    vcfHeaderOutPtr = std::shared_ptr<bcf_hdr_t>(bcf_hdr_dup(vcfHeaderIn), bcf_hdr_destroy);
    vcfHeaderOut = vcfHeaderOutPtr.get();

    bcf_hdr_printf(vcfHeaderOut,
        ("##INFO=<ID=SVTYPE,Type=String,Number=1,Description="
         "\"Type of called structural variant\">"));
    bcf_hdr_printf(vcfHeaderOut,
        ("##INFO=<ID=CIPOS,Type=Integer,Number=2,Description="
         "\"Confidence interval around POS\">"));
    bcf_hdr_printf(vcfHeaderOut,
        ("##INFO=<ID=CIEND,Type=Integer,Number=2,Description="
         "\"Confidence interval around END\">"));
    bcf_hdr_printf(vcfHeaderOut,
        ("##ALT=<ID=CNV,Type=String,Number=1,Description="
         "\"Value for ALT field for describing CNV called by CNVetti\">"));
    bcf_hdr_printf(vcfHeaderOut,
        ("##FORMAT=<ID=CN,Type=Integer,Number=1,Description="
         "\"Discrete copy number data\">"));
    bcf_hdr_printf(vcfHeaderOut,
        ("##FORMAT=<ID=BC,Type=Integer,Number=1,Description="
         "\"Base count for the CNV call\">"));
    bcf_hdr_printf(vcfHeaderOut,
        ("##FORMAT=<ID=RC,Type=Integer,Number=1,Description="
         "\"Read count for the CNV call\">"));

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeaderOut, "##cnvetti_segmentCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOutPtr.get(), vcfHeaderOut);
}

void CnvettiSegmentApp::checkFiles()
{
    // Check that the file looks like created by "cnvetti background"
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

void CnvettiSegmentApp::processChromosomes()
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

// TODO: this needs a lot of cleanup!
void CnvettiSegmentApp::processRegion(std::string const & contig, int beginPos, int endPos)
{
    if (options.verbosity >= 1)
    {
        if (beginPos != 0 && endPos != 0)
            std::cerr << "Processing " << contig << ":" << (beginPos + 1) << "-" << endPos << "\n";
        else
            std::cerr << "Processing " << contig << "\n";
    }

    if (bcf_sr_seek(vcfReaderPtr.get(), contig.c_str(), 0) != 0)
        throw std::runtime_error("Could not seek to beginning of region!");

    std::vector<int> pos;
    std::vector<std::vector<float>> values;
    values.resize(bcf_hdr_nsamples(vcfHeaderIn));

    float * valsPtr = nullptr;
    int32_t sizeVals = 0;

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);

        // Weed out the noisy windows
        float gcContent, mapability, iqrRc0, iqrCov0;

        float * ptrTmp = nullptr;
        int32_t sizeTmp = 0;
        if (bcf_get_info_float(vcfHeaderIn, line, "GC", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/GC");
        }
        else
        {
            gcContent = *ptrTmp;
        }
        if (bcf_get_info_float(vcfHeaderIn, line, "MAPABILITY", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/MAPABILITY");
        }
        else
        {
            mapability = *ptrTmp;
        }
        if (bcf_get_info_float(vcfHeaderIn, line, "IQR_RC0", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/IQR_RC0");
        }
        else
        {
            iqrRc0 = *ptrTmp;
        }
        if (bcf_get_info_float(vcfHeaderIn, line, "IQR_COV0", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/IQR_COV0");
        }
        else
        {
            iqrCov0 = *ptrTmp;
        }
        free(ptrTmp);
        int32_t * ptrFlagTmp = nullptr;
        sizeTmp = 0;
        bool isGap = false;
        if (bcf_get_info_flag(vcfHeaderIn, line, "GAP", &ptrFlagTmp, &sizeTmp))
            isGap = true;
        free(ptrFlagTmp);

        if (bcf_has_filter(vcfHeaderIn, line, const_cast<char *>("FEW_GCWINDOWS")) == 1 ||
                isGap ||
                gcContent < options.minGC ||
                gcContent > options.maxGC ||
                mapability < options.minMapability ||
                iqrRc0 > options.maxIqrRC ||
                iqrCov0 > options.maxIqrCov)
            continue;  // Skip noisy window

        // Extract the metric values
        if (!bcf_get_format_float(vcfHeaderIn, line, options.metric.c_str(), &valsPtr, &sizeVals))
        {
            free(valsPtr);
            throw std::runtime_error(std::string("Could not get value of FORMAT/") + options.metric);
        }

        pos.push_back(line->pos);
        for (int i = 0; i < sizeVals; ++i)
            values[i].push_back(valsPtr[i]);
    }

    if (options.verbosity >= 1)
        std::cerr << "# windows: " << pos.size() << "\n";

    int sampleID = 0;
    for (auto & measurements : values)
    {
        std::cerr << "Segmenting sample " << vcfHeaderIn->samples[sampleID++] << "\n";
        doSegmentation(pos, measurements);
    }

    if (bcf_sr_seek(vcfReaderPtr.get(), contig.c_str(), 0) != 0)
        throw std::runtime_error("Could not seek to beginning of region!");

    int idx = 0;
    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);

        // Weed out the noisy windows
        float gcContent, mapability, iqrRc0, iqrCov0;

        float * ptrTmp = nullptr;
        int32_t sizeTmp = 0;
        if (bcf_get_info_float(vcfHeaderIn, line, "GC", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/GC");
        }
        else
        {
            gcContent = *ptrTmp;
        }
        if (bcf_get_info_float(vcfHeaderIn, line, "MAPABILITY", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/MAPABILITY");
        }
        else
        {
            mapability = *ptrTmp;
        }
        if (bcf_get_info_float(vcfHeaderIn, line, "IQR_RC0", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/IQR_RC0");
        }
        else
        {
            iqrRc0 = *ptrTmp;
        }
        if (bcf_get_info_float(vcfHeaderIn, line, "IQR_COV0", &ptrTmp, &sizeTmp) != 1)
        {
            free(ptrTmp);
            throw std::runtime_error("Could not get value of INFO/IQR_COV0");
        }
        else
        {
            iqrCov0 = *ptrTmp;
        }
        free(ptrTmp);
        int32_t * ptrFlagTmp = nullptr;
        sizeTmp = 0;
        bool isGap = false;
        if (bcf_get_info_flag(vcfHeaderIn, line, "GAP", &ptrFlagTmp, &sizeTmp))
            isGap = true;
        free(ptrFlagTmp);

        if (bcf_has_filter(vcfHeaderIn, line, const_cast<char *>("FEW_GCWINDOWS")) == 1 ||
                isGap ||
                gcContent < options.minGC ||
                gcContent > options.maxGC ||
                mapability < options.minMapability ||
                iqrRc0 > options.maxIqrRC ||
                iqrCov0 > options.maxIqrCov)
            continue;  // Skip noisy window

        // Remove all metrics except for the one to segment (cannot mark for remove an update
        // at the same time).
        if (options.metric != "COV")
            bcf_update_format_float(vcfHeaderIn, line, "COV", NULL, 0);
        if (options.metric != "COV0")
            bcf_update_format_float(vcfHeaderIn, line, "COV0", NULL, 0);
        if (options.metric != "RC")
            bcf_update_format_float(vcfHeaderIn, line, "RC", NULL, 0);
        if (options.metric != "RC0")
            bcf_update_format_float(vcfHeaderIn, line, "RC0", NULL, 0);

        // Write out segment values
        std::vector<float> outVals(values.size(), 0.0);
        for (size_t i = 0; i < outVals.size(); ++i)
            outVals[i] = values[i][idx];
        if (bcf_update_format_float(
                vcfHeaderIn, line, options.metric.c_str(), &outVals[0], outVals.size()) != 0)
            throw std::runtime_error("Could not write out segmented value");

        if (bcf_translate(vcfHeaderOut, vcfHeaderIn, line))
            throw std::runtime_error("Could not translate from old to new header");
        bcf_write(vcfFileOutPtr.get(), vcfHeaderOut, line);
        ++idx;
    }

    free(valsPtr);
}

void CnvettiSegmentApp::doSegmentation(std::vector<int> const & pos, std::vector<float> & vals) const
{
    std::transform(vals.begin(), vals.end(), vals.begin(), log2);  // log2-transform values
    seg_haar(vals, options.haarBreaksFqdrQ);
    std::transform(vals.begin(), vals.end(), vals.begin(), exp2);  // transform back
}

}  // anonymous namespace


int mainSegment(CnvettiSegmentOptions const & options)
{
    CnvettiSegmentApp app(options);
    app.run();

    return 0;
}
