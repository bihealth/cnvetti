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

// Epsilon to add normalized metric to prevent -nan through log2()
float const PSEUDO_EPSILON = 1e-10;


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

    // Parse genomic regions if any
    void parseGenomicRegions();

    // Parse genomic regions
    void processGenomicRegions();
    // Perform chromosome-wise processing
    void processChromosomes();
    // Process a region
    void processRegion(std::string const & contig, int beginPos, int endPos);
    // Perform segmentation for one sample
    std::vector<float> doSegmentation(
        std::vector<int> const & pos, std::vector<float> & values) const;

    CnvettiSegmentOptions const & options;

    // Genomic regions to process
    std::vector<GenomicRegion> genomicRegions;

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
        std::cerr << "\nPARSING GENOMIC REGIONS\n\n";
    parseGenomicRegions();

    if (options.verbosity >= 1)
        std::cerr << "\nPROCESSING\n\n";
    if (genomicRegions.empty())
        processChromosomes();
    else
        processGenomicRegions();

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiSegmentApp::parseGenomicRegions()
{
    for (auto const & regionString : options.genomeRegions)
    {
        GenomicRegion region = parseGenomicRegion(regionString);

        // Look up numeric ID
        int contigID = 0;
        for (int i = 0; i < vcfHeaderIn->nhrec && region.rID == -1; ++i)
        {
            bcf_hrec_t * hrec = vcfHeaderIn->hrec[i];
            if (hrec->type == BCF_HL_CTG)
            {
                for (int j = 0; j < hrec->nkeys; ++j)
                {
                    if (strcmp("ID", hrec->keys[j]) == 0 && region.contig == hrec->vals[j])
                    {
                        region.rID = contigID;
                        break;
                    }
                }
                ++contigID;
            }
        }
        if (region.rID == -1)
            throw std::runtime_error(
                std::string("Contig/chromosome ") + region.contig + " not see in BAM file");

        // Add to genomic regions
        genomicRegions.push_back(region);
    }

    std::sort(genomicRegions.begin(), genomicRegions.end(),
            [](GenomicRegion const & lhs, GenomicRegion const & rhs) {
                return (std::make_pair(lhs.rID, lhs.beginPos) <
                    std::make_pair(rhs.rID, rhs.beginPos));
            });
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
        ("##FORMAT=<ID=SEG,Type=Float,Number=1,Description="
         "\"Segment value\">"));
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

void CnvettiSegmentApp::processGenomicRegions()
{
    for (auto const region : genomicRegions)
        processRegion(region.contig, region.beginPos, region.endPos);

    if (options.verbosity >= 1)
        std::cerr << " DONE\n";
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

    if (bcf_sr_seek(vcfReaderPtr.get(), contig.c_str(), beginPos) != 0)
    {
        std::cerr << "WARNING: could not seek to beginning of " << contig << ":" << beginPos << "-" << endPos << "\n";
        return;
    }

    std::vector<int> pos;
    std::vector<std::vector<float>> values;
    values.resize(bcf_hdr_nsamples(vcfHeaderIn));

    float * valsPtr = nullptr;
    int32_t sizeVals = 0;

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (!(beginPos == 0 && endPos == 0) && line->pos >= endPos)
            break;

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
    std::vector<std::vector<float>> segmented(values.size());
    for (auto & measurements : values)
    {
        std::cerr << "Segmenting sample " << vcfHeaderIn->samples[sampleID] << "\n";
        segmented[sampleID] = doSegmentation(pos, measurements);
        ++sampleID;
    }

    if (bcf_sr_seek(vcfReaderPtr.get(), contig.c_str(), beginPos) != 0)
        throw std::runtime_error("Could not seek to beginning of region!");

    int idx = 0;
    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (!(beginPos == 0 && endPos == 0) && line->pos >= endPos)
            break;

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

        // Translate record from input to output file header
        if (bcf_translate(vcfHeaderOut, vcfHeaderIn, line))
            throw std::runtime_error("Could not translate from old to new header");

        // Put segmentation values into record to write
        std::vector<float> outVals(segmented.size(), 0.0);
        for (size_t i = 0; i < outVals.size(); ++i)
            outVals[i] = segmented[i][idx];
        if (bcf_update_format_float(
                vcfHeaderOut, line, "SEG", &outVals[0], outVals.size()) != 0)
            throw std::runtime_error("Could not write out segmented value");

        // Write out the record
        bcf_write(vcfFileOutPtr.get(), vcfHeaderOut, line);
        ++idx;
    }

    free(valsPtr);
}

std::vector<float> CnvettiSegmentApp::doSegmentation(
    std::vector<int> const & pos, std::vector<float> & vals) const
{
    std::transform(vals.begin(), vals.end(), vals.begin(),
                   [](float x) { return x + PSEUDO_EPSILON; });
    std::transform(vals.begin(), vals.end(), vals.begin(), log2);  // log2-transform values

    std::vector<size_t> breakpoints = segmentHaarSeg(
        vals, options.haarSegBreaksFdrQ, nullptr, nullptr, options.haarSegLmin,
        options.haarSegLmax);
    std::vector<float> segmented = replaceWithSegmentMedians(vals, breakpoints);

    std::transform(segmented.begin(), segmented.end(), segmented.begin(), exp2);  // transform back
    std::transform(segmented.begin(), segmented.end(), segmented.begin(),
                   [](float x) { return x - PSEUDO_EPSILON; });
    return segmented;
}

}  // anonymous namespace


int mainSegment(CnvettiSegmentOptions const & options)
{
    CnvettiSegmentApp app(options);
    app.run();

    return 0;
}
