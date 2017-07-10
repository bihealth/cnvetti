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

#include <cstdlib>
#include <iostream>
#include <memory>
#include <set>
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
// Class CnvettiRatioApp
// ---------------------------------------------------------------------------

class CnvettiRatioApp
{
public:
    CnvettiRatioApp(CnvettiRatioOptions const & options) :
            options(options), vcfHeaderIn(nullptr), vcfHeaderOut(nullptr)
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

    CnvettiRatioOptions const & options;

    std::shared_ptr<bcf_srs_t> vcfReaderPtr;
    std::shared_ptr<vcfFile> vcfFileOutPtr;
    std::shared_ptr<bcf_hdr_t> vcfHeaderOutPtr;
    bcf_hdr_t * vcfHeaderIn;
    bcf_hdr_t * vcfHeaderOut;

    int nSamplesOut;
    std::vector<int> idMap;
};


void CnvettiRatioApp::run()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "cnvetti ratio\n"
                  << "=============\n";
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

void CnvettiRatioApp::openFiles()
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

    // Copy input to output header, removing tumor sample.
    nSamplesOut = bcf_hdr_nsamples(vcfHeaderIn) - 1;
    std::vector<char *> samplesOut;
    for (int i = 0; i < bcf_hdr_nsamples(vcfHeaderIn); ++i) {
        if (options.normalSample != vcfHeaderIn->samples[i]) {
            samplesOut.push_back(vcfHeaderIn->samples[i]);
        }
    }
    idMap.resize(bcf_hdr_nsamples(vcfHeaderIn), -1);
    vcfHeaderOutPtr = std::shared_ptr<bcf_hdr_t>(
            bcf_hdr_subset(vcfHeaderIn, nSamplesOut, &samplesOut[0], &idMap[0]),
            bcf_hdr_destroy);
    vcfHeaderOut = vcfHeaderOutPtr.get();

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeaderOut, "##cnvetti_ratioCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOutPtr.get(), vcfHeaderOut);
}

void CnvettiRatioApp::checkFiles()
{
    // Check that the file contains the two samples.
    std::set<std::string> samples;
    for (int i = 0; i < bcf_hdr_nsamples(vcfHeaderIn); ++i) {
        samples.insert(vcfHeaderIn->samples[i]);
    }
    if (samples.count(options.tumorSample) == 0u) {
        throw std::runtime_error(
                "Could not find tumor sample " + options.tumorSample + " in input file.");
    }
    if (samples.count(options.normalSample) == 0u) {
        throw std::runtime_error(
                "Could not find normal sample " + options.normalSample + " in input file.");
    }
    if (options.normalSample == options.tumorSample) {
        throw std::runtime_error("Normal sample name == tumor sample name!");
    }

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

void CnvettiRatioApp::processChromosomes()
{
    for (int i = 0; i < vcfHeaderIn->nhrec; ++i)
    {
        bcf_hrec_t * hrec = vcfHeaderIn->hrec[i];
        if (hrec->type == BCF_HL_CTG)
        {
            for (int j = 0; j < hrec->nkeys; ++j)
            {
                if (strcmp("ID", hrec->keys[j]) == 0)
                {
                    processRegion(hrec->vals[j], 0, 0);
                }
            }
        }
    }
}

void CnvettiRatioApp::processRegion(std::string const & contig, int beginPos, int endPos)
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

    int tumorSampleId = -1;
    int normalSampleId = -1;
    for (int i = 0; i < bcf_hdr_nsamples(vcfHeaderIn); ++i) {
        if (options.tumorSample == vcfHeaderIn->samples[i]) {
            tumorSampleId = i;
        } else if (options.normalSample == vcfHeaderIn->samples[i]) {
            normalSampleId = i;
        }
    }

    while (bcf_sr_next_line(vcfReaderPtr.get()))
    {
        bcf1_t * line = bcf_sr_get_line(vcfReaderPtr.get(), 0);
        if (contig != vcfReaderPtr.get()->regions->seq_names[line->rid])
            break;  // on different contig

        // Compute COV ratios
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderIn, line, "COV", &values, &sizeVals);
            if (!res) {
                free(values);
                throw std::runtime_error("Could not get FORMAT/COV value");
            }
            if (values[normalSampleId] == 0.0) {
                bcf_float_set_missing(values[tumorSampleId]);
            } else {
                values[tumorSampleId] /= values[normalSampleId];
            }
            bcf_update_format_float(vcfHeaderIn, line, "COV", values, sizeVals);
            free(values);
        }

        // Compute COV0 ratios
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderIn, line, "COV0", &values, &sizeVals);
            if (!res) {
                free(values);
                throw std::runtime_error("Could not get FORMAT/COV0 value");
            }
            if (values[normalSampleId] == 0.0) {
                bcf_float_set_missing(values[tumorSampleId]);
            } else {
                values[tumorSampleId] /= values[normalSampleId];
            }
            bcf_update_format_float(vcfHeaderIn, line, "COV0", values, sizeVals);
            free(values);
        }

        // Compute WINSD ratios
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderIn, line, "WINSD", &values, &sizeVals);
            if (!res) {
                free(values);
                throw std::runtime_error("Could not get FORMAT/WINSD value");
            }
            if (values[normalSampleId] == 0.0) {
                bcf_float_set_missing(values[tumorSampleId]);
            } else {
                values[tumorSampleId] /= values[normalSampleId];
            }
            bcf_update_format_float(vcfHeaderIn, line, "WINSD", values, sizeVals);
            free(values);
        }

        // Compute WINSD0 ratios
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderIn, line, "WINSD0", &values, &sizeVals);
            if (!res) {
                free(values);
                throw std::runtime_error("Could not get FORMAT/WINSD0 value");
            }
            if (values[normalSampleId] == 0.0) {
                bcf_float_set_missing(values[tumorSampleId]);
            } else {
                values[tumorSampleId] /= values[normalSampleId];
            }
            bcf_update_format_float(vcfHeaderIn, line, "WINSD0", values, sizeVals);
            free(values);
        }

        // Compute RC ratios
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderIn, line, "RC", &values, &sizeVals);
            if (!res) {
                free(values);
                throw std::runtime_error("Could not get FORMAT/RC value");
            }
            if (values[normalSampleId] == 0.0) {
                bcf_float_set_missing(values[tumorSampleId]);
            } else {
                values[tumorSampleId] /= values[normalSampleId];
            }
            bcf_update_format_float(vcfHeaderIn, line, "RC", values, sizeVals);
            free(values);
        }

        // Compute RC0 ratios
        {
            float * values = nullptr;
            int sizeVals = 0;
            int res = bcf_get_format_float(vcfHeaderIn, line, "RC0", &values, &sizeVals);
            if (!res) {
                free(values);
                throw std::runtime_error("Could not get FORMAT/RC0 value");
            }
            if (values[normalSampleId] == 0.0) {
                bcf_float_set_missing(values[tumorSampleId]);
            } else {
                values[tumorSampleId] /= values[normalSampleId];
            }
            bcf_update_format_float(vcfHeaderIn, line, "RC0", values, sizeVals);
            free(values);
        }

        // Subset record before writing out.
        if (bcf_subset(vcfHeaderOut, line, nSamplesOut, &idMap[0]) != 0) {
            throw std::runtime_error("Could not subset record.");
        }

        bcf_write(vcfFileOutPtr.get(), vcfHeaderOut, line);
    }
}

}  // anonymous namespace

int mainRatio(CnvettiRatioOptions const & options)
{
    CnvettiRatioApp app(options);
    app.run();

    return 0;
}
