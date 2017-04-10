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

#include "cnvetti/program_options.h"

#include <iostream>
#include <memory>
#include <vector>

#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include "cnvetti/program_options.h"
#include "cnvetti/version.h"
#include "cnvetti/utils.h"

namespace {  // anonymous namespace

// ---------------------------------------------------------------------------
// Class CnvettiNormalizeApp
// ---------------------------------------------------------------------------

class CnvettiNormalizeApp
{
public:
    CnvettiNormalizeApp(CnvettiNormalizeOptions const & options) :
            options(options)
    {}

    void run();

private:

    // Open the input files and store them in bamFilesIn and open vcfFileOut;
    void openFiles();
    // Check that the references FASTA file fits to the BAM file.
    void checkFiles();
    // Perform chromosome-wise processing
    void processChromosomes();

    CnvettiNormalizeOptions const & options;

    std::shared_ptr<bcf_srs_t> vcfReader;
    std::shared_ptr<vcfFile> vcfFileOut;
    bcf_hdr_t * vcfHeaderIn;
    std::shared_ptr<bcf_hdr_t> vcfHeaderOut;
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
        std::cerr << "\nPROCESSING\n\n";
    processChromosomes();

    if (options.verbosity >= 1)
        std::cerr << "\nDone. Have a nice day!\n";
}

void CnvettiNormalizeApp::openFiles()
{
    // Open input VCF file and read header
    if (options.verbosity >= 1)
        std::cerr << "Opening " << options.inputFileName << " ...\n";
    vcfReader = std::shared_ptr<bcf_srs_t>(bcf_sr_init(), bcf_sr_destroy);
    bcf_sr_set_opt(vcfReader.get(), BCF_SR_REQUIRE_IDX);
    if (bcf_sr_set_threads(vcfReader.get(), options.numIOThreads) < 0)
        throw std::runtime_error("Failed to create threads");
    if (!bcf_sr_add_reader(vcfReader.get(), options.inputFileName.c_str()))
        throw std::runtime_error(std::string("Failed to open ") + options.inputFileName);
    vcfHeaderIn = vcfReader.get()->readers[0].header;

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
    vcfFileOut = std::shared_ptr<vcfFile>(
        hts_open(options.outputFileName.c_str(), openMode.c_str()),
        hts_close);
    hts_set_threads(vcfFileOut.get(), options.numIOThreads);
    vcfFileOut->format.format = openFormat;
    vcfFileOut->format.compression = openCompression;
    // if (vcfFileOut->format.compression == bgzf)
    //     bgzf_index_build_init(vcfFileOut->fp.bgzf);

    // Copy input to output header
    vcfHeaderOut = std::shared_ptr<bcf_hdr_t>(bcf_hdr_dup(vcfHeaderIn), bcf_hdr_destroy);
    // Remove pseudo-contig with statistics and FILTER
    bcf_hdr_remove(vcfHeaderOut.get(), BCF_HL_CTG, "__cnvetti_stats");
    bcf_hdr_remove(vcfHeaderOut.get(), BCF_HL_FLT, "CNVETTI_STATS");
    bcf_hdr_remove(vcfHeaderOut.get(), BCF_HL_INFO, "GCWINDOWS");

    std::string cmdLine = options.argv[1];
    for (int i = 2; i < options.argc; ++i) {
        cmdLine += " ";
        cmdLine += options.argv[i];
    }
    bcf_hdr_printf(vcfHeaderOut.get(), "##cnvetti_normalizeCommand=%s", cmdLine.c_str());

    bcf_hdr_write(vcfFileOut.get(), vcfHeaderOut.get());
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

void CnvettiNormalizeApp::processChromosomes()
{

}

}  // anonymous namespace

int mainNormalize(CnvettiNormalizeOptions const & options)
{
    CnvettiNormalizeApp app(options);
    app.run();

    return 0;
}
