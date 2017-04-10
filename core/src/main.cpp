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
#include "cnvetti/version.h"

#include <iostream>

#include <CLI/CLI.hpp>
#include <CLI/Timer.hpp>


// ----------------------------------------------------------------------------
// Forward Declarations
// ----------------------------------------------------------------------------

int mainCoverage(CnvettiCoverageOptions const & options);
int mainNormalize(CnvettiNormalizeOptions const & options);


// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    CLI::App app{"CNV calling from NGS data"};

    // Global, top-level options
    bool print_version = false;
    int verbosity = 0;

    app.add_flag("--version", print_version, "Print version and exit");
    app.add_flag("-v,--verbose", verbosity, "Verbosity, pass N times for level N");

    // Add sub command `cnvetti coverage`

    CnvettiCoverageOptions covOptions;
    covOptions.argc = argc;
    covOptions.argv = argv;

    CLI::App * covApp = app.add_subcommand(
        "coverage", "Collect per-sample coverage information");

    covApp->add_option(
        "-r,--reference", covOptions.referencePath,
        "Path to FAI-indexed FASTA file (required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    covApp->add_option(
        "-i,--input", covOptions.inputFileNames,
        "Path to input BAM file (required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    covApp->add_option(
        "-o,--output", covOptions.outputFileName,
        "Path to output VCF/BCF file (required)"
    )->required()->group("Input / Output");
    covApp->add_option(
        "--genome-regions", covOptions.genomeRegions,
        "Genome regions to process"
    )->group("Input / Output");
    covApp->add_option(
        "--mapability-bed", covOptions.mapabilityBedFileName,
        "Path to mapability BED file (required)"
    )->check(CLI::ExistingFile)->group("Input / Output");
    covApp->add_option(
        "--num-io-threads", covOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    covApp->add_option(
        "--window-length", covOptions.windowLength,
    "Window length to use"
    )->group("Algorithm Parameters");
    covApp->add_option(
        "--min-unclipped", covOptions.minUnclipped,
        "Minimal unclipped fraction of reads to keep, in percent"
    )->group("Algorithm Parameters");

    // Add sub command `cnvetti normalize`

    CnvettiNormalizeOptions normOptions;
    normOptions.argc = argc;
    normOptions.argv = argv;

    CLI::App * cnvettiNormalize = app.add_subcommand(
        "normalize", "Perform per-sample normalization on cnvetti coverage output");
    cnvettiNormalize->add_option(
        "-i,--input", normOptions.inputFileName,
        "Path to input BAM file (required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    cnvettiNormalize->add_option(
        "-o,--output", normOptions.outputFileName,
        "Path to output VCF/BCF file (required)"
    )->required()->group("Input / Output");
    cnvettiNormalize->add_option(
        "--num-io-threads", normOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    try {
        app.parse(argc, argv);

        if (print_version) {
            std::cout << "cnvetti " << GIT_VERSION << "\n";
            return 0;
        }

        if (app.got_subcommand("coverage")) {
            CLI::AutoTimer timer("running time");
            mainCoverage(covOptions);
        } else if (app.got_subcommand("normalize")) {
            CLI::AutoTimer timer("running time");
            mainNormalize(normOptions);
        } else {
            throw CLI::CallForHelp();
        }
    } catch (const CLI::ParseError &e) {
        // Will make the app exit with 0 on --help
        return app.exit(e);
    }

    return 0;
}
