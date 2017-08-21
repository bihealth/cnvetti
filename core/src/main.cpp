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
#include "cnvetti/version.h"

#include <iostream>

#include <CLI/CLI.hpp>
#include <CLI/Timer.hpp>


// ----------------------------------------------------------------------------
// Forward Declarations
// ----------------------------------------------------------------------------

int mainCoverage(CnvettiCoverageOptions const & options);
int mainJoin(CnvettiJoinOptions const & options);
int mainPeaks(CnvettiPeaksOptions const & options);
int mainNormalize(CnvettiNormalizeOptions const & options);
int mainRatio(CnvettiRatioOptions const & options);
int mainBackground(CnvettiBackgroundOptions const & options);
int mainSegment(CnvettiSegmentOptions const & options);

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
        "Path to mapability BED file"
    )->check(CLI::ExistingFile)->group("Input / Output");
    covApp->add_option(
        "--target-bed", covOptions.targetBedFile,
        "Path to target BED file (optional)"
    )->check(CLI::ExistingFile)->group("Input / Output");
    covApp->add_option(
        "--num-io-threads", covOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    covApp->add_flag(
        "--ignore-off-target", covOptions.ignoreOffTarget,
        "Ignore off-target bins (--target-bed must be given)"
    )->group("Algorithm Parameters");
    covApp->add_flag(
        "--ignore-discordant-pairs", covOptions.ignoreDiscordantPairs,
        "Ignore discordant pairs (for second pass, targeted seq.)"
    )->group("Algorithm Parameters");
    covApp->add_flag(
        "--compute-bin-stdevs", covOptions.computeBinStdevs,
        "Whether or not to compute per-bin stdevs in coverage"
    )->group("Algorithm Parameters");
    covApp->add_option(
        "--window-length", covOptions.windowLength,
        "Window length to use"
    )->group("Algorithm Parameters");
    covApp->add_option(
        "--min-alignment-quality", covOptions.minAlignmentQuality,
        "Minimal alignment quality to enforce"
    )->group("Algorithm Parameters");
    covApp->add_option(
        "--gc-step-size", covOptions.gcStepSize,
        "GC step size to use"
    )->group("Algorithm Parameters");
    covApp->add_option(
        "--min-unclipped", covOptions.minUnclipped,
        "Minimal unclipped fraction of reads to keep, in percent"
    )->group("Algorithm Parameters");

    covApp->add_option(
        "--peaks-bed-file", covOptions.peaksBedFile,
        "Peaks BED file from 'cnvetti peaks' output"
    )->group("Off-Target Read Processing");
    covApp->add_option(
        "--peak-boundary", covOptions.peakBoundary,
        "Number of base pairs to ignore adjacent to peak windows"
    )->group("Off-Target Read Processing");

    // Add sub command `cnvetti peaks`

    CnvettiPeaksOptions peaksOptions;
    peaksOptions.argc = argc;
    peaksOptions.argv = argv;

    CLI::App * peaksApp = app.add_subcommand(
        "peaks", "Compute peaks (amplified) regions BED file");

    peaksApp->add_option(
        "-i,--input", peaksOptions.inputFileName,
        "Path to input BCF file (required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    peaksApp->add_option(
        "-o,--output", peaksOptions.outputFileName,
        "Path to output BED file (required)"
    )->required()->group("Input / Output");
    peaksApp->add_option(
        "--num-io-threads", peaksOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    peaksApp->add_option(
        "--percentile", peaksOptions.percentile,
        "Percentile to select for threshold, in percent"
    )->group("Algorithm Parameters");
    peaksApp->add_option(
        "--thresh-factor", peaksOptions.threshFactor,
        "Factor to multiply threshold with after selecting by percentile"
    )->group("Algorithm Parameters");

    // Add sub command `cnvetti join`

    CnvettiJoinOptions joinOptions;
    joinOptions.argc = argc;
    joinOptions.argv = argv;

    CLI::App * joinApp = app.add_subcommand(
        "join", "Join adjacent windows (for on-target coverage)");

    joinApp->add_option(
        "-i,--input", joinOptions.inputFileName,
        "Path to input BCF file (required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    joinApp->add_option(
        "-o,--output", joinOptions.outputFileName,
        "Path to output BED file (required)"
    )->required()->group("Input / Output");
    joinApp->add_option(
        "--num-io-threads", joinOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    // Add sub command `cnvetti normalize`

    CnvettiNormalizeOptions normOptions;
    normOptions.argc = argc;
    normOptions.argv = argv;

    CLI::App * cnvettiNormalize = app.add_subcommand(
        "normalize", "Perform per-sample normalization on cnvetti coverage output");
    cnvettiNormalize->add_option(
        "-i,--input", normOptions.inputFileName,
        "Path to input VCF/BCF file (from cnvetti coverage; required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    cnvettiNormalize->add_option(
        "-o,--output", normOptions.outputFileName,
        "Path to output VCF/BCF file (required)"
    )->required()->group("Input / Output");
    cnvettiNormalize->add_option(
        "--num-io-threads", normOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    cnvettiNormalize->add_option(
        "--min-gc-window-count", normOptions.minGCWindowCount,
        "Flag GC contents with fewer windows as FEW_GCWINDOWS"
    )->group("Algorithm Parameters");

    // Add sub command `cnvetti ratio`

    CnvettiRatioOptions ratioOptions;
    ratioOptions.argc = argc;
    ratioOptions.argv = argv;

    CLI::App * cnvettiRatio = app.add_subcommand(
        "ratio", ("Compute ratios from cnvetti normalize output (for cancer matched cancer "
                  "samples)"));
    cnvettiRatio->add_option(
        "-i,--input", ratioOptions.inputFileName,
        "Path to input VCF/BCF file (from cnvetti normalize; required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    cnvettiRatio->add_option(
        "-o,--output", ratioOptions.outputFileName,
        "Path to output VCF/BCF file (required)"
    )->required()->group("Input / Output");
    cnvettiRatio->add_option(
        "-n,--normal-sample", ratioOptions.normalSample,
        "Name of the matched normal sample in BCF file"
    )->required()->group("Input / Output");
    cnvettiRatio->add_option(
        "-t,--tumor-sample", ratioOptions.tumorSample,
        "Name of the matched tumor sample in BCF file"
    )->required()->group("Input / Output");
    cnvettiRatio->add_option(
        "--num-io-threads", ratioOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    // Add sub command `cnvetti background`

    CnvettiBackgroundOptions bgOptions;
    bgOptions.argc = argc;
    bgOptions.argv = argv;

    CLI::App * cnvettiBackground = app.add_subcommand(
        "background", "Build background from cross-cohort normalized data");
    cnvettiBackground->add_option(
        "-i,--input", bgOptions.inputFileName,
        "Path to input VCF/BCF file (from cnvetti normalize; required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    cnvettiBackground->add_option(
        "-o,--output", bgOptions.outputFileName,
        "Path to output VCF/BCF file (required)"
    )->required()->group("Input / Output");
    cnvettiBackground->add_option(
        "--num-io-threads", bgOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");
    cnvettiBackground->add_flag(
        "--write-samples", bgOptions.writeSamples,
        "Write out samples, if omitted only INFO fields will be written"
    )->group("Input / Output");

    // Add sub command `cnvetti segment`

    CnvettiSegmentOptions segOptions;
    segOptions.argc = argc;
    segOptions.argv = argv;

    CLI::App * cnvettiSegment = app.add_subcommand(
        "segment", "Perform segmentation from normalized coverage and background");
    cnvettiSegment->add_option(
        "-i,--input", segOptions.inputFileName,
        "Path to input VCF/BCF file (cnvettig background and 1+ sample(s); required)"
    )->required()->check(CLI::ExistingFile)->group("Input / Output");
    cnvettiSegment->add_option(
        "-o,--output", segOptions.outputFileName,
        "Path to output VCF/BCF file (required)"
    )->required()->group("Input / Output");
    cnvettiSegment->add_option(
        "--genome-regions", segOptions.genomeRegions,
        "Genome regions to process"
    )->group("Input / Output");
    cnvettiSegment->add_option(
        "--num-io-threads", segOptions.numIOThreads,
        "Number of threads to use for de-/compression in I/O"
    )->group("Input / Output");

    cnvettiSegment->add_option(
        "--metric", segOptions.metric,
        "Metric to use for segmentation (allowed: COV, COV0, RC, RC0; default: COV0)"
    )->group("Segmentation Algorithm")->check([](std::string const & val) {
            std::vector<std::string> const ALLOWED_METRICS { "COV", "COV0", "RC", "RC0" };
            return std::find(ALLOWED_METRICS.begin(), ALLOWED_METRICS.end(), val) != ALLOWED_METRICS.end();
    });
    cnvettiSegment->add_option(
        "--min-mapability", segOptions.minMapability,
        "Minimal mapability value of windows to consider (required)"
    )->required()->group("Segmentation Algorithm");
    cnvettiSegment->add_option(
        "--min-gc", segOptions.minGC,
        "Minimal GC value of windows to consider"
    )->group("Segmentation Algorithm");
    cnvettiSegment->add_option(
        "--max-gc", segOptions.maxGC,
        "Maxmial GC value of windows to consider"
    )->group("Segmentation Algorithm");
    cnvettiSegment->add_option(
        "--max-iqr-rc", segOptions.maxIqrRC,
        "Minimal read count IQR of windows to consider"
    )->group("Segmentation Algorithm");
    cnvettiSegment->add_option(
        "--max-iqr-cov", segOptions.maxIqrCov,
        "Maximal coverage IQR of windows to consider"
    )->group("Segmentation Algorithm");
    cnvettiSegment->add_option(
        "--haar-seg-breaks-fdr-q", segOptions.haarSegBreaksFdrQ,
        "FDR q value for Haar wavelet--based segmentation"
    )->group("Segmentation Algorithm");

    try {
        app.parse(argc, argv);

        if (print_version) {
            std::cout << "cnvetti " << GIT_VERSION << "\n";
            return 0;
        }

        if (app.got_subcommand("coverage")) {
            CLI::AutoTimer timer("running time");
            mainCoverage(covOptions);
        } else if (app.got_subcommand("peaks")) {
            CLI::AutoTimer timer("running time");
            mainPeaks(peaksOptions);
        } else if (app.got_subcommand("join")) {
            CLI::AutoTimer timer("running time");
            mainJoin(joinOptions);
        } else if (app.got_subcommand("normalize")) {
            CLI::AutoTimer timer("running time");
            mainNormalize(normOptions);
        } else if (app.got_subcommand("background")) {
            CLI::AutoTimer timer("running time");
            mainBackground(bgOptions);
        } else if (app.got_subcommand("ratio")) {
            CLI::AutoTimer timer("running time");
            mainRatio(ratioOptions);
        } else if (app.got_subcommand("segment")) {
            CLI::AutoTimer timer("running time");
            mainSegment(segOptions);
        } else {
            throw CLI::CallForHelp();
        }
    } catch (const CLI::ParseError &e) {
        // Will make the app exit with 0 on --help
        return app.exit(e);
    }

    return 0;
}
