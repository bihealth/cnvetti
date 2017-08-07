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

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
// Class CnvettiCoverageOptions
// ----------------------------------------------------------------------------

// Program options for computing per-sample coverage BCF File

class CnvettiCoverageOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to reference.
    std::string referencePath;

    // Paths to input file names, at least one.
    std::vector<std::string> inputFileNames;

    // Genomic regions to process
    std::vector<std::string> genomeRegions;

    // Path to output file.
    std::string outputFileName;

    // Path to tabix-indexed, bgzip-ed mapability file.
    std::string mapabilityBedFileName;

    // Number of threads to use for decompression/compression on I/O
    int numIOThreads;

    // Length of the windows in bp.
    int windowLength;

    // Minimal unclipped fraction of reads to consider, in percent.
    int minUnclipped;

    // Whether or not to compute bin stdev values
    bool computeBinStdevs;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiCoverageOptions() :
        verbosity(1), numIOThreads(1), windowLength(500), minUnclipped(80),
        computeBinStdevs(false), argc(0), argv(nullptr)
    {}

    void print(std::ostream & out) const;
};


// ----------------------------------------------------------------------------
// Class CnvettiPeaksOptions
// ----------------------------------------------------------------------------

class CnvettiPeaksOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to input file name
    std::string inputFileName;

    // Path to output file.
    std::string outputFileName;

    // Percentile for selecting threshold by;
    double percentile;

    // Factor to multiply threshold with;
    double threshFactor;

    // Number of threads to use for decompression/compression on I/O
    int numIOThreads;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiPeaksOptions() :
        verbosity(1), argc(0), argv(nullptr), percentile(90.0),
        threshFactor(2.0), numIOThreads(1)
    {}

    void print(std::ostream & out) const;
};


// ----------------------------------------------------------------------------
// Class CnvettiNormalizeOptions
// ----------------------------------------------------------------------------

// Program options for the normalization of `cnvetti coverage` results

class CnvettiNormalizeOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to input file name
    std::string inputFileName;

    // Path to output file.
    std::string outputFileName;

    // Flag windows with GC content for which there are fewer than this number of windows
    int minGCWindowCount;

    // Number of threads to use for decompression/compression on I/O
    int numIOThreads;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiNormalizeOptions() :
        verbosity(1), argc(0), argv(nullptr), minGCWindowCount(100), numIOThreads(1)
    {}

    void print(std::ostream & out) const;
};


// ----------------------------------------------------------------------------
// Class CnvettiRatioOptions
// ----------------------------------------------------------------------------

// Program options for the computation of ratios from the `cnvetti normalize`
// results.

class CnvettiRatioOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to input file name
    std::string inputFileName;

    // Path to output file.
    std::string outputFileName;

    // Name of the matched tumor sample.
    std::string tumorSample;

    // Name of the matched normal sample.
    std::string normalSample;

    // Number of threads to use for decompression/compression on I/O.
    int numIOThreads;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiRatioOptions() :
        verbosity(1), numIOThreads(1), argc(0), argv(nullptr)
    {}

    void print(std::ostream & out) const;
};


// ----------------------------------------------------------------------------
// Class CnvettiBackgroundOptions
// ----------------------------------------------------------------------------

// Program options for computing background from multiple `cnvetti normalize`
// results combined with `bcftools merge`

class CnvettiBackgroundOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to input file name
    std::string inputFileName;

    // Path to output file.
    std::string outputFileName;

    // Write out samples with additional INFO values
    bool writeSamples;

    // Number of threads to use for decompression/compression on I/O
    int numIOThreads;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiBackgroundOptions() :
        verbosity(1), argc(0), argv(nullptr), writeSamples(false), numIOThreads(1)
    {}

    void print(std::ostream & out) const;
};


// ----------------------------------------------------------------------------
// Class CnvettiSegmentOptions
// ----------------------------------------------------------------------------

// Program options for segmentation based on the normalized coverage of
// individuals with a given background

class CnvettiSegmentOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

    // Path to input file
    std::string inputFileName;

    // Path to output file
    std::string outputFileName;

    // Genomic regions to process
    std::vector<std::string> genomeRegions;

    // Number of threads to use for decompression/compression on I/O
    int numIOThreads;

    // Metric to base the segmentation on
    std::string metric;

    // Minimal mapability in window
    double minMapability;

    // Minimal GC content in window
    double minGC;

    // Maximal GC content in window
    double maxGC;

    // Maximal normalized IQR of normalized read count
    double maxIqrRC;

    // Maximal normalized IQR of normalized coverage
    double maxIqrCov;

    // Parameters for HaarSeg algorithm -----------------------------------------------------------

    // The L_MIN parameter for HaarSeg, should be ceil(log2(k))
    int haarSegLmin;

    // The L_MAX parameter for HaarSeg, should be large enough
    int haarSegLmax;

    // FDR parameter Q
    double haarSegBreaksFdrQ;

    // Whether to use the background or not.
    bool useBackground;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiSegmentOptions() :
        verbosity(1), argc(0), argv(nullptr), numIOThreads(1), metric("COV0"),
        minMapability(0.99), minGC(20), maxGC(70), maxIqrRC(0.25), maxIqrCov(0.2),
        haarSegBreaksFdrQ(1e-7), haarSegLmin(1), haarSegLmax(5), useBackground(false)
    {}

    void print(std::ostream & out) const;
};
