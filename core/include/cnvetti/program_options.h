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

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiCoverageOptions() :
        verbosity(1), numIOThreads(1), windowLength(1000), minUnclipped(80),
        argc(0), argv(nullptr)
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

    // Skip any GC content with fewer than this number of windows
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

    // FQDR q parameter
    double haarBreaksFqdrQ;

    // argc and argv from command line
    int argc;
    char ** argv;

    CnvettiSegmentOptions() :
        verbosity(1), argc(0), argv(nullptr), numIOThreads(1), metric("COV0"),
        minMapability(0.99), minGC(20), maxGC(70), maxIqrRC(0.25), maxIqrCov(0.2),
        haarBreaksFqdrQ(1e-3)
    {}

    void print(std::ostream & out) const;
};
