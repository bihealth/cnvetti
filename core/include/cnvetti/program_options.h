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
// Class InvalidArgumentsException
// ----------------------------------------------------------------------------

class InvalidCommandLineArgumentsException : public std::runtime_error
{
public:
    InvalidCommandLineArgumentsException() : std::runtime_error("")
    {}
};

// ----------------------------------------------------------------------------
// Enum CnvettiCommand
// ----------------------------------------------------------------------------

// Enum for describing the available commands

enum class CnvettiCommand
{
    // Print version and exit(0)
    VERSION,
    // Print top-level use help and exit(0)
    TOPLEVEL_HELP,
    // No command given, print help and exit
    NONE,
    // Run "cnvetti coverage"
    COVERAGE,
    // Run "cnvetti summaries"
    SUMMARIES
};

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

    // Length of the windows in bp.
    int windowLength;

    // Minimal unclipped fraction of reads to consider, in percent.
    int minUnclipped;

    CnvettiCoverageOptions() : verbosity(1), windowLength(1000), minUnclipped(80)
    {}

    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Class CnvettiSummariesOptions
// ----------------------------------------------------------------------------

// Program options for summarising merged results of `convetti coverage`

class CnvettiSummariesOptions
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

    // Length of the windows in bp.
    int windowLength;

    // Minimal unclipped fraction of reads to consider, in percent.
    int minUnclipped;

    CnvettiSummariesOptions() : verbosity(1), windowLength(1000), minUnclipped(80)
    {}

    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Function printTopLevelHelp()
// ----------------------------------------------------------------------------

void printTopLevelHelp(std::ostream & out);

// ----------------------------------------------------------------------------
// Function parseTopLevelCommandLine()
// ----------------------------------------------------------------------------

CnvettiCommand parseTopLevelCommandLine(int argc, char ** argv);

// ----------------------------------------------------------------------------
// Function parseCoverageCommandLine()
// ----------------------------------------------------------------------------

CnvettiCoverageOptions parseCoverageCommandLine(int argc, char ** argv);
