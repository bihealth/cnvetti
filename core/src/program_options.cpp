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

#include <ostream>
#include <cstring>

void CnvettiCoverageOptions::print(std::ostream & out) const
{
    out
        << "options:\n"
        << "    verbosity:             " << verbosity << "\n"
        << "    inputFileNames:";

    if (inputFileNames.empty()) {
        out << " []\n";
    } else {
        out << "\n";
        for (std::string fname : inputFileNames)
            out << "    - '" << fname << "'\n";
    }

    out << "    genomeRegions";
    if (genomeRegions.empty()) {
        out << " []\n";
    } else {
        out << "\n";
        for (std::string region : genomeRegions)
            out << "    - '" << region << "'\n";
    }

    out << "    outputFileName:        '" << outputFileName << "'\n"
        << "    mapabilityBedFileName: '" << mapabilityBedFileName << "'\n"
        << "    windowLength:          " << windowLength << "\n"
        << "    minUnclipped:          " << minUnclipped << " # percent\n"
        << "    numIOThreads:          " << numIOThreads << "\n"
        << "\n";
}


void CnvettiNormalizeOptions::print(std::ostream & out) const
{
    out
        << "options:\n"
        << "    verbosity:      " << verbosity << "\n"
        << "    inputFileName:  '" << inputFileName << "'\n"
        << "    outputFileName: '" << outputFileName << "'\n"
        << "    numIOThreads:          " << numIOThreads << "\n"
        << "\n";
}
