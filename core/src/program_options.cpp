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

#include <ostream>
#include <cstring>

namespace {  // anonymous

char const * trueFalse(bool b) {
    if (b)
        return "true";
    else
        return "false";
}

}  // anonymous namespace


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
        << "    computeBinStdevs:      " << trueFalse(computeBinStdevs) << "\n"
        << "    windowLength:          " << windowLength << "\n"
        << "    minUnclipped:          " << minUnclipped << " # percent\n"
        << "    numIOThreads:          " << numIOThreads << "\n"
        << "\n";
}


void CnvettiNormalizeOptions::print(std::ostream & out) const
{
    out
        << "options:\n"
        << "    verbosity:        " << verbosity << "\n"
        << "    inputFileName:   '" << inputFileName << "'\n"
        << "    outputFileName:  '" << outputFileName << "'\n"
        << "    numIOThreads:     " << numIOThreads << "\n"
        << "    minGCWindowCount: " << minGCWindowCount << "\n"
        << "\n";
}


void CnvettiBackgroundOptions::print(std::ostream & out) const
{
    out
        << "options:\n"
        << "    verbosity:        " << verbosity << "\n"
        << "    inputFileName:   '" << inputFileName << "'\n"
        << "    outputFileName:  '" << outputFileName << "'\n"
        << "    writeSamples:     " << trueFalse(writeSamples) << "\n"
        << "\n";
}


void CnvettiSegmentOptions::print(std::ostream & out) const
{
    out
        << "options:\n"
        << "    verbosity:        " << verbosity << "\n"
        << "    inputFileName:   '" << inputFileName << "'\n"
        << "    outputFileName:  '" << outputFileName << "'\n";

    out << "    genomeRegions";
    if (genomeRegions.empty()) {
        out << " []\n";
    } else {
        out << "\n";
        for (std::string region : genomeRegions)
            out << "    - '" << region << "'\n";
    }

    out << "    metric:            '" << metric << "'\n"
        << "    minMapability:     " << "\n"
        << "    minGC:             " << minGC << "\n"
        << "    maxGC:             " << maxGC << "\n"
        << "    maxIqrRC:          " << maxIqrRC << "\n"
        << "    maxIqrCov:         " << maxIqrCov << "\n"
        << "    haarSegBreaksFdrQ: " << haarSegBreaksFdrQ << "\n"
        << "    haarSegLmin:       " << haarSegLmin << "\n"
        << "    haarSegLmax:       " << haarSegLmax << "\n"
        << "\n";
}
