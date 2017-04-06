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

#include <tclap/CmdLine.h>


// Forward Declarations -------------------------------------------------------

class CnvettiSummariesOptions;

int mainCoverage(CnvettiCoverageOptions const & options);
int mainSummaries(CnvettiSummariesOptions const & options);


// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    try {
        CnvettiCommand cmd = parseTopLevelCommandLine(argc, argv);
        switch (cmd) {
        case CnvettiCommand::NONE:
            std::cout
                << "Usage: cnvetti <command> [option]\n"
                << "Try `cnvetti --help` for more information\n";
            return 1;

        case CnvettiCommand::VERSION:
            std::cout << "cnvetti " << GIT_VERSION << "\n";
            return 0;

        case CnvettiCommand::TOPLEVEL_HELP:
            printTopLevelHelp(std::cerr);
            return 0;

        case CnvettiCommand::COVERAGE:
            return mainCoverage(parseCoverageCommandLine(argc, argv));
        }
    } catch (TCLAP::ArgException & e) {
        std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << "\n\n";
        printTopLevelHelp(std::cerr);
        return 1;
    }

    return 0;
}
