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

#include <cstring>

#include <tclap/CmdLine.h>

#include "cnvetti/version.h"

// ----------------------------------------------------------------------------
// Function printTopLevelHelp()
// ----------------------------------------------------------------------------

void printTopLevelHelp(std::ostream & out)
{
    out
        << "Program: CNVetti (CNV calling from WGS data)\n"
        << "Version: " << GIT_VERSION << "\n"
        << "Contact: Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>\n"
        << "\n"
        << "Usage:   cnvetti <command> [options]\n"
        << "\n"
        << "Command: coverage   Compute coverage on single BAM file\n"
        << "         summaries  Summarize multi-sample BCF file\n"
        << "\n"
        << "Note: Use `bcftools merge` for merging multiple BCF files resulting from the\n"
        << "      coverage step for the input of the `cnvetti summaries` step.\n"
        << "\n";
}

// ----------------------------------------------------------------------------
// Function parseTopLevelCommandLine()
// ----------------------------------------------------------------------------

CnvettiCommand parseTopLevelCommandLine(int argc, char ** argv)
{
    if (argc < 2)
        return CnvettiCommand::NONE;

    std::string argv1(argv[1]);

    if (argv1 == "--help")
        return CnvettiCommand::TOPLEVEL_HELP;
    else if (argv1 == "--version")
        return CnvettiCommand::VERSION;
    else if (argv1 == "coverage")
        return CnvettiCommand::COVERAGE;
    else if (argv1 == "summaries")
        return CnvettiCommand::SUMMARIES;
    else
        throw TCLAP::ArgException(std::string("Invalid <command>: ") + argv1, "<command>");
}

// ----------------------------------------------------------------------------
// Function parseCoverageCommandLine()
// ----------------------------------------------------------------------------

CnvettiCoverageOptions parseCoverageCommandLine(int argc, char ** argv)
{
    CnvettiCoverageOptions result;

    return result;
}
