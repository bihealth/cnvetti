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

#include "cnvetti/contig_selection.h"

#include <algorithm>
#include <vector>

namespace {  // anonymous namespace

// TODO: Make configurable
static std::vector<std::string> CONTIG_WHITELIST {
    "1", "chr1",
    "2", "chr2",
    "3", "chr3",
    "4", "chr4",
    "5", "chr5",
    "6", "chr6",
    "7", "chr7",
    "8", "chr8",
    "9", "chr9",
    "10", "chr10",
    "11", "chr11",
    "12", "chr12",
    "13", "chr13",
    "14", "chr14",
    "15", "chr15",
    "16", "chr16",
    "17", "chr17",
    "18", "chr18",
    "19", "chr19",
    "20", "chr20",
    "21", "chr21",
    "22", "chr22",
    "X", "chrX",
    "Y", "chrY"
};

}  // anonymous namespace

bool contigWhitelisted(std::string contig)
{
    return std::find(
        CONTIG_WHITELIST.begin(), CONTIG_WHITELIST.end(), contig) != CONTIG_WHITELIST.end();
}
