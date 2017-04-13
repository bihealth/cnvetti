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
