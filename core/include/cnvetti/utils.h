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

#include <string>
#include <stdexcept>

bool hasSuffix(std::string const & str, std::string const & suffix);

// ---------------------------------------------------------------------------
// Struct GenomicRegion
// ---------------------------------------------------------------------------

struct GenomicRegion
{
    std::string contig;
    int beginPos;
    int endPos;
    int rID;

    GenomicRegion() = default;

    GenomicRegion(std::string const & contig, int beginPos, int endPos) :
        contig(contig), beginPos(beginPos), endPos(endPos), rID(-1)
    {}

    GenomicRegion(std::string const & contig, int beginPos, int endPos, int rID) :
        contig(contig), beginPos(beginPos), endPos(endPos), rID(rID)
    {}

    std::string toString() const;
};


GenomicRegion parseGenomicRegion(std::string const & str);
