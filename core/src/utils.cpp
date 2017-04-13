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

#include "cnvetti/utils.h"

#include <sstream>

bool hasSuffix(std::string const & str, std::string const & suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string GenomicRegion::toString() const
{
    std::stringstream ss;
    ss << contig << ":" << (beginPos + 1) << "-" << endPos;
    return ss.str();
}

GenomicRegion parseGenomicRegion(std::string const & str)
{
    std::string sContig, sBegin, sEnd;
    enum { CONTIG, BEGIN, END } state = CONTIG;
    for (char c : str) {
        switch (state) {
            case CONTIG:
                if (c != ':')
                    sContig.push_back(c);
                else
                    state = BEGIN;
                break;
            case BEGIN:
                if (c == '-')
                    state = END;
                else if (c == ',')
                    continue;  // ignore
                else if (isdigit(c))
                    sBegin.push_back(c);
                else
                        throw std::runtime_error(std::string("Invalid region: ") + str);
                break;
            case END:
                if (c == ',')
                    continue;  // ignore
                else if (isdigit(c))
                    sEnd.push_back(c);
                else
                        throw std::runtime_error(std::string("Invalid region: ") + str);
                break;
        }
    }
    if (sContig.empty() || sBegin.empty() || sEnd.empty())
        throw std::runtime_error(std::string("Invalid region: ") + str);
    return GenomicRegion(sContig, atoi(sBegin.c_str()) - 1, atoi(sEnd.c_str()));
}
