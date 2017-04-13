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

#include <memory>
#include <vector>

// ---------------------------------------------------------------------------
// Forward Declarations
// ---------------------------------------------------------------------------

class HistoStatsHandlerImpl;

// ---------------------------------------------------------------------------
// Enum StatsMetric
// ---------------------------------------------------------------------------

enum class StatsMetric
{
    COV = 0,
    COV0 = 1,
    RC = 2,
    RC0 = 3,
    MAX = 3
};

// ---------------------------------------------------------------------------
// Class HistoStatsHandler
// ---------------------------------------------------------------------------

class HistoStatsHandler
{
public:
    HistoStatsHandler(unsigned sampleCount);
    ~HistoStatsHandler();

    void registerValue(unsigned sampleID, StatsMetric metric, int gcContent, double value);

    double getMedian(unsigned sampleID, StatsMetric metric, int gcContent) const;

    int getNumWindows(int gcContent) const;

private:
    std::unique_ptr<HistoStatsHandlerImpl> impl;
};
