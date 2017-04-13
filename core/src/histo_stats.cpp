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

#include "cnvetti/histo_stats.h"

#include <cmath>
#include <map>
#include <vector>

// ---------------------------------------------------------------------------
// Class HistoStatsHandlerImpl
// ---------------------------------------------------------------------------

class HistoStatsHandlerImpl
{
public:
    HistoStatsHandlerImpl(unsigned sampleCount)
    {
        // Pre-allocate the histograms (histograms are ordered std::map)
        histograms.resize(sampleCount);
        for (auto & forSample : histograms)
        {
            forSample.resize(static_cast<int>(StatsMetric::MAX) + 1);
            for (auto & forMetric : forSample)
                forMetric.resize(101);
        }
    }

    void registerValue(
        unsigned sampleID, StatsMetric metric, int gcContent, double value)
    {
        histograms.at(sampleID).at(static_cast<int>(metric)).at(gcContent)[
            int(round(value * factor))] += 1;
    }

    double getMedian(unsigned sampleID, StatsMetric metric, int gcContent) const
    {
        std::map<int, int> const & histo = histograms.at(
            sampleID).at(static_cast<int>(metric)).at(gcContent);
        // Determine total number of value registered in histogram
        int totalCount = 0;
        for (auto const & pair : histo)
            totalCount += pair.second;
        if (totalCount == 0u)
            return 0.0;  // Short-circuit if empty

        // Compute median
        int low;
        int up;
        if (totalCount % 2 == 0)
        {
            low = (totalCount / 2) - 1;
            up = low + 1;
        }
        else
        {
            low = totalCount / 2;
            up = low;
        }

        int a = -1, b = -1;
        for (auto const & pair : histo)
        {
            int const val = pair.first;
            int const num = pair.second;
            if (a == -1 && num > low)
                a = val;
            else if (a == -1 && num <= low)
                low -= num;
            if (b == -1 && num > up)
                b = val;
            else if (b == -1 && num <= up)
                up -= num;
            if (a != -1 && b != -1)
                break;
        }

        if (a == -1 && b == -1)
            return 0.0;
        else
            return static_cast<double>(a + b) / 2.0 / factor;
    }

    int getNumWindows(int gcContent) const
    {
        int result = 0;
        std::map<int, int> const & histo = histograms.at(
            0).at(static_cast<int>(StatsMetric::COV)).at(gcContent);
        for (auto const & pair : histo)
            result += pair.second;
        return result;
    }

private:
    // histograms[sampleID][metric][gcContent][value] = count
    std::vector<
        std::vector<
            std::vector<
                std::map<int, int>
            >
        >
    > histograms;

    // Factor to use for rounding (number of significant figures)
    int factor { 100 };
};


// ---------------------------------------------------------------------------
// Class HistoStatsHandler
// ---------------------------------------------------------------------------

HistoStatsHandler::HistoStatsHandler(unsigned sampleCount) :
    impl(new HistoStatsHandlerImpl(sampleCount))
{}

HistoStatsHandler::~HistoStatsHandler() = default;

void HistoStatsHandler::registerValue(
    unsigned sampleID, StatsMetric metric, int gcContent, double value)
{
    impl->registerValue(sampleID, metric, gcContent, value);
}

double HistoStatsHandler::getMedian(
    unsigned sampleID, StatsMetric metric, int gcContent) const
{
    return impl->getMedian(sampleID, metric, gcContent);
}

int HistoStatsHandler::getNumWindows(int gcContent) const
{
    return impl->getNumWindows(gcContent);
}
