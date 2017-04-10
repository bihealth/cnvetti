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

#include "cnvetti/histo_stats.h"

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
        unsigned sampleID, StatsMetric metric, int gcContent, int value)
    {
        histograms.at(sampleID).at(static_cast<int>(metric)).at(gcContent)[value] += 1;
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
            return static_cast<double>(a + b) / 2.0;
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
};


// ---------------------------------------------------------------------------
// Class HistoStatsHandler
// ---------------------------------------------------------------------------

HistoStatsHandler::HistoStatsHandler(unsigned sampleCount) :
    impl(new HistoStatsHandlerImpl(sampleCount))
{}

HistoStatsHandler::~HistoStatsHandler() = default;

void HistoStatsHandler::registerValue(
    unsigned sampleID, StatsMetric metric, int gcContent, int value)
{
    impl->registerValue(sampleID, metric, gcContent, value);
}

double HistoStatsHandler::getMedian(
    unsigned sampleID, StatsMetric metric, int gcContent) const
{
    return impl->getMedian(sampleID, metric, gcContent);
}
