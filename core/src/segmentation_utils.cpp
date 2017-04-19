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

#include "cnvetti/segmentation.h"

#include <algorithm>
#include <limits>

#include <armadillo>

#define DEBUG 1


std::vector<float> replaceWithSegmentMedians(
    std::vector<float> const & values,
    std::vector<size_t> const & breakpoints)
{
    // Replace segment values by their medians
    std::vector<float> result;
    if (values.empty())
        return result;

    arma::fvec fvals(values);

    for (int i = 0; i + 1 < breakpoints.size(); ++i)
        std::fill_n(
            std::back_inserter(result),
            breakpoints[i + 1] - breakpoints[i],
            median(fvals.rows(breakpoints[i], breakpoints[i + 1] - 1)));

    return result;
}


std::vector<size_t> refineSegmentationShiftBreakpoints(
    std::vector<float> const & values,
    std::vector<size_t> const & breakpoints)
{
    return breakpoints;
}

namespace {  // anonymous

struct HeapElem
{
    HeapElem() = default;
    HeapElem(float oldVal, float newVal, size_t breakpoint) :
        oldVal(oldVal), newVal(newVal), breakpoint(breakpoint)
    {}

    float oldVal { 0.0 };
    float newVal { 0.0 };
    size_t breakpoint { 0 };

    bool operator<(HeapElem const & other) const
    {
        return (weight() < other.weight());
    }

private:

    float weight() const
    {
        return (oldVal / newVal);
    }
};

size_t INVALID = std::numeric_limits<size_t>::max();


std::pair<float, float>
computeError(
    arma::fvec const & values, arma::fvec const & errorSquares,
    size_t bpPrev, size_t bpThis, size_t bpNext)
{
    size_t const leftLen = bpThis - bpPrev;
    size_t const rightLen = bpNext - bpThis;

    // Compute old error
    float const oldLeft = sqrt(accu(errorSquares.rows(bpPrev, bpThis - 1)) / leftLen);
    float const oldRight = sqrt(accu(errorSquares.rows(bpThis, bpNext - 1)) / rightLen);
    float const oldVal = (oldLeft + oldRight) / 2.0;

    // Compute potential change for selecting breakpoint at idx for deletion
    float const newMedian = median(values.rows(bpPrev, bpNext - 1));
    float const newLeft = sqrt(accu(square(values.rows(bpPrev, bpThis - 1) - newMedian)) / leftLen);
    float const newRight = sqrt(accu(square(values.rows(bpThis, bpNext - 1) - newMedian))  / rightLen);
    float const newVal = (newLeft + newRight) / 2.0;

    return std::make_pair(oldVal, newVal);
}

}  // anonymous namespace


std::vector<size_t> refineSegmentationDeleteBreakpoints(
    std::vector<float> const & values,
    std::vector<size_t> breakpoints,
    float mergeDelta)
{
    if (DEBUG)
    {
        std::cerr << "BREAKPOINTS\n";
        for (size_t i = 0; i < breakpoints.size(); ++i)
            std::cerr << i << "\t" << breakpoints[i] << "\n";
        std::cerr << "\n";
    }

    arma::fvec fValues(values);
    arma::fvec segMedians(replaceWithSegmentMedians(values, breakpoints));
    arma::fvec errorSquares = square(fValues - segMedians);

    // Mapping from input breakpoint index to output, INVALID if already removed.  This is the
    // simplest way to emulate a heap with "DECREATE-KEY()" using the STL.
    std::vector<size_t> bpIdx(breakpoints.size());
    std::iota(bpIdx.begin(), bpIdx.end(), 0);

    // Max-heap that stores candidate breakpoints to delete
    std::vector<HeapElem> heap;

    // Initialize heap by computing the effect in error change
    for (size_t idx = 1; idx + 1 < bpIdx.size(); ++idx)
    {
        size_t const leftLen = breakpoints[idx] - breakpoints[idx - 1];
        size_t const rightLen = breakpoints[idx + 1] - breakpoints[idx];

        // Compute old error
        float const oldLeft = sqrt(accu(errorSquares.rows(breakpoints[idx - 1], breakpoints[idx] - 1)) / leftLen);
        float const oldRight = sqrt(accu(errorSquares.rows(breakpoints[idx], breakpoints[idx + 1] - 1)) / rightLen);
        float const oldVal = (oldLeft + oldRight) / 2.0;

        // Compute potential change for selecting breakpoint at idx for deletion
        float const newMedian = median(fValues.rows(breakpoints[idx - 1], breakpoints[idx + 1] - 1));
        float const newLeft = sqrt(accu(square(fValues.rows(breakpoints[idx - 1], breakpoints[idx] - 1) - newMedian)) / leftLen);
        float const newRight = sqrt(accu(square(fValues.rows(breakpoints[idx], breakpoints[idx + 1] - 1) - newMedian))  / rightLen);
        float const newVal = (newLeft + newRight) / 2.0;

        heap.push_back(HeapElem(oldVal, newVal, idx));
        if (DEBUG)
            std::cerr
                << "HeapElem(" << heap.back().oldVal << ", " << heap.back().newVal << ", "
                << heap.back().breakpoint << ")\n\tnewMedian=" << newMedian << ", newLeft="
                << newLeft << ", newRight=" << newRight << "; oldLeft=" << oldLeft
                << ", oldRight=" << oldRight << "\n";
        std::push_heap(heap.begin(), heap.end());
    }

    // Consider all elements in heap, if removing the currently best candidate is valid then
    // do it, otherwise skip
    while (!heap.empty())
    {
        // Pop best candidate from heap
        std::pop_heap(heap.begin(), heap.end());
        HeapElem const elem = heap.back();
        heap.pop_back();

        if (DEBUG)
            std::cerr
                << "Considering HeapElem(" << elem.oldVal << ", " << elem.newVal << ", "
                << elem.breakpoint << "); " << bpIdx[elem.breakpoint] << "\n"
                << "elem.newVal = " << elem.newVal << ", elem.oldVal  = "
                << elem.oldVal << ", 1 + mergeDelta = " << (1 + mergeDelta) << ", oldVal / newVal = "
                << (elem.oldVal / elem.newVal) << ")\n";

        // Skip if breakpoint is already removed
        if (bpIdx[elem.breakpoint] == INVALID)
            continue;
        // Break if error increase would be too large, no more suitable candidates
        if (elem.newVal / elem.oldVal > (1.0 + mergeDelta))
            break;

        // Store indices in breakpoint for the breakpoints left and right of to be removed
        // breakpoint for later
        size_t bpLeft = elem.breakpoint - 1;
        while (bpLeft > 0u && bpIdx[bpLeft] == INVALID)
            bpLeft--;
        size_t bpRight = elem.breakpoint + 1;
        while (bpRight + 1 < bpIdx.size() && bpIdx[bpRight] == INVALID)
            bpRight++;

        // Remove breakpoint, shift indices of right breakpoints to the left
        if (DEBUG)
            std::cerr << "Removing breakpoint " << breakpoints[bpIdx[elem.breakpoint]] << "\n";
        breakpoints.erase(breakpoints.begin() + bpIdx[elem.breakpoint]);
        for (auto it = bpIdx.begin() + elem.breakpoint; it != bpIdx.end(); ++it)
            if (*it != INVALID)
                *it -= 1;
        bpIdx[elem.breakpoint] = INVALID;

        // Update median and standard error of the new segment
        float newMedian = median(fValues.rows(breakpoints[bpIdx[bpLeft]], breakpoints[bpIdx[bpRight]] - 1));
        if (DEBUG)
            std::cerr << "newMedian = " << newMedian << "\n";
        std::fill(
            segMedians.begin() + breakpoints[bpIdx[bpLeft]],
            segMedians.begin() + breakpoints[bpIdx[bpRight]],
            newMedian);
        errorSquares.rows(breakpoints[bpIdx[bpLeft]], breakpoints[bpIdx[bpRight]] - 1) =
            square(fValues.rows(breakpoints[bpIdx[bpLeft]], breakpoints[bpIdx[bpRight]] - 1) -
                   segMedians.rows(breakpoints[bpIdx[bpLeft]], breakpoints[bpIdx[bpRight]] - 1));

        // Push breakpoints bpLeft and bpRight into heap, except if first/last
        if (bpIdx[bpLeft] > 0u && bpIdx[bpLeft] + 1 < breakpoints.size())
        {
            std::pair<float, float> error = computeError(
                fValues, errorSquares, breakpoints[bpIdx[bpLeft] - 1], breakpoints[bpIdx[bpLeft]],
                breakpoints[bpIdx[bpLeft] + 1]);
            heap.push_back(HeapElem(error.first, error.second, bpLeft));
            if (DEBUG)
                std::cerr
                    << "pushing new LEFT element HeapElement(" << heap.back().oldVal << ", " << heap.back().newVal
                    << ", " << heap.back().breakpoint << ")\n";
            push_heap(heap.begin(), heap.end());
        }
        if (bpIdx[bpRight] + 1 < breakpoints.size())
        {
            std::pair<float, float> error = computeError(
                fValues, errorSquares, breakpoints[bpIdx[bpRight] - 1], breakpoints[bpIdx[bpRight]],
                breakpoints[bpIdx[bpRight] + 1]);
            heap.push_back(HeapElem(error.first, error.second, bpRight));
            if (DEBUG)
                std::cerr
                    << "pushing new RIGHT element HeapElement(" << heap.back().oldVal << ", " << heap.back().newVal
                    << ", " << heap.back().breakpoint << ")\n";
            push_heap(heap.begin(), heap.end());
        }
        if (DEBUG)
            std::cerr << "\n";
    }

    return breakpoints;
}
