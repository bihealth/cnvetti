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
#include <type_traits>

// ---------------------------------------------------------------------------
// Class Interval
// ---------------------------------------------------------------------------

struct Interval
{
    Interval() : beginPos(), endPos() {}
    Interval(int beginPos, int endPos) : beginPos(beginPos), endPos(endPos) {}

    int beginPos { 0 };
    int endPos { 0 };

    int overlapLength(Interval const & other) const
    {
        if (!overlapsWith(other))
        {
            return 0;
        }
        else
        {
            return std::min(endPos, other.endPos) - std::max(beginPos, other.beginPos);
        }
    }

    int overlapsWith(Interval const & other) const
    {
        return (other.beginPos < endPos && beginPos < other.endPos);
    }
};

// ---------------------------------------------------------------------------
// Class IntervalWithLabel
// ---------------------------------------------------------------------------

template <typename TLabel>
struct IntervalWithLabel
{
    int beginPos { 0 };
    int endPos { 0 };
    TLabel label;
};

namespace impl {

// ---------------------------------------------------------------------------
// Class IntervalTreeEntry
// ---------------------------------------------------------------------------

template <typename TInterval>
class IntervalTreeEntry
{
public:
    IntervalTreeEntry() : interval(), maxEnd()
    {}

    IntervalTreeEntry(TInterval const & interval, int maxEnd) : interval(interval), maxEnd(maxEnd)
    {}

    TInterval const & getInterval() const
    {
        return interval;
    }

    int getBeginPos() const
    {
        return interval.beginPos;
    }

    int getEndPos() const
    {
        return interval.endPos;
    }

    int getMaxEnd() const
    {
        return maxEnd;
    }

    int setMaxEnd(int maxEnd)
    {
        this->maxEnd = maxEnd;
    }

    bool allLeftOf(int pos) const
    {
        return (maxEnd < pos);
    }

    bool isLeftOf(int pos) const
    {
        return (interval.endPos <= pos);
    }

    bool isRightOf(int pos) const
    {
        return (pos < interval.beginPos);
    }

    bool contains(int pos) const
    {
        return (interval.beginPos <= pos && pos < interval.endPos);
    }

    bool overlapsWith(int beginPos, int endPos) const
    {
        return (beginPos < interval.endPos && interval.beginPos < endPos);
    }

private:
    TInterval interval;
    int maxEnd;
};

}  // namespace impl

// ---------------------------------------------------------------------------
// Class IntervalTree
// ---------------------------------------------------------------------------

template <typename TLabel>
class IntervalTree
{
public:
    typedef typename std::conditional<
        std::is_void<TLabel>::value,
        Interval,
        IntervalWithLabel<TLabel>>::type TInterval;

private:
    typedef impl::IntervalTreeEntry<TInterval> TEntry;

public:
    IntervalTree() {}

    template <typename TIt>
    IntervalTree(TIt itBegin, TIt itEnd)
    {
        // Fill entries array.
        entries.reserve(itEnd - itBegin);
        for (TIt it = itBegin; it != itEnd; ++it)
        {
            entries.push_back(TEntry(*it, it->endPos));
        }

        // Sort by begin position.
        std::sort(
            entries.begin(), entries.end(),
            [](TEntry const & lhs, TEntry const & rhs) {
                return (lhs.getBeginPos() < rhs.getBeginPos());
            });

        // Compute maxEnd properties.
        computeMaxEndProperties(0, entries.size());
    }

    std::vector<TInterval> findOverlapping(int beginPos, int endPos)
    {
        std::vector<TInterval> result;
        findOverlapping(beginPos, endPos, result);
        return result;
    }

    void findOverlapping(int beginPos, int endPos, std::vector<TInterval> & result)
    {
        result.clear();
        findOverlapping(beginPos, endPos, 0, entries.size(), entries.size() / 2, result);
    }

    bool empty() const { return entries.empty(); }

private:

    void findOverlapping(int beginPos, int endPos, size_t beginIdx, size_t endIdx,
                         size_t centerIdx, std::vector<TInterval> & result)
    {
        if (beginIdx >= endIdx)
        {
            return;  // base case
        }

        TEntry const & entry = entries[centerIdx];

        if (entry.allLeftOf(beginPos))
        {
            // posBegin is right of the rightmost point of any interval in this node
            return;
        }

        if (beginIdx < centerIdx)  // recurse left
        {
            findOverlapping(
                beginPos, endPos, beginIdx, centerIdx,
                beginIdx + (centerIdx - beginIdx) / 2, result);
        }

        if (entry.overlapsWith(beginPos, endPos))  // check this node
        {
            result.push_back(entry.getInterval());
        }

        if (entry.isRightOf(endPos - 1))
        {
            // Last interval entry is left of the start of the interval,
            // cannot go to the right.
            return;
        }

        if (centerIdx + 1 < endIdx)  // recurse right
        {
            findOverlapping(
                beginPos, endPos, centerIdx + 1, endIdx,
                (centerIdx + 1) + (endIdx - (centerIdx + 1)) / 2,
                result);
        }
    }

    int computeMaxEndProperties(size_t beginIdx, size_t endIdx)
    {
        if (beginIdx == endIdx)
        {
            return -1;
        }

        size_t centerIdx = (endIdx + beginIdx) / 2;
        TEntry & entry = entries[centerIdx];

        if (beginIdx + 1 != endIdx)
        {
            entry.setMaxEnd(std::max(
                    entry.getMaxEnd(),
                    std::max(computeMaxEndProperties(beginIdx, centerIdx),
                             computeMaxEndProperties(centerIdx + 1, endIdx))));
        }

        return entry.getMaxEnd();
    }

    std::vector<TEntry> entries;
};
