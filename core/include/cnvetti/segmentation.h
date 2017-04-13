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

#include <vector>
#include <cstddef>


// HaarSeq segmentation algorithm
//
// Single outliers that have a chi square value > outlierThresh are ignored
// for breakpoint computation
std::vector<size_t> segmentHaarSeg(
    std::vector<float> & vals,
    float breaksFdrQ,
    std::vector<float> const * quals = nullptr,
    std::vector<float> const * rawI = nullptr,
    int haarStartLevel = 1,
    int haarEndLevel = 5
);


// Replace values by segment medians
std::vector<float> replaceWithSegmentMedians(
    std::vector<float> const & values,
    std::vector<size_t> const & breakpoints);


// Segmentation refinement based on shifting breakpoints to reduce overall
// segment error.
std::vector<size_t> refineSegmentationShiftBreakpoints(
    std::vector<float> const & values,
    std::vector<size_t> const & breakpoints);


// Segmentation refinement based on deleting breakpoints, affected segments'
// total error must be less than the given mergeDelta times original error
std::vector<size_t> refineSegmentationDeleteBreakpoints(
    std::vector<float> const & values,
    std::vector<size_t> breakpoints,
    float mergeDelta = 0.01);
