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

#include <cstdlib>
#include <iostream>
#include <random>

#include "cnvetti/segmentation.h"

int const SEED = 41;

int main(int argc, char ** argv)
{
    std::vector<float> origVals(1500, 0.0);
    std::fill(origVals.begin() + 500, origVals.begin() + 1000, 1.0);

    std::mt19937 gen(SEED);
    float noise = 0.4;
    if (argc >= 2)
        noise = atof(argv[1]);
    std::normal_distribution<> d(0.0, noise);
    std::vector<float> noiseVals;
    for (auto & x : origVals)
    {
        noiseVals.push_back(d(gen));
        x += noiseVals.back();
    }

    float breaksFdrQ = 1e-3;
    if (argc >= 3)
        breaksFdrQ = atof(argv[2]);

    int lMin = 1;
    if (argc >= 4)
        lMin = atoi(argv[3]);
    int lMax = 5;
    if (argc >= 5)
        lMax = atoi(argv[4]);

    float denoise = 0;
    if (argc >= 6)
        denoise = atof(argv[5]);

    std::vector<size_t> breakpoints = segmentHaarSeg(origVals, breaksFdrQ, NULL, NULL, lMin, lMax);

    if (denoise != 0)
    {
        std::cerr << "DENOISING DATA\n";
        breakpoints = refineSegmentationShiftBreakpoints(origVals, breakpoints);
        breakpoints = refineSegmentationDeleteBreakpoints(origVals, breakpoints, 0.01);
    }

    std::vector<float> segVals = replaceWithSegmentMedians(origVals, breakpoints);

    std::cout << "pos\tvalue\tnoise\tsegmented\n";
    int pos = 0;
    for (auto it = origVals.begin(), it2 = noiseVals.begin(), it3 = segVals.begin();
         it != origVals.end(); ++it, ++it2, ++it3, ++pos)
        std::cout << pos << '\t' << *it << '\t' << *it2 << '\t' << *it3 << '\n';

    return 0;
}
