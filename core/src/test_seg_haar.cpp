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

int const SEED = 42;

int main(int argc, char ** argv)
{
    std::vector<float> origVals(1500, 0.0);
    std::fill(origVals.begin() + 500, origVals.begin() + 1000, 1.0);

    std::mt19937 gen(SEED);
    float noise = 0.4;
    if (argc >= 2)
        noise = atof(argv[1]);
    std::normal_distribution<> d(0.0, noise);
    for (auto & x : origVals)
        x += d(gen);

    float breaksFdrQ = 1e-3;
    if (argc >= 3)
        breaksFdrQ = atof(argv[2]);

    int lMin = 1;
    if (argc >= 4)
        lMin = atoi(argv[3]);
    int lMax = 5;
    if (argc >= 5)
        lMax = atoi(argv[4]);

    std::vector<float> segVals = segmentHaarSeg(origVals, breaksFdrQ, NULL, NULL, lMin, lMax);

    std::cerr
        << "noise = " << noise << "\n"
        << "breaks FDR Q = " << breaksFdrQ << "\n"
        << "lMin = " << lMin << "\n"
        << "lMax = " << lMax << "\n";

    int pos = 0;
    for (auto it = origVals.begin(), it2 = segVals.begin(); it != origVals.end(); ++it, ++it2, ++pos)
        std::cout << pos << '\t' << *it << '\t' << *it2 << '\n';

    return 0;
}
