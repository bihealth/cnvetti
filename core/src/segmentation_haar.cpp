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

// Based on HaarSeg.c for which the license header is reproduced below

/*
 *    Copyright (C) 2008  Erez Ben-Yaacov
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    <http://www.gnu.org/licenses/>
 *
 */

#include "cnvetti/segmentation.h"

#include <cmath>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <armadillo>

#define DEBUG 0


namespace { // anonymous

arma::fvec haarConvolution(
    arma::fvec const & signal, arma::fvec const & weight, int stepHalfSize)
{
    int const signalSize = signal.n_elem;
    bool const haveWeight = (weight.n_elem != 0u);
    if (DEBUG > 1)
        std::cerr << signalSize << " ~~ " << signal << "\n";
    arma::fvec result(signalSize, arma::fill::zeros);  // prepare all-zero result

    if (stepHalfSize > signalSize)
    {
        std::cerr
            << "Error?: stepHalfSize(" << stepHalfSize << ") > signalSize ("
            << signalSize << ")\n";
        return result;
    }

    // Only required when weights are given
    float highWeightSum, highNonNormed, lowWeightSum, lowNonNormed;

    if (haveWeight)
    {
        // Initialize weight sums
        auto const weightTail = weight.tail_rows(weight.n_elem - stepHalfSize);
        auto const signalTail = signal.tail_rows(signal.n_elem - stepHalfSize);
        highWeightSum = accu(weightTail);
        highNonNormed = accu(weightTail * signalTail);
        // Circular padding
        lowWeightSum = highWeightSum;
        lowNonNormed = -highNonNormed;
    }

    // This loop is the hotspot
    for (int k = 1; k < signalSize; ++k)
    {
        int highEnd = k + stepHalfSize - 1;
        if (highEnd >= signalSize)
            highEnd = signalSize - 1 - (highEnd - signalSize);
        int lowEnd = k - stepHalfSize - 1;
        if (lowEnd < 0)
            lowEnd = -lowEnd - 1;

        if (!haveWeight)
        {
            result[k] = result[k - 1] + signal[highEnd] + signal[lowEnd] - 2 * signal[k - 1];
        }
        else
        {
            lowNonNormed += signal[lowEnd] * weight[lowEnd] - signal[k - 1] * weight[k - 1];
            highNonNormed += signal[highEnd] * weight[highEnd] - signal[k - 1] * weight[k - 1];
            lowWeightSum += weight[k - 1] - weight[lowEnd];
            highWeightSum += weight[highEnd] - weight[k - 1];
            result[k] = sqrt(stepHalfSize / 2) * (
                lowNonNormed / lowWeightSum + highNonNormed / highWeightSum);
        }
    }

    if (!haveWeight)
    {
        float const stepNorm = sqrt(2 * stepHalfSize);
        result.rows(1, signalSize - 2) /= stepNorm;
    }

    return result;
}


std::vector<unsigned long long> findLocalPeaks(arma::fvec const & signal)
{
    int minSuspect = -1;
    int maxSuspect = -1;
    std::vector<unsigned long long> result;

    for (int k = 1; k + 1 < signal.n_elem; ++k)  // first and last not considered
    {
        float sigPrev = signal[k - 1];
        float sigCurr = signal[k];
        float sigNext = signal[k + 1];
        if (sigCurr > 0)
        {
            // Look for local maxima
            if (sigCurr > sigPrev && sigCurr > sigNext)
            {
                result.push_back(k);
            }
            else if (sigCurr > sigPrev && sigCurr == sigNext)
            {
                maxSuspect = k;
            }
            else if (sigCurr == sigPrev && sigCurr > sigNext)
            {
                if (maxSuspect != -1)
                {
                    result.push_back(maxSuspect);
                    maxSuspect = -1;
                }
            }
            else if (sigCurr == sigPrev && sigCurr < sigNext)
            {
                maxSuspect = -1;
            }
        }
        else if (sigCurr < 0)
        {
            // Look for local maxima
            if (sigCurr < sigPrev && sigCurr < sigNext)
            {
                result.push_back(k);
            }
            else if (sigCurr < sigPrev && sigCurr == sigNext)
            {
                minSuspect = k;
            }
            else if (sigCurr == sigPrev && sigCurr < sigNext)
            {
                if (minSuspect != -1)
                {
                    result.push_back(minSuspect);
                    minSuspect = -1;
                }
            }
            else if (sigCurr == sigPrev && sigCurr > sigNext)
            {
                minSuspect = -1;
            }
        }
    }

    return result;
}


// Perform merging of breakpoints as described in HaarSeg paper
std::vector<size_t> mergeBreakpoints(
    std::vector<size_t> const & prevBreakpoints,
    arma::uvec const & newBreakpoints,
    size_t windowSize)
{
    if (newBreakpoints.is_empty())
        return prevBreakpoints;

    // Merge all addon items outside a window around each base item
    std::vector<size_t> result;

    auto itNew = newBreakpoints.begin();
    for (auto itPrev = prevBreakpoints.begin(); itPrev != prevBreakpoints.end(); ++itPrev)
    {
        for (/*nop*/; (itNew != newBreakpoints.end()) && (*itNew <= *itPrev + windowSize); ++itNew)
        {
            if (*itNew + windowSize < *itPrev)
                result.push_back(*itNew);
        }
        result.push_back(*itPrev);
    }

    for (/*nop*/; itNew != newBreakpoints.end(); ++itNew)
        result.push_back(*itNew);

    // Sort breakpoints by position and return
    std::sort(result.begin(), result.end());
    return result;
}


arma::fvec pulseConvolution(arma::fvec const & signal, int pulseSize)
{
    int signalSize = signal.n_elem;
    if (pulseSize > signalSize)
        throw std::runtime_error("pulseSize > signalSize");
    float pulseHeight = 1.0 / pulseSize;

    // Circular padding init
    arma::fvec result(signalSize, arma::fill::zeros);
    for (int k = 0; k < (pulseSize + 1) / 2; ++k)
        result[0] += signal[k];
    for (int k = 0; k < pulseSize / 2; ++k)
        result[0] += signal[k];
    result[0] *= pulseHeight;

    int n = 1;
    for (int k = pulseSize / 2; k < signalSize + (pulseSize / 2) - 1; ++k)
    {
        int tail = k - pulseSize;
        if (tail < 0)
            tail = -tail - 1;
        int head = k;
        if (head >= signalSize)
            head = signalSize - 1 - (head - signalSize);
        result[n] = result[n - 1] + ((signal[head] - signal[tail]) * pulseHeight);
        n += 1;
    }

    return result;
}


// Compute median of absolute values (use differences from mean for MAD)
float medAbs(arma::fvec const & vals)
{
    if (vals.is_empty())
        return 0.0;
    else
        return median(abs(vals)) / 0.6745;
}


// Cumulative density function of the normal distribution
float pnorm(float mean, float sd)
{
    return 0.5 * (1 + erf(mean / sd / sqrt(2)));
}


// Compute FDR threshold
float fdrThresh(arma::fvec const & x, float q, float stdev)
{
    int const M = x.n_elem;
    if (M < 2)
        return 0.0;

    arma::fvec m(M, arma::fill::zeros);
    std::iota(m.begin(), m.end(), 1);
    m /= M;

    arma::fvec xSorted(sort(abs(x), "descend"));
    arma::fvec p(xSorted);
    std::transform(  // like R pnorm
        p.begin(), p.end(), p.begin(),
        [stdev](float x) { return 2 * (1 - pnorm(x, stdev)); });

    // Get the largest index for which p <= m * q
    arma::umat isSmaller(p <= m * q);
    int idx = -1;
    auto ptr = isSmaller.begin();
    for (int i = 0; i < isSmaller.n_elem; ++i, ++ptr)
        if (*ptr)
            idx = i;
    if (idx != -1)
    {
        return xSorted[idx];
    }
    else
    {
        std::cerr
            << "No passing p-values: min p= " << p[0] << ", min m=" << m[0]
            << ", q=" << q << "\n";
        return xSorted[0] + 1e-16;  // ~= 2^-52, like MATLAB "eps"
    }

}


// Signal standard deviation estimates and data for non-stationary variance estimation
struct SigmaEstimates
{
    // Estimate of standard deviatiation in signal for peaks
    float peakEst;

    // Estimate of standard deviation in signal for noise, only used when using
    // the raw values for correcting non-stationary variance
    float noiseEst;

    // Threshold used for creating binary mask
    float const tNSV;

    // Mask where the raw values are < tNSV (b[n] in the HaarSeg paper)
    arma::fvec const varMask;
};


// Perform estimation of sigmas and compute other values used when raw values are available
// for compensation of non-stationary variance compensation
SigmaEstimates estimateSigmas(
    arma::fvec const & vals,
    std::vector<float> const * rawValsPtr)
{
    float const THRESH_NSV = 50.0;  // empirically using 50.0 as in HaarSeg paper
    arma::fvec diffVals(haarConvolution(vals, arma::fvec(), 1));

    float peakEst = 0.0;
    float noiseEst = 0.0;
    arma::fvec varMask;

    if (rawValsPtr)
    {
        arma::uvec const uMask = arma::fvec(*rawValsPtr) < THRESH_NSV;
        varMask = arma::fvec(uMask.n_elem, arma::fill::zeros);
        std::copy(uMask.begin(), uMask.end(), varMask.begin());

        int const pulseSize = 2;
        arma::uvec const diffMask = (pulseConvolution(varMask, pulseSize) >= 0.5);
        peakEst = medAbs(diffVals.elem(find(diffMask)));
        noiseEst = medAbs(diffVals.elem(find(diffMask == 0u)));
    }
    else
    {
        peakEst = medAbs(diffVals);
    }

    return SigmaEstimates { peakEst, noiseEst, THRESH_NSV, varMask };
}


// Perform compensation for non-stationary variance
void performNsvCompensation(SigmaEstimates & sigmaEst, arma::fvec & convRes, int stepHalfSize)
{
    int const pulseSize = 2 * stepHalfSize;
    arma::uvec const convMask = (pulseConvolution(sigmaEst.varMask, pulseSize) >= 0.5);
    arma::fvec fConvMask(convMask.n_elem, arma::fill::zeros);
    std::copy(convMask.begin(), convMask.end(), fConvMask.begin());

    arma::fvec sigmaEstVec =
        (1 - fConvMask) * sigmaEst.peakEst + fConvMask * sigmaEst.noiseEst;
    convRes /= sigmaEstVec;
    sigmaEst.peakEst = 1.0;
}

}  // anonymous namespace


std::vector<size_t> segmentHaarSeg(
    std::vector<float> & vals,
    float breaksFdrQ,
    std::vector<float> const * qualsPtr,
    std::vector<float> const * rawValsPtr,
    int haarStartLevel,
    int haarEndLevel)
{
    // Convert input std::vector of float values into an Armadillo float vector
    arma::fvec fvals(vals);

    // Estimate sigma and compute the data necessary for non-stationary variance compensation
    SigmaEstimates sigmaEst = estimateSigmas(fvals, rawValsPtr);

    // Create Armadillo float vector for quality values, empty if qualsPtr == nullptr
    arma::fvec quals;
    if (qualsPtr)
        quals = *qualsPtr;

    // 0-based positions with breakpoints, built incrementally
    std::vector<size_t> breakpoints;

    // Collect the breakpoints for each of the considered levels
    for (int level = haarStartLevel; level <= haarEndLevel; ++level)
    {
        int const stepHalfSize = (1 << level);
        // Perform Haar convolution and find local peaks
        arma::fvec convRes = haarConvolution(fvals, quals, stepHalfSize);
        arma::uvec const localPeaks = findLocalPeaks(convRes);
        if (DEBUG > 1)
            std::cerr << "Found " << localPeaks.n_elem << " peaks at level " << level << "\n";

        // Perform compenstation of non-stationary variance if raw values have been given
        if (rawValsPtr)
            performNsvCompensation(sigmaEst, convRes, stepHalfSize);

        // Compute FDR threshold
        float const T = fdrThresh(convRes.elem(localPeaks), breaksFdrQ, sigmaEst.peakEst);

        // Keep only the peak values where the signal amplitude is large enough
        arma::fvec localPeaksConvRes(convRes.elem(localPeaks));
        arma::uvec addonPeaks(localPeaks.elem(find(abs(localPeaksConvRes) >= T)));
        breakpoints = mergeBreakpoints(breakpoints, addonPeaks, (1 << (level - 1)));
        if (DEBUG > 0)
        {
            std::cerr << "level " << level << "\n\nbreakpoints\n";
            for (int bp : breakpoints)
                std::cerr << bp << "\n";
            std::cerr << "\n";
        }
    }

    if (DEBUG > 1)
    {
        std::cerr
            << "Found " << breakpoints.size() << " breakpoints \n"
            << "Breakpoints\n";
        for (auto bp : breakpoints)
            std::cerr << "\t" << bp << "\n";
        std::cerr << "(#vals = " << vals.size() << ")\n";
    }

    // Add start and end marker to breakpoints, for writing out segmented values
    breakpoints.insert(breakpoints.begin(), 0u);
    breakpoints.push_back(vals.size());

    return breakpoints;
}
