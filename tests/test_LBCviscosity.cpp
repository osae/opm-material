// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This test makes sure that mandated API is adhered to by all component classes
 */
#include "config.h"

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include "checkComponent.hpp"

#include <opm/material/fluidsystems/DecaneCO2FluidSystem.hpp>

#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include <iostream>
#include <iomanip>

namespace Opm {
namespace ComponentsTest {
#include <opm/material/components/co2tables.inc>
}}

#include <dune/common/parallel/mpihelper.hh>

template <class Scalar, class Evaluation, bool useLBCmod>
void testDecaneCO2()
{
    //These tables have been verified useng an independent implementation
    const std::vector<Scalar> lbcModOil = {3.64E-04, 2.72E-04, 1.97E-04, 3.67E-04, 2.76E-04, 2.01E-04,
        3.70E-04, 2.80E-04, 2.05E-04, 3.73E-04, 2.83E-04, 2.08E-04, 3.61E-04, 2.69E-04, 1.93E-04,
        3.64E-04, 2.73E-04, 1.97E-04, 3.67E-04, 2.76E-04, 2.01E-04, 3.70E-04, 2.80E-04, 2.05E-04,
        3.58E-04, 2.65E-04, 1.90E-04, 3.61E-04, 2.69E-04, 1.94E-04, 3.64E-04, 2.73E-04, 1.98E-04,
        3.67E-04, 2.76E-04, 2.01E-04, 3.54E-04, 2.62E-04, 1.86E-04, 3.58E-04, 2.65E-04, 1.90E-04,
        3.61E-04, 2.69E-04, 1.94E-04, 3.64E-04, 2.73E-04, 1.98E-04, 3.51E-04, 2.58E-04, 1.82E-04,
        3.54E-04, 2.62E-04, 1.86E-04, 3.58E-04, 2.65E-04, 1.90E-04, 3.61E-04, 2.69E-04, 1.94E-04};
    const std::vector<Scalar> lbcModGas = {1.49E-05, 1.72E-05, 1.95E-05, 1.52E-05, 1.74E-05, 1.97E-05,
        1.55E-05, 1.77E-05, 1.99E-05, 1.60E-05, 1.80E-05, 2.01E-05, 1.46E-05, 1.69E-05, 1.91E-05,
        1.49E-05, 1.71E-05, 1.93E-05, 1.52E-05, 1.73E-05, 1.95E-05, 1.59E-05, 1.77E-05, 1.97E-05,
        1.43E-05, 1.65E-05, 1.87E-05, 1.46E-05, 1.67E-05, 1.89E-05, 1.50E-05, 1.70E-05, 1.91E-05,
        1.59E-05, 1.75E-05, 1.94E-05, 1.40E-05, 1.61E-05, 1.83E-05, 1.43E-05, 1.64E-05, 1.85E-05,
        1.48E-05, 1.67E-05, 1.87E-05, 1.61E-05, 1.73E-05, 1.91E-05, 1.37E-05, 1.58E-05, 1.79E-05,
        1.40E-05, 1.60E-05, 1.81E-05, 1.47E-05, 1.65E-05, 1.84E-05, 1.73E-05, 1.71E-05, 1.88E-05};
    const std::vector<Scalar> lbcStdOil = {4.41E-04, 3.19E-04, 2.21E-04, 4.45E-04, 3.23E-04, 2.26E-04,
        4.50E-04, 3.28E-04, 2.31E-04, 4.54E-04, 3.33E-04, 2.35E-04, 4.45E-04, 3.20E-04, 2.20E-04,
        4.50E-04, 3.25E-04, 2.25E-04, 4.54E-04, 3.30E-04, 2.30E-04, 4.59E-04, 3.35E-04, 2.35E-04,
        4.49E-04, 3.21E-04, 2.20E-04, 4.53E-04, 3.26E-04, 2.25E-04, 4.58E-04, 3.31E-04, 2.30E-04,
        4.62E-04, 3.36E-04, 2.35E-04, 4.51E-04, 3.21E-04, 2.18E-04, 4.56E-04, 3.26E-04, 2.23E-04,
        4.61E-04, 3.31E-04, 2.28E-04, 4.65E-04, 3.36E-04, 2.33E-04, 4.53E-04, 3.20E-04, 2.16E-04,
        4.58E-04, 3.25E-04, 2.22E-04, 4.63E-04, 3.30E-04, 2.27E-04, 4.67E-04, 3.35E-04, 2.32E-04};
    const std::vector<Scalar> lbcStdGas = {1.49E-05, 1.72E-05, 1.95E-05, 1.52E-05, 1.74E-05, 1.97E-05,
        1.55E-05, 1.77E-05, 1.99E-05, 1.60E-05, 1.80E-05, 2.01E-05, 1.46E-05, 1.69E-05, 1.91E-05,
        1.49E-05, 1.71E-05, 1.93E-05, 1.53E-05, 1.74E-05, 1.95E-05, 1.59E-05, 1.77E-05, 1.98E-05,
        1.43E-05, 1.65E-05, 1.87E-05, 1.46E-05, 1.67E-05, 1.89E-05, 1.51E-05, 1.71E-05, 1.91E-05,
        1.60E-05, 1.75E-05, 1.95E-05, 1.40E-05, 1.61E-05, 1.83E-05, 1.43E-05, 1.64E-05, 1.85E-05,
        1.49E-05, 1.68E-05, 1.88E-05, 1.62E-05, 1.73E-05, 1.92E-05, 1.37E-05, 1.58E-05, 1.79E-05,
        1.41E-05, 1.61E-05, 1.81E-05, 1.48E-05, 1.65E-05, 1.85E-05, 1.76E-05, 1.72E-05, 1.89E-05};

    const std::vector<Scalar>& lbcOil = useLBCmod ? lbcModOil : lbcStdOil;
    const std::vector<Scalar>& lbcGas = useLBCmod ? lbcModGas : lbcStdGas;

    typedef Opm::MathToolbox<Evaluation> EvalToolbox;

    // Decane -- CO2
    typedef Opm::DecaneCO2FluidSystem<Scalar, useLBCmod> FluidSystem;
    typename FluidSystem::template ParameterCache<typename FluidSystem::Scalar> paramCache;

    Opm::CompositionalFluidState<Scalar, FluidSystem> fs;

    FluidSystem::init(293.15,393.15,1e5,40e6);

    //Phases indices:
    const int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    const int gasPhaseIdx = FluidSystem::gasPhaseIdx;

    //Components indices:
    const int decaneIdx = FluidSystem::DecaneIdx;
    const int co2Idx = FluidSystem::CO2Idx;

    if (!EvalToolbox::isSame(FluidSystem::criticalTemperature(co2Idx), 304.1, /*tolerance=*/1e-1))
       throw std::logic_error("testDecaneCO2(): Deviation in CO2 critical temperature");
    if (!EvalToolbox::isSame(FluidSystem::criticalPressure(co2Idx)/1e6, 7.38, /*tolerance=*/1e-2))
       throw std::logic_error("testDecaneCO2(): Deviation in CO2 critical pressure");
    if (!EvalToolbox::isSame(FluidSystem::criticalDensity(co2Idx), 467.6, /*tolerance=*/1e-1))
       throw std::logic_error("testDecaneCO2(): Deviation in CO2 critical density");
    if (!EvalToolbox::isSame(FluidSystem::molarMass(co2Idx)*1000, 44, /*tolerance=*/1e-1))
       throw std::logic_error("testDecaneCO2(): Deviation in CO2 molar mass");

    if (!EvalToolbox::isSame(FluidSystem::criticalTemperature(decaneIdx), 617.7, /*tolerance=*/1e-1))
       throw std::logic_error("testDecaneCO2(): Deviation in decane critical temperature");
    if (!EvalToolbox::isSame(FluidSystem::criticalPressure(decaneIdx)/1e6, 2.11, /*tolerance=*/1e-2))
       throw std::logic_error("testDecaneCO2(): Deviation in decane critical pressure");
    if (!EvalToolbox::isSame(FluidSystem::criticalDensity(decaneIdx), 228.0, /*tolerance=*/1e-1))
       throw std::logic_error("testDecaneCO2(): Deviation in decane critical density");
    if (!EvalToolbox::isSame(FluidSystem::molarMass(decaneIdx)*1000, 142.3, /*tolerance=*/1e-1))
       throw std::logic_error("testDecaneCO2(): Deviation in decane molar mass");

    int idx = 0;
    const int numMF = 5;
    for (int iMF=0; iMF<numMF; ++iMF) {
        fs.setMoleFraction(oilPhaseIdx, decaneIdx, 0.99-iMF*0.02);
        fs.setMoleFraction(oilPhaseIdx, co2Idx, 0.01+iMF*0.02);

        fs.setMoleFraction(gasPhaseIdx, decaneIdx, 0.01+iMF*0.02);
        fs.setMoleFraction(gasPhaseIdx, co2Idx, 0.99-iMF*0.02);

        const int numP = 4;
        const Scalar p0=1.0e5;
        for (int iP = 0; iP < numP; ++iP) {
            fs.setPressure(gasPhaseIdx, p0+iP*1e6);
            fs.setPressure(oilPhaseIdx, p0+iP*1e6);

            const int numT = 3;
            const Scalar T0=293.15;
            for (int iT = 0; iT < numT; ++iT) {
                fs.setTemperature(T0+iT*50.0);

                paramCache.updatePhase(fs, oilPhaseIdx);
                paramCache.updatePhase(fs, gasPhaseIdx);

                fs.setDensity(oilPhaseIdx, FluidSystem::density(fs, paramCache, oilPhaseIdx));
                fs.setDensity(gasPhaseIdx, FluidSystem::density(fs, paramCache, gasPhaseIdx));

                fs.setViscosity(oilPhaseIdx, FluidSystem::viscosity(fs, paramCache, oilPhaseIdx));
                fs.setViscosity(gasPhaseIdx, FluidSystem::viscosity(fs, paramCache, gasPhaseIdx));

                if (!EvalToolbox::isSame(fs.viscosity(oilPhaseIdx), lbcOil[idx], /*tolerance=*/5e-5)) {
                   std::cout << "gasvisc: " << fs.viscosity(oilPhaseIdx) << " (" << lbcOil[idx] << ") iMF:" << iMF << " iP:" << iP << " iT:" << iT << std::endl;
                   const std::string modelMsg = (useLBCmod ? "modified LBC model" : "standard LBC model");
                   throw std::logic_error("testDecaneCO2(): Deviation in oil viscosity: " + modelMsg);
                }

                if (!EvalToolbox::isSame(fs.viscosity(gasPhaseIdx), lbcGas[idx], /*tolerance=*/5e-6)) {
                   std::cout << "gasvisc: " << fs.viscosity(gasPhaseIdx) << " (" << lbcGas[idx] << ") iMF:" << iMF << " iP:" << iP << " iT:" << iT << std::endl;
                   const std::string modelMsg = (useLBCmod ? "modified LBC model" : "standard LBC model");
                   throw std::logic_error("testDecaneCO2(): Deviation in gas viscosity: " + modelMsg);
                }
                ++idx;
            }
        }
    }
}

template <class Scalar, class Evaluation>
void testAllComponents()
{
}

template <class Scalar>
inline void testAll()
{
    typedef Opm::DenseAd::Evaluation<Scalar, 3> Evaluation;

    testDecaneCO2<Scalar, Evaluation, true>();
    testDecaneCO2<Scalar, Evaluation, false>();

}


int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    testAll<double>();
    //testAll<float>();

    return 0;
}
