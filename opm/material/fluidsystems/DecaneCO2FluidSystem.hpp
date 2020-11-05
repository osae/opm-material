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
 * \copydoc Opm::DecaneCO2FluidSystem
 */
#ifndef DECANE_CO2_FLUIDSYSTEM_HH
#define DECANE_CO2_FLUIDSYSTEM_HH

#include <iostream>
#include <cassert>
#include <stdexcept>  // invalid_argument
#include <sstream>
#include <iostream>
#include <string>
#include <random>    // mt19937, normal_distribution
#include <limits>    // epsilon
#include <boost/format.hpp>  // boost::format

#include <opm/common/Exceptions.hpp>
#include <opm/material/IdealGas.hpp>

#include <opm/material/components/Component.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/Brine.hpp>
#include <opm/material/eos/PengRobinsonMixture.hpp>
#include <opm/material/eos/PengRobinsonParamsMixture.hpp>

#include <opm/material/fluidsystems/DecaneCO2ParameterCache.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LBCviscosity.hpp>
//#include "ChiParameterCache.hpp"

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>

namespace Opm {

template <class Scalar>
class ChiwomsCO2 : public Opm::SimpleCO2<Scalar>
{
public:
    /// Acentric factor
    static Scalar acentricFactor() { return 0.225; }

    /// Critical density
    static Scalar criticalDensity() { return 467.6; }
};

template <class Scalar>
class Decane : public Opm::Component<Scalar, Decane<Scalar> >
{
public:
        /// copied from cool probleps
        /// Chemical name
        static const char* name() { return "C10"; }

        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.1422817; }

        /// Critical temperature in \f$\mathrm[K]}\f$
        static Scalar criticalTemperature() { return 617.7; }

        /// Critical pressure in \f$\mathrm[Pa]}\f$
        static Scalar criticalPressure() { return 2.11e6; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.4884; }

        /// Critical density
        static Scalar criticalDensity() { return 228; }
};
/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with brine and decane as the main components
 * in each their phase, and CO2 as solvent in both.
 */
template <class Scalar>
class DecaneCO2FluidSystem
        : public Opm::BaseFluidSystem<Scalar, DecaneCO2FluidSystem<Scalar> >
{
    typedef DecaneCO2FluidSystem<Scalar> ThisType;
    typedef Opm::BaseFluidSystem<Scalar, ThisType> Base;
    typedef typename Opm::PengRobinson<Scalar> PengRobinson;
    typedef typename Opm::PengRobinsonMixture<Scalar, ThisType> PengRobinsonMixture;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    
    //template <class Evaluation>
    //using ParameterCache = Opm::ChiParameterCache<Evaluation, ThisType>;    
    template <class Evaluation>
    using ParameterCache = Opm::DecaneCO2ParameterCache<Evaluation, ThisType>;

    /****************************************
     * Fluid phase related static parameters
    ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the liquid phase
    static const int oilPhaseIdx = 0;
    static const int gasPhaseIdx = 1;
    static const int waterPhaseIdx = -17; //osae: hack for Spe5Cache ...

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {"o",  // oleic phase
                                    "g"};  // gas phase

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned phaseIdx)
    {
        if (phaseIdx == oilPhaseIdx)
            return true;

        // CO2 have associative effects with Decane
        return true;
    }


    /****************************************
    * Component related static parameters
    ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;  // Decane, co2

    //! The component index of the oil; Decane
    static const int DecaneIdx = 0;

    //! The component index of the solvent; co2
    static const int CO2Idx = 1;

    //! The component for pure oil
    typedef Opm::Decane<Scalar> Decane;

    //! The component for pure solvent
    typedef Opm::ChiwomsCO2<Scalar> CO2;


    static void init(Scalar minT = 273.15,
                     Scalar maxT = 373.15,
                     Scalar minP = 1e4,
                     Scalar maxP = 100e6)
    {
        Opm::PengRobinsonParamsMixture<Scalar, ThisType, oilPhaseIdx, /*useSpe5=*/true> prParams;

        // find envelopes of the 'a' and 'b' parameters for the range
        // minT <= T <= maxT and minP <= p <= maxP. For
        // this we take advantage of the fact that 'a' and 'b' for
        // mixtures is just a convex combination of the attractive and
        // repulsive parameters of the pure components

        Scalar minA = 1e30, maxA = -1e30;
        Scalar minB = 1e30, maxB = -1e30;

        prParams.updatePure(minT, minP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, minP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(minT, maxP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, maxP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };
        PengRobinson::init(/*aMin=*/minA, /*aMax=*/maxA, /*na=*/100,
                           /*bMin=*/minB, /*bMax=*/maxB, /*nb=*/200);
    }

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
            static const char* name[] = {
                    Decane::name(),
                    CO2::name()
            };
            assert(0 <= compIdx && compIdx < numComponents);
            return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        return (compIdx == DecaneIdx)
            ? Decane::molarMass()
            : (compIdx == CO2Idx)
            ? CO2::molarMass()
            : throw std::invalid_argument("Molar mass component index");
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
            return (compIdx == DecaneIdx)
                    ? Decane::criticalTemperature()
                    : (compIdx == CO2Idx)
                    ? CO2::criticalTemperature()
                    : throw std::invalid_argument("Critical temperature component index");
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
            return (compIdx == DecaneIdx)
                    ? Decane::criticalPressure()
                    : (compIdx == CO2Idx)
                    ? CO2::criticalPressure()
                    : throw std::invalid_argument("Critical pressure component index");
    }

    static Scalar criticalDensity(unsigned compIdx)
    {
        return (compIdx == DecaneIdx)
                ? Decane::criticalDensity()
                : (compIdx == CO2Idx)
                ? CO2::criticalDensity()
                : throw std::invalid_argument("Critical density component index");
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
            return (compIdx == DecaneIdx)
                    ? Decane::acentricFactor()
                    : (compIdx == CO2Idx)
                    ? CO2::acentricFactor()
                    : throw std::invalid_argument("Molar mass component index");
    }

    /*!
     * \copydoc BaseFluidSystem::isLiquid
     */
    static bool isLiquid(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx != gasPhaseIdx;
    }

    /*!
     * \copydoc BaseFluidSystem::isIdealGas
     */
    static bool isIdealGas(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == gasPhaseIdx)
            return CO2::gasIsIdeal();
        return false;
    }

    /*!
     * \copydoc BaseFluidSystem::isCompressible
     */
    static bool isCompressible(unsigned phaseIdx OPM_OPTIM_UNUSED)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }


    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::density
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& paramCache,
                           unsigned phaseIdx)
    {
        if (false)// set to true if you want constant density
        {
            if(phaseIdx == oilPhaseIdx) {
                return 670; 
            } else {
                return 1.7; 
            }
        } else {
            LhsEval result = Opm::decay<LhsEval>(fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx));
            return result;
        }
    }

        //! \copydoc BaseFluidSystem::viscosity
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval viscosity(const FluidState& fluidState,
                                 const ParameterCache<ParamCacheEval>& paramCache,
                                 unsigned phaseIdx)
        {
            assert(0 <= phaseIdx && phaseIdx < numPhases);

            //#warning We use constant viscosity. These needs to be checked
            //#warning We use the same as for octane
            //std::cout << x << " " << LBC(fluidState,paramCache,phaseIdx) << std::endl;
            return Opm::decay<LhsEval>(Opm::LBCviscosity<Scalar, ThisType>::LBC(fluidState,paramCache,phaseIdx));
            //if(phaseIdx == oilPhaseIdx) {
            //    return 5e-4;
            //} else {
            //    return 1e-5;
            //}
        }

        //! \copydoc BaseFluidSystem::enthalpy
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval enthalpy(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
        {
            const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
            const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));
            throw std::runtime_error("We don't use the enthalpy for non-isothermal runs");
        }


        //! \copydoc BaseFluidSystem::fugacityCoefficient
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& paramCache,
                                           unsigned phaseIdx,
                                           unsigned compIdx)
        {
            assert(0 <= phaseIdx && phaseIdx < numPhases);
            assert(0 <= compIdx && compIdx < numComponents);

            if (phaseIdx == oilPhaseIdx) {
#if 1
#warning HACK We use henry's law
#warning Use value for octane
                if (compIdx == DecaneIdx)
                    return 26.6e3/Opm::decay<LhsEval>(fluidState.pressure(oilPhaseIdx));
                else
                return 40e3/Opm::decay<LhsEval>(fluidState.pressure(oilPhaseIdx));
#else
                if (compIdx == CO2Idx)
                    return 500e3/Opm::decay<LhsEval>(fluidState.pressure(oilPhaseIdx));
                else {
                    return PengRobinsonMixture::computeFugacityCoefficient(fluidState,
								     paramCache,
								     phaseIdx,
								     compIdx);
                }
#endif
            } else if (phaseIdx == gasPhaseIdx) {
                return 1.0;
            } else {
                throw std::invalid_argument("expects oil or gas phase!");
            }
        }

        //! \copydoc BaseFluidSystem::diffusionCoefficient
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval diffusionCoefficient(const FluidState& /*fluidState*/,
                                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                            unsigned /*phaseIdx*/,
                                            unsigned compIdx)
        {
            // The diffusionCoefficient. From  Cadogan, S. P., Mistry, B., Wong, Y. et al. (2016). Diffusion Coefficients of Carbon Dioxide in Eight Hydrocarbon Liquids at Temperatures between (298.15 and 423.15) K at Pressures up to 69 MPa. Journal of Chemical & Engineering Data 61 (11): 3922-3932. https://doi.org/10.1021/acs.jced.6b00691
            // n-decane 6*e-9
            // n-octane 7.44*e-9
            //if (compIdx == CO2Idx)
            return 6.0e-9;

            //return 0.8e-9;
        }

    /*!
     * \brief Returns the interaction coefficient for two components.
     *
     * The values are from Ivar
     */
    static Scalar interactionCoefficient(unsigned comp1Idx, unsigned comp2Idx)
    {
        unsigned i = std::min(comp1Idx, comp2Idx);
        unsigned j = std::max(comp1Idx, comp2Idx);
        if (i == DecaneIdx && j == CO2Idx)
            return 0.1032;

        return 0;
    }

};

}//namespace opm

#endif // DECANE_CO2_FLUIDSYSTEM_HH
