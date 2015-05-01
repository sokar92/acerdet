#ifndef _AcerDet_Core_Smearing_
#define _AcerDet_Core_Smearing_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"
#include <random>

namespace AcerDet {
	namespace core {
		
		/**
		 * Utility class with methods for energy smearing computation.
		 */
		class Smearing {
		private:
			static std::default_random_engine generator;             /*!< C++ Standard library random variable generator. */
			static std::normal_distribution<Real64_t> distribution;  /*!< Variable from normal distribution over doubles. */
		public:
			
			/**
			 * Computes energy smearing for electrons.
			 * \param ene electron energy.
			 * \return electron energy after smearing.
			 */
			static Real64_t forElectron(Real64_t ene);
			
			 /**
			 * Computes energy smearing for hadrons.
			 * \param ene hadron energy.
			 * \param eta hadronic particle angle.
			 * \param caloth calorimetric energy deposition factor.
			 * \return hadron energy after smearing.
			 */
			static Real64_t forHadron(Real64_t ene, Real64_t eta, Real64_t caloth);
			
			/**
			 * Computes energy smearing for muons.
			 * \param pt muon energy.
			 * \return muon energy after smearing.
			 */
			static Real64_t forMuon(Real64_t pt);
			
			/**
			 * Computes energy smearing for photons.
			 * \param ene photon energy.
			 * \return photon energy after smearing.
			 */
			static Real64_t forPhoton(Real64_t ene);
		};
	}
}

#endif
