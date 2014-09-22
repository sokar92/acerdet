#ifndef _AcerDet_Core_Smearing_
#define _AcerDet_Core_Smearing_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"
#include <random>

namespace AcerDet {
	namespace core {
		class Smearing {
		private:
			static std::default_random_engine generator; /*!< detailed description  */
			static std::normal_distribution<Real64_t> distribution; /*!< detailed description  */
		public:
			/*
			 * parametrizes electron resolution 
			 */
			static Real64_t forElectron(Real64_t ene);
			
			/*
			 * parametrizes smearing for calorymetric energy deposition
			 */
			static Real64_t forHadron(Real64_t ene, Real64_t eta, Real64_t caloth);
			
			/*
			 * parametrizes muon resolution 
			 */
			static Real64_t forMuon(Real64_t pt);
			
			/*
			 * parametrizes photon resolution 
			 */
			static Real64_t forPhoton(Real64_t ene);
		};
	}
}

#endif
