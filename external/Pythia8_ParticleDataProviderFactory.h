#ifndef _AcerDet_External_Pythia8_ParticleDataProviderFactory_
#define _AcerDet_External_Pythia8_ParticleDataProviderFactory_

#pragma once

#include "../src/core/IParticleDataProviderFactory.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace external {
		
		//! Pythia8 specialization of AcerDET core::IParticleDataProvider
		class Pythia8_ParticleDataProviderFactory 
			: public IParticleDataProviderFactory {
		public:
		
			/**
			 * Creates new ParticleDataProvider instance.
			 * \result new ParticleDataProvider instance.
			 */
			ParticleDataProvider create();
		};

	}
}

#endif
