#ifndef _AcerDet_External_Pythia8_ParticleDataProviderFactory_
#define _AcerDet_External_Pythia8_ParticleDataProviderFactory_

#pragma once

#include "../core/IParticleDataProviderFactory.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace external {
		
		class Pythia8_ParticleDataProviderFactory 
			: public IParticleDataProviderFactory {
		public:
			ParticleDataProvider create();
		};

	}
}

#endif
