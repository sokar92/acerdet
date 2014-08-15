#ifndef _AcerDet_Core_IParticleDataProviderFactory_
#define _AcerDet_Core_IParticleDataProviderFactory_

#include "ParticleDataProvider.h"

#pragma once

namespace AcerDet {
	namespace core {
		
		/**
		 * Interface for external ParticleDataProvider factory.
		 */
		class IParticleDataProviderFactory {
		public:
			/**
			 * Create new ParticleDataProvider.
			 * Has to be implemented in subclass.
			 */
			virtual ParticleDataProvider create() = 0;
		};

	}
}

#endif
