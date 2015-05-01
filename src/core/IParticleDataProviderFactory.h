#ifndef _AcerDet_Core_IParticleDataProviderFactory_
#define _AcerDet_Core_IParticleDataProviderFactory_

#include "ParticleDataProvider.h"

#pragma once

namespace AcerDet {
	namespace core {
		
		//! Abstract interface describing contracts for any external implementation of particle data provider supported by AcerDET.
		/**
		 * Abstract interface describing contracts for any external implementation of particle data provider
		 * supported by AcerDET.
		 */
		class IParticleDataProviderFactory {
		public:
		
			/**
			 * Virtual destructor.
			 */
			virtual ~IParticleDataProviderFactory() {}
			
			/**
			 * Create new ParticleDataProvider.
			 * Has to be implemented in subclass.
			 */
			virtual ParticleDataProvider create() = 0;
		};

	}
}

#endif
