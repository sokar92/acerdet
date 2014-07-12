#ifndef _AcerDet_Core_ParticleFactory_
#define _AcerDet_Core_ParticleFactory_

#include "Particle.h"

namespace AcerDet {
	namespace core {
		
		class ParticleFactory {
		public:
			/* Ctors */
			ParticleFactory() {}
			
			static Particle createBeam( Real64_t, Real64_t, Real64_t, Real64_t );
		};

	}
}

#endif
