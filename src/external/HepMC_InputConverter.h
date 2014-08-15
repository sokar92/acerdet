#ifndef _AcerDet_External_HepMC_InputConverter_
#define _AcerDet_External_HepMC_InputConverter_

// HepMC event record
#include <HepMC/GenEvent.h>
using namespace HepMC;

#include "../io/InputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace external {
		/*
		 * Static converter between HepMC record and AcerDET InputRecord.
		 */
		class HepMC_InputConverter {
		private:
			/*
			 * Convert pdg_id code from HepMC to ParticleType
			 */
			static ParticleType getParticleType(int hepmc_code);

			/*
			 */
			static ParticleState getParticleStatus(GenParticle*);

		public:
			static InputRecord convert(const GenEvent&);
		};

	}
}

#endif
