#ifndef _AcerDet_IO_InputReader_
#define _AcerDet_IO_InputReader_

// HepMC event record
#include <HepMC/GenEvent.h>
using namespace HepMC;

#include "InputRecord.h"

namespace AcerDet {
	namespace io {
		/*
		 * Static converter between HepMC record and AcerDET InputRecord.
		 */
		class InputReader {
		private:
			/*
			 * Convert pdg_id code from HepMC to ParticleType
			 */
			static ParticleType getParticleType(int hepmc_code);

			/*
			 */
			static ParticleState getParticleStatus(GenParticle*);

		public:
			static InputRecord computeEvent(const GenEvent&);
		};

	}
}

#endif
