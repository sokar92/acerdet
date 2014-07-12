#ifndef _AcerDet_IO_InputReader_
#define _AcerDet_IO_InputReader_

// HepMC event record
#include <HepMC/GenEvent.h>
using namespace HepMC;

#include "InputRecord.h"

namespace AcerDet {
	namespace io {
		
		/*
		 * Set of codes describing particles
		 * Montecarlo particle numbering scheme
		 * from page http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
		 */
		enum PDGcode {
			QUARK_D = 1,
			QUARK_U = 2,
			QUARK_S = 3,
			QUARK_C = 4,
			QUARK_B = 5,
			QUARK_T = 6,
			QUARK_Bp = 7,
			QUARK_Tp = 8,

			LEPT_ELECTRON = 11,
			LEPT_MUON = 13,
			LEPT_TAU = 15
		};

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
			static ParticleState getParticleStatus(int hepmc_status);

		public:
			static InputRecord computeEvent(const GenEvent&);
		};

	}
}

#endif
