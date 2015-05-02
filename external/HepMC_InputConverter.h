#ifndef _AcerDet_External_HepMC_InputConverter_
#define _AcerDet_External_HepMC_InputConverter_

// HepMC event record
#include <HepMC/GenEvent.h>
using namespace HepMC;

#include "../src/io/InputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace external {
		
		//! Converter from HepMC event record to AcerDET input record.
		/**
		 * Static converter between HepMC record and AcerDET InputRecord.
		 */
		class HepMC_InputConverter {
		private:
		
			/**
			 * Convert pdg_id code from HepMC to ParticleType.
			 * \param hepmc_code HepMC internal particle code.
			 * \return Enumerated AcerDET type of particle.
			 */
			static ParticleType getParticleType(int hepmc_code);

			/**
			 * Convert HepMC part to AcerDET particle and reads it's status.
			 * \param part HepMC particle description.
			 * \return Enumerated AcerDET particle status.
			 */
			static ParticleStatus getParticleStatus(GenParticle* part);

		public:
		
			/**
			 * Converts HepMC event to AcerDET input record.
			 * \param event HepMC event.
			 * \result AcerDET input record.
			 */
			static InputRecord convert(const GenEvent& event);
		};

	}
}

#endif
