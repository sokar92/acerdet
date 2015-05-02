#ifndef _AcerDet_IO_InputRecord_
#define _AcerDet_IO_InputRecord_

#include "../core/Particle.h"
using namespace AcerDet::core;

#include <vector>
using namespace std;

namespace AcerDet {
	namespace io {

		//! Data record describing AcerDET analyse input.
		/**
		 * Data record containing all particles from MonteCarlo event generator.
		 */
		class InputRecord {
		public:
			vector<Particle> parts; /*!< A list of incomming Particles. */

		public:
			
			/**
			 * Default constructor.
			 * \param vec list of particles from MonteCarlo event generator.
			 * \return new AcerDET input data record instance.
			 */
			InputRecord(const vector<Particle>& vec);
			
			/**
			 * Get list of all particles stored in current input record.
			 */
			const vector<Particle>& particles() const;

		};

	}
}

#endif
