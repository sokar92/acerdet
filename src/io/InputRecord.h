#ifndef _AcerDet_IO_InputRecord_
#define _AcerDet_IO_InputRecord_

#include "../core/Particle.h"
using namespace AcerDet::core;

#include <vector>
using namespace std;

namespace AcerDet {
	namespace io {

		class InputRecord {
		public:
			vector<Particle> parts; /*!< A list of incomming Particles */

		public:
			/*
			 * Constructor
			 * creates an input record from given particle list
			 */
			InputRecord(const vector<Particle>&);
			
			/*
			 * Access the particle list
			 */
			const vector<Particle>& particles() const;

		};

	}
}

#endif
