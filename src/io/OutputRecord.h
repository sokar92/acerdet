#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace io {

		class OutputRecord {
		public:
			vector<Particle> vDum; // used in cell & cluster
			vector<Particle> vCell;
		
			vector<Particle> vJets;
		
			OutputRecord();
		};

	}
}

#endif
