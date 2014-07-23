#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
#include "../core/Dum.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace io {
		
		class OutputRecord {
		public:
			vector<Dum> vDum;
			
			vector<Dum> vCell;
			vector<Particle> vCluster;
			vector<Particle> vJet;
			vector<Particle> vElectron;
			vector<Particle> vPhoton;
			vector<Particle> vMuon;
		
			OutputRecord();
			
			void clear();
		};

	}
}

#endif
