#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
#include "../core/Dum.h"
#include "../core/CellData.h"
#include "../core/ClusterData.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace io {
		
		class OutputRecord {
		public:			
			vector<CellData> cells;
			vector<ClusterData> clusters;
			
			vector<Particle> muons;
			vector<Particle> isolatedMuons;
			
			vector<Particle> electrons;
			vector<Particle> isolatedElectrons;
			
			vector<Particle> photons;
			vector<Particle> isolatedPhotons;
			
			vector<Particle> jets;
			vector<Particle> Bjets;
			vector<Particle> Cjets;
		
			OutputRecord();
			
			void clear();
		};

	}
}

#endif
