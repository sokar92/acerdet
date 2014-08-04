#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
#include "../core/CellData.h"
#include "../core/ClusterData.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace io {
		
		class OutputRecord {
		public:			
			vector<CellData> Cells;
			vector<ClusterData> Clusters;
			
			vector<Particle> Muons;
			vector<Particle> IsolatedMuons;
			
			vector<Particle> Electrons;
			vector<Particle> IsolatedElectrons;
			
			vector<Particle> Photons;
			vector<Particle> IsolatedPhotons;
			
			vector<Particle> Jets;
			vector<Particle> BJets;
			vector<Particle> CJets;
		
			OutputRecord();
			
			void clear();
		};

	}
}

#endif
