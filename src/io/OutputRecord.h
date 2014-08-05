#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
#include "../core/CellData.h"
#include "../core/ClusterData.h"
#include "../core/PartData.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace io {
		
		class OutputRecord {
		public:			
			vector<CellData> Cells;
			vector<ClusterData> Clusters;
			
			vector<PartData> Muons;
			vector<PartData> Electrons;
			vector<PartData> Photons;
			
			vector<Particle> Jets;
			vector<Particle> BJets;
			vector<Particle> CJets;
		
			OutputRecord();
			
			void clear();
		};

	}
}

#endif
