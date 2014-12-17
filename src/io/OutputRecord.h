#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
#include "../core/CellData.h"
#include "../core/ClusterData.h"
#include "../core/JetData.h"
#include "../core/PartData.h"
#include "../core/MisData.h"
using namespace AcerDet::core;

namespace AcerDet {
	namespace io {
		
		class OutputRecord {
		public:			
			vector<CellData> Cells; /*!< detailed description  */
			vector<ClusterData> Clusters; /*!< detailed description  */
			
			vector<PartData> Muons; /*!< detailed description  */
			vector<PartData> NonisolatedMuons; /*!< detailed description  */

			vector<PartData> Electrons; /*!< detailed description  */
			vector<PartData> Photons; /*!< detailed description  */
			
			vector<JetData> Jets; /*!< detailed description  */
			MisData Miss; /*!< detailed description  */
			
			OutputRecord();
			
			void clear();
		};

	}
}

#endif
