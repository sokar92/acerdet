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
		
		//! Data record describing AcerDET analyse output.
		/**
		 * Data record containing all reconstructed partons by AcerDET analyse algorithms.
		 */
		class OutputRecord {
		public:			
			vector<CellData> Cells;            /*!< List of created cells. */
			vector<ClusterData> Clusters;      /*!< List of created clusters. */
			
			vector<PartData> Muons;            /*!< List of reconstructed muons. */
			vector<PartData> NonisolatedMuons; /*!< List of reconstructed and non-isolated muons. */

			vector<PartData> Electrons;        /*!< List of reconstructed electrons. */
			vector<PartData> Photons;          /*!< List of reconstructed photons. */
			
			vector<JetData> Jets;              /*!< List of reconstructed Jets (including b-jets and c-jets). */
			MisData Miss;                      /*!< Missing energy statistics. */
			
			/**
			 * Default constructor.
			 * \return new AcerDET output data record instance.
			 */
			OutputRecord();
			
			/**
			 * Clear whole record.
			 * Deletes all data stored in current output record.
			 */
			void clear();
		};

	}
}

#endif
