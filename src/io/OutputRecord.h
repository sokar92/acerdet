#ifndef _AcerDet_IO_OutputRecord_
#define _AcerDet_IO_OutputRecord_

#include "../core/Particle.h"
#include "../core/CellData.h"
#include "../core/ClusterData.h"
#include "../core/JetData.h"
#include "../core/ObjectData.h"
#include "../core/MissData.h"
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
			
			vector<ObjectData> Muons;            /*!< List of reconstructed muons. */
			vector<ObjectData> NonisolatedMuons; /*!< List of reconstructed and non-isolated muons. */

			vector<ObjectData> Electrons;        /*!< List of reconstructed electrons. */
			vector<ObjectData> Photons;          /*!< List of reconstructed photons. */
			
			vector<JetData> Jets;              /*!< List of reconstructed Jets (including b-jets and c-jets). */
			MissData Miss;                      /*!< Missing energy statistics. */
			
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
