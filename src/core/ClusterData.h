#ifndef _AcerDet_Core_ClusterData_
#define _AcerDet_Core_ClusterData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		//! AcerDET internal representation of detector cell grouping unit - Cluster.
		/**
		 * AcerDET internal representation of detector cell grouping unit - Cluster.
		 */
		class ClusterData {
		public:
			Int32_t cellID;   /*!< Unique id describing cluster position in detector (inherited from Cell) */
			Int32_t hits;     /*!< Number of energy peeks accumulated in cluster. */
			Int32_t status;   /*!< Cluster status. */
			Real64_t eta;     /*!< Cluster eta coordinate in detector. */
			Real64_t phi;     /*!< Cluster phi coordinate in detector. */
			Real64_t eta_rec; /*!< Custom total accumulated eta. */
			Real64_t phi_rec; /*!< Custom total accumulated phi. */
			Real64_t pT;      /*!< Total accumulated energy in custer. */
			Bool_t alreadyUsed; /*!< If used in jet collection. */
			
			/**
			 * Default constructor.
			 * \return new structure describing detector cluster.
			 */
			ClusterData();
			
			/**
			 * Sorts clusters by total accumulated energy.
			 * \param v list of cluster to sort.
			 */
			static void sortBy_pT(vector<ClusterData>& v);
			
		private:
			static bool comparator_pT(const ClusterData&, const ClusterData&);
		};
	}
}

#endif
