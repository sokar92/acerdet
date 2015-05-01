#ifndef _AcerDet_Core_ClusterData_
#define _AcerDet_Core_ClusterData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		/**
		 * AcerDET internal representation of detector cell grouping unit - Cluster.
		 */
		class ClusterData {
		public:
			/*
			 * K
			 * 0 -
			 * 1 -
			 * 2 - cellID
			 * 3 - hits
			 * 4 - status
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - eta_rec
			 * 3 - phi_rec 
			 * 4 - pT
			 */
			
			Int32_t cellID;   /*!< Unique id describing cluster position in detector (inherited from Cell) */
			Int32_t hits;     /*!< Number of energy peeks accumulated in cluster. */
			Int32_t status;   /*!< Cluster status. */
			Real64_t eta;     /*!< Cluster eta coordinate in detector. */
			Real64_t phi;     /*!< Cluster phi coordinate in detector. */
			Real64_t eta_rec; /*!< Custom total accumulated eta. */
			Real64_t phi_rec; /*!< Custom total accumulated phi. */
			Real64_t pT;      /*!< Total accumulated energy in custer. */
			
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
