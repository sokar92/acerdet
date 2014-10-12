#ifndef _AcerDet_Core_ClusterData_
#define _AcerDet_Core_ClusterData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
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
			
			Int32_t cellID; /*!< unique id describing cluster position in detector (inherited from Cell) */
			Int32_t hits; /*!< sum of energy peeks accumulated by cluster */
			Int32_t status; /*!< cluster status */
			Real64_t eta; /*!< eta */
			Real64_t phi; /*!< phi */
			Real64_t eta_rec; /*!< custom accumulated eta */
			Real64_t phi_rec; /*!< custom accumulated phi */
			Real64_t pT; /*!< accumulated energy */
			
			ClusterData();
			
			static void sortBy_pT(vector<ClusterData>&);
			
		private:
			static bool comparator_pT(const ClusterData&, const ClusterData&);
		};
	}
}

#endif
