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
			 * 4 - state
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - eta_rec
			 * 3 - phi_rec 
			 * 4 - eT
			 */
			
			Int32_t cellID; /*!< detailed description  */
			Int32_t hits; /*!< detailed description  */
			Int32_t state; /*!< detailed description  */
			Real64_t eta; /*!< detailed description  */
			Real64_t phi; /*!< detailed description  */
			Real64_t eta_rec; /*!< detailed description  */
			Real64_t phi_rec; /*!< detailed description  */
			Real64_t pT; /*!< detailed description  */
			
			ClusterData();
			
			static void sortBy_pT(vector<ClusterData>&);
			
		private:
			static bool comparator_pT(const ClusterData&, const ClusterData&);
		};
	}
}

#endif
