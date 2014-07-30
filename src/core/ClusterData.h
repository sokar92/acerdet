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
			 * 2 - IEPTH
			 * 3 - hits
			 * 4 - state
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - ? TODO : eta_rec
			 * 3 - ? TODO : phi_rec 
			 * 4 - eT
			 */
			
			Int32_t iepth, hits, state;
			Real64_t eta, phi, pT;
			
			ClusterData();
			ClusterData(const ClusterData&);
			
			ClusterData& operator = (const ClusterData&);
			
			static void sortBy_pT(vector<ClusterData>&);
			
		private:
			static bool comparator_pT(const ClusterData&, const ClusterData&);
		};
	}
}

#endif
