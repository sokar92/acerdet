#ifndef _AcerDet_Core_ClusterData_
#define _AcerDet_Core_ClusterData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		struct ClusterData {
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
			 * 2 -
			 * 3 -
			 * 4 - eT
			 */
			
			Int32_t iepth, hits, state;
			Real64_t eta, phi, eT;
			
			ClusterData();
			ClusterData(const ClusterData&);
			
			ClusterData& operator = (const ClusterData&);
			
			static void sortByE(vector<ClusterData>&);
		};
	}
}

#endif
