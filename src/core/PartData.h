#ifndef _AcerDet_Core_PartData_
#define _AcerDet_Core_PartData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		class PartData {
		public:
			/*
			 * K
			 * 0 - num
			 * 1 - state
			 * 2 - particleID
			 * 3 - motherState
			 * 4 - 
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - eta
			 * 3 - phi
			 * 4 - pT
			 */
			
			Int32_t num, particleID, motherState, state;
			Real64_t eta, phi, pT;
			
			PartData();
			
			static void sortBy_pT(vector<PartData>&);
			
		private:
			static bool comparator_pT(const PartData&, const PartData&);
		};
	}
}

#endif
