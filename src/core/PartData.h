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
			 * 0 -
			 * 1 -
			 * 2 - 
			 * 3 - 
			 * 4 - 
			 */
			
			/*
			 * P
			 * 0 - 
			 * 1 - 
			 * 2 - 
			 * 3 -  
			 * 4 - pT
			 */
			
			Real64_t pT;
			
			PartData();
			
			static void sortBy_pT(vector<PartData>&);
			
		private:
			static bool comparator_pT(const PartData&, const PartData&);
		};
	}
}

#endif
