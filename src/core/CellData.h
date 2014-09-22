#ifndef _AcerDet_Core_CellData_
#define _AcerDet_Core_CellData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		class CellData {
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
			 * 2 -
			 * 3 -
			 * 4 - eT
			 */
			
			Int32_t cellID; /*!< detailed description  */
			Int32_t hits; /*!< detailed description  */
			Int32_t state; /*!< detailed description  */
			Real64_t eta; /*!< detailed description  */
			Real64_t phi; /*!< detailed description  */
			Real64_t pT; /*!< detailed description  */
			
			CellData();
			
			static void sortBy_pT(vector<CellData>&);
		
		private:
			static bool comparator_pT(const CellData&, const CellData&);
		};
	}
}

#endif
