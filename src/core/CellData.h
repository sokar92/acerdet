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
			 * 4 - status
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 -
			 * 3 -
			 * 4 - pT
			 */
			Int32_t cellID; /*!< unique cell id describing it's position in detector */
			Int32_t status; /*!< cell status */
			Int32_t hits; /*!< number of energy peeks contained by this cell */
			Real64_t eta; /*!< eta */
			Real64_t phi; /*!< phi */
			Real64_t pT; /*!< energy accumulated by cell */
			
			CellData();
			
			static void sortBy_pT(vector<CellData>&);
		
		private:
			static bool comparator_pT(const CellData&, const CellData&);
		};
	}
}

#endif
