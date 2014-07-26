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
			Real64_t eta, phi, pT;
			
			CellData();
			CellData(const CellData&);
			
			CellData& operator = (const CellData&);
			
			static void sortBy_pT(vector<CellData>&);
		
		private:
			static bool comparator_pT(const CellData&, const CellData&);
		};
	}
}

#endif
