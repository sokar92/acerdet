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
			 * 1 - status
			 * 2 - barcode
			 * 3 - motherStatus
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
			
			Int32_t num; /*!< particle unique id (sequential number in event) */
			Int32_t barcode; /*!< particle barcode */
			Int32_t motherStatus; /*!< status of mother particle (if exists) */
			Int32_t status; /*!< particle status */
			Real64_t eta; /*!< eta */
			Real64_t phi; /*!< phi */
			Real64_t pT; /*!< energy */
			
			PartData();
			
			static void sortBy_pT(vector<PartData>&);
			
		private:
			static bool comparator_pT(const PartData&, const PartData&);
		};
	}
}

#endif
