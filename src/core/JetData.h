#ifndef _AcerDet_Core_JetData_
#define _AcerDet_Core_JetData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		class JetData {
		public:
			/*
			 * K
			 * 0 -
			 * 1 - type
			 * 2 - 
			 * 3 - 
			 * 4 - 
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - eta_rec
			 * 3 - phi_rec 
			 * 4 - eT
			 */
			
			Int32_t type; /*!< detailed description  */
			Real64_t eta; /*!< detailed description  */
			Real64_t phi; /*!< detailed description  */
			Real64_t eta_rec; /*!< detailed description  */
			Real64_t phi_rec; /*!< detailed description  */
			Real64_t pT; /*!< detailed description  */
			
			JetData();
			
			Bool_t isBJet() const;
			Bool_t isCJet() const;
			
			static void sortBy_pT(vector<JetData>&);
			
		private:
			static bool comparator_pT(const JetData&, const JetData&);
		};
	}
}

#endif
