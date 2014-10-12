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
			
			enum JetType {
				UNKNOWN = 0,
				C_JET = 4,
				B_JET = 5
			} type; /*!< A type of jet (if recognised) */
			
			Real64_t eta; /*!< eta (inherited from ?) */
			Real64_t phi; /*!< phi (inherited from ?) */
			Real64_t eta_rec; /*!< custom accumulated eta (inherited from ?) */
			Real64_t phi_rec; /*!< custom accumulated phi (inherited from ?) */
			Real64_t pT; /*!< accumulated energy */
			
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
