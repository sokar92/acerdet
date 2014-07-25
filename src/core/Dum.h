#ifndef _AcerDet_Core_Dum_
#define _AcerDet_Core_Dum_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		struct Dum {
			/*
			 * 0 -
			 * 1 -
			 * 2 - IEPTH
			 * 3 - hits
			 * 4 - state
			 */
			Int32_t K[5];
			
			/*
			 * 0 - eta
			 * 1 - phi
			 * 2 -
			 * 3 -
			 * 4 - eT
			 */
			Real64_t P[5];
			
			Dum();
			Dum(const Dum&);
			
			Dum& operator = (const Dum&);
			
			static vector<Dum> sortByE(const vector<Dum>&);
		};
	}
}

#endif
