#ifndef _AcerDet_Core_Dum_
#define _AcerDet_Core_Dum_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		struct Dum {
			Int32_t K[5];
			Real64_t P[5];
			
			Dum();
			Dum(const Dum&);
			
			Dum& operator = (const Dum&);
		};
	}
}

#endif
