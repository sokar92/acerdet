#ifndef _AcerDet_Core_MisData_
#define _AcerDet_Core_MisData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		class MisData {
		public:
			Real64_t PXREC, PYREC;
			Real64_t PXSUM, PYSUM;
			Real64_t PXXCALO, PYYCALO;
			Real64_t SUMET;
			
			MisData() :
				PXREC (0.0), PYREC (0.0), PXSUM (0.0), PYSUM(0.0),
				PXXCALO (0.0), PYYCALO (0.0), SUMET(0.0);
		};
	}
}

#endif
