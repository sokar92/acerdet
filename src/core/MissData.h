#ifndef _AcerDet_Core_MissData_
#define _AcerDet_Core_MissData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		//! AcerDET internal representation of missing energy.
		/**
		 * AcerDET internal representation of missing energy.
		 */
		class MissData {
		public:
			Real64_t PXREC, PYREC;
			Real64_t PXSUM, PYSUM;
			Real64_t PXXCALO, PYYCALO;
			Real64_t SUMET;
			
			/**
			 * Default constructor.
			 * \return new structure describing missing energy.
			 */
			MissData();
			
			/**
			 * Clears all stored data about missing energy.
			 */
			void clear();
		};
	}
}

#endif
