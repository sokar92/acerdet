#ifndef _AcerDet_Core_Functions_
#define _AcerDet_Core_Functions_

#pragma once

#include "Typedefs.h"

namespace AcerDet {
	namespace core {
		
		/*
		 * Compute radius from given coordinates
		 */
		inline double radius(double, double);
		
		/*
		 * Returns angle from given coordinates
		 * in range [-PI,PI]
		 */
		double angle(double, double);
		
		Int32_t kfcomp(Real64_t);
		Int32_t kuchge(Real64_t);
	}
}

#endif
