#ifndef _AcerDet_Core_Functions_
#define _AcerDet_Core_Functions_

#pragma once

#include "Typedefs.h"

namespace AcerDet {
	namespace core {
		
		/*
		 * Compute radius from given coordinates
		 */
		inline double radius(double x, double y);
		
		/*
		 * Returns angle from given coordinates
		 * in range [-PI,PI]
		 */
		double angle(double x, double y);
		
		Int32_t kfcomp(Real64_t kf);
		
	}
}

#endif
