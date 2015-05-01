#ifndef _AcerDet_Core_Functions_
#define _AcerDet_Core_Functions_

#pragma once

#include "Typedefs.h"
#include "Consts.h"
#include "Particle.h"

#include <vector>
#include <algorithm>
using namespace std;

namespace AcerDet {
	namespace core {
		
		/**
		 * Saturates given angle to range [-pi,pi].
		 * \param angle angle in radians
		 * \return equivalent angle to given but in range [-pi, pi].
		 */
		template<typename T>
		inline T saturatePi(T angle) {
			double res = angle;
			if (angle > PI)
				res = angle - 2.0 * PI;
			if (angle < -PI)
				res = angle + 2.0 * PI;
			return res;
		}
		
		/**
		 * Compute radius from given coordinates in XY plane.
		 * \param x coordinate
		 * \param y coordinate
		 * \return vector radius (distance from (0,0) to (x,y)).
		 */
		template<typename T>
		inline T radius(T x, T y) {
			return (T)sqrt(x*x + y*y);
		}
		
		/**
		 * Signum function
		 * \param val value
		 * \return -1 if value is less than zero, +1 if value us greater than zero, 0 otherwise
		 */
		template<typename T>
		T sgn(T val) {
			if (val > 0) return (T)(+1);
			if (val < 0) return (T)(-1);
			return (T)(0);
		}
		
		/**
		 * Return first value with the sign of second
		 * \return |arg1| * sgn(arg2)
		 */
		template<typename T, typename Q>
		T sign(T val, Q coef) {
			Q sg = sgn(coef);
			return (T)(abs(val) * sg);
		}
		
		/**
		 * Computes angle from given coordinates in range [-PI,PI]
		 * \param x coordinate
		 * \param y coordinate
		 * \return (x,y) vector angle in XY plane saturated to range [-pi,pi].
		 */
		template<typename T>
		T angle(T x, T y) {
			Real64_t epsi = 1E-20;
			T r = radius(x,y);
			
			if (r < epsi)
				return (T)(0);
				
			if (abs(x)/r < 0.8) {
				return sign(acos(x/r),y);
			} else {
				T angle = asin(y/r);
				if (x < 0 && angle >= 0)
					return PI - angle;
				else if (x < 0)
					return -PI - angle;
				else
					return angle;
			}
			
			return radius(x,y);
		}
		
		/**
		 * Check if given particle is a hard process praticle.
		 * \param parts list of particles from event.
		 * \param i index of particle to check.
		 * \returns true if particle is final and has mother from set {23,24,25}.
		 */
		Bool_t isHardProcess(const vector<Particle>& parts, int i);
	}
}

#endif
