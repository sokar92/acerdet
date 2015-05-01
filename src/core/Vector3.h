#ifndef _AcerDet_Core_Vector3_
#define _AcerDet_Core_Vector3_

#pragma once

#include <cmath>
#include "Typedefs.h"

namespace AcerDet {
	namespace core {

		//! A 3-Vector.
		/**
		 * A 3-Vector.
		 * Composition of 3 independent coordinates.
		 */
		template<typename T>
		class Vector3 {
		public:
			T x; /*!< 3-vector first (X) coordinate. */
			T y; /*!< 3-vector second (Y) coordinate. */
			T z; /*!< 3-vector third (Z) coordinate. */

			/**
			 * Default constructor.
			 * \return a zero 3-vector [0,0,0].
			 */
			Vector3() : x(0), y(0), z(0) {};
			
			/**
			 * Creates a new 3-vector with a given coordinates.
			 * \param x first cooridinate of new 3-vector.
			 * \param y second cooridinate of new 3-vector.
			 * \param z third cooridinate of new 3-vector.
			 * \return new four-vector [x,y,z].
			 */
			Vector3(T X, T Y, T Z) : x(X), y(Y), z(Z) {}

			/**
			 * Negation operator.
			 * Converts 3-vector [x,y,z] to 3-vector [-x,-y,-z].
			 * \return negated 3-vector.
			 */
			Vector3<T>& operator - ();
			
			/**
			 * Adds a given 3-vector by coordinates to current 3-vector.
			 * \param v 3-vector to be added to current one.
			 * \return 3-vector being a sum of current 3-vector and the given one.
			 */
			Vector3<T>& operator += (const Vector3<T>& v);
			
			/**
			 * Substracts a given 3-vector by coordinates from current 3-vector.
			 * \param v 3-vector to be substracted from current one.
			 * \return 3-vector being result of substraction 3-vector v from current one.
			 */
			Vector3<T>& operator -= (const Vector3<T>& v);
			
			/**
			 * Multiplies current 3-vector by coordinates by given scalar value c.
			 * \param c a skalar value.
			 * \return current 3-vector multiplied by scalar value.
			 */
			Vector3<T>& operator *= (T c);
			
			/**
			 * Divides current 3-vector by coordinates by given scalar value c.
			 * Causes error when dividing by zero.
			 * \param a skalar non-zero value.
			 * \return current 3-vector divided by scalar value.
			 */
			Vector3<T>& operator /= (T c);

			/**
			 * Adds a given 3-vector by coordinates to current 3-vector.
			 * \param v 3-vector to be added to current one.
			 * \return 3-vector being a sum of current 3-vector and the given one.
			 */
			Vector3<T> operator + (const Vector3<T>& v);
			
			/**
			 * Substracts a given 3-vector by coordinates from current 3-vector.
			 * \param v 3-vector to be substracted from current one.
			 * \return 3-vector being result of substraction 3-vector v from current one.
			 */
			Vector3<T> operator - (const Vector3<T>& v);
			
			/**
			 * Multiplies current 3-vector by coordinates by given scalar value c.
			 * \param c a skalar value.
			 * \return current 3-vector multiplied by scalar value.
			 */
			Vector3<T> operator * (T c);
			
			/**
			 * Divides current 3-vector by coordinates by given scalar value c.
			 * Causes error when dividing by zero.
			 * \param a skalar non-zero value.
			 * \return current 3-vector divided by scalar value.
			 */
			Vector3<T> operator / (T c);
			
			/**
			 * Equality operator.
			 * \param v 3-vector to compare with.
			 * \return true if 3-vectors are equal (in respect to they coordinates), false otherwise.
			 */
			Bool_t operator == (const Vector3<T>& v);
			
			/**
			 * Inequality operator.
			 * \param v 3-vector to compare with.
			 * \return false if 3-vectors are equal (in respect to they coordinates), true otherwise.
			 */
			Bool_t operator != (const Vector3<T>& v); 

			/**
			 * Computes length (in Euklidian sense) of current 3-vector.
			 * \return Length of current 3-vector. ( sqrt(x*x + y*y + z*z) )
			 */
			inline T length() const;
			
			/**
			 * Normalizes current 3-vector.
			 * After normalization length of current vector is 1.
			 * \throw error when zero 3-vector.
			 */
			Vector3<T>& normalize ();
			
			/**
			 * Normalizes current 3-vector.
			 * Length of returned 3-vector is 1.
			 * \throw error when zero 3-vector.
			 * \return normalized form of current 3-vector.
			 */
			Vector3<T> getNormalized () const;

			/**
			 * Computes dot product between current 3-vector and given one.
			 * \param v second argument of dot product (first is current 3-vector).
			 * \return dot product as a scalar value.
			 */
			inline T dotProduct (const Vector3<T>& v) const;
			
			/**
			 * Computes cross product between current 3-vector and given one.
			 * \param v second argument of cross product (first is current 3-vector).
			 * \return cross product as new 3-vector perpendicular to current and given one.
			 */
			inline Vector3<T> crossProduct (const Vector3<T>& v) const;

			/**
			 * Computes angle (in radians) between current vector and given one.
			 * \param v 3-vector as second arm of the angle.
			 * \return angle in radians between current 3-vector and v.
			 */
			inline T angle (const Vector3<T>& v) const;
			
			/**
			 * Computes angle (in radians) between current 3-vector and OX axis.
			 * \return angle in radians between current 3-vector and OX axis.
			 */
			inline T angleOX () const;
			
			/**
			 * Computes angle (in radians) between current 3-vector and OY axis.
			 * \return angle in radians between current 3-vector and OY axis.
			 */
			inline T angleOY () const;
			
			/**
			 * Computes angle (in radians) between current 3-vector and OZ axis.
			 * \return angle in radians between current 3-vector and OZ axis.
			 */
			inline T angleOZ () const;
		};

		/* Typedefs */
		typedef Vector3<Real64_t> Vector3f;

		/* Mutable operators */
		template<typename T>
		Vector3<T>& Vector3<T>::operator - () {
			x = -x; y = -y; z= -z;
			return *this;
		}

		template<typename T>
		Vector3<T>& Vector3<T>::operator += (const Vector3<T>& v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}

		template<typename T>
		Vector3<T>& Vector3<T>::operator -= (const Vector3<T>& v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}

		template<typename T>
		Vector3<T>& Vector3<T>::operator *= (T coef) {
			x *= coef;
			y *= coef;
			z *= coef;
			return *this;
		}

		template<typename T>
		Vector3<T>& Vector3<T>::operator /= (T coef) {
			x /= coef;
			y /= coef;
			z /= coef;
			return *this;
		}

		/* Static operators */
		template<typename T>
		Vector3<T> Vector3<T>::operator + (const Vector3<T>& r) {
			return Vector3<T>(x + r.x, y + r.y, z + r.z);
		}

		template<typename T>
		Vector3<T> Vector3<T>::operator - (const Vector3<T>& r) {
			return Vector3<T>(x - r.x, y - r.y, z - r.z);
		}

		template<typename T>
		Vector3<T> Vector3<T>::operator * (T coef) {
			return Vector3<T>(x * coef, y * coef, z * coef);
		}

		template<typename T>
		Vector3<T> Vector3<T>::operator / (T coef) {
			return Vector3<T>(x / coef, y / coef, z / coef);
		}
		
		/* Comparison operators */
		template<typename T>
		Bool_t Vector3<T>::operator == (const Vector3<T>& v) {
			return x == v.x && y == v.y && z == v.z;
		}
		
		template<typename T>
		Bool_t Vector3<T>::operator != (const Vector3<T>& v) {
			return !(*this == v);
		}

		/* Normalization */
		template<typename T>
		inline T Vector3<T>::length() const {
			return sqrt(x*x + y*y + z*z);
		}

		template<typename T>
		Vector3<T>& Vector3<T>::normalize() {
			T len = length();
			x /= len;
			y /= len;
			z /= len;
			return *this;
		}

		template<typename T>
		Vector3<T> Vector3<T>::getNormalized() const {
			T len = length();
			return Vector3<T>(x / len, y / len, z / len);
		}

		/* Vector multiplication */
		template<typename T>
		inline T Vector3<T>::dotProduct(const Vector3<T>& v) const {
			return x*v.x + y*v.y + z*v.z;
		}

		template<typename T>
		inline Vector3<T> Vector3<T>::crossProduct(const Vector3<T>& v) const {
			return Vector3<T>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
		}

		/* Angles */
		template<typename T>
		inline T Vector3<T>::angle(const Vector3<T>& axis) const {
			T dot = dotProduct(axis);
			T axisLen = axis.length();
			T thisLen = length();
			T cosine = dot / (axisLen * thisLen);
			return acos (cosine);
		}

		template<typename T>
		inline T Vector3<T>::angleOX () const {
			return acos (x / length());
		}

		template<typename T>
		inline T Vector3<T>::angleOY () const {
			return acos (y / length());
		}

		template<typename T>
		inline T Vector3<T>::angleOZ () const {
			return acos (z / length());
		}

	}
}

#endif
