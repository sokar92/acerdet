#ifndef _AcerDet_Core_Vector4_
#define _AcerDet_Core_Vector4_

#pragma once

#include "Vector3.h"
#include <ostream>
using namespace std;

namespace AcerDet {
	namespace core {

		/**
		 * A FourVector.
		 * Composition of 3-Vector and e/t as it's fourth coordinate.
		 */
		template<typename T>
		class Vector4 {
		public:
			/*! A 3-vector representing first 3 coordinates of four-vector. */
			Vector3<T> p;
			
			/*! A fourth coordinate of four-vector. */
			T e;
			
			/**
			 * Default constructor.
			 * \return a zero four-vector [0,0,0,0].
			 */
			Vector4() : p(), e(0) {}
			
			/**
			 * A copy constructor.
			 * \param v a four-vector to copy.
			 * \return a new four-vector filled with the same values as a given four-vector.
			 */
			Vector4(const Vector4<T>& v) : p(v.p), e(v.e) {}
			
			/**
			 * \param v a 3-vector representing first 3 coordinates of new four-vector.
			 * \param e a fourth coordinate of new four-vector.
			 * \return new four-vector with first 3 coordinates copied from 3-vector and fourth equal to parameter e.
			 */
			Vector4(const Vector3<T>& v, T val) : p(v), e(val) {}
			
			/**
			 * Creates a new four-vector with a given coordinates.
			 * \param x first cooridinate of new four-vector.
			 * \param y second cooridinate of new four-vector.
			 * \param z third cooridinate of new four-vector.
			 * \param w fourth cooridinate of new four-vector.
			 * \return new four-vector [x,y,z,w].
			 */
			Vector4(T x, T y, T z, T w) : p(Vector3<T>(x,y,z)), e(w) {}
			
			/**
			 * Equality operator.
			 * \param v four-vector to compare with.
			 * \return true if four-vectors are equal (in respect to they coordinates), false otherwise.
			 */
			Bool_t operator == (const Vector4<T>& v);
			
			/**
			 * Inequality operator.
			 * \param v four-vector to compare with.
			 * \return false if four-vectors are equal (in respect to they coordinates), true otherwise.
			 */
			Bool_t operator != (const Vector4<T>& v);
			
			/**
			 * Negation operator.
			 * Converts four-vector [x,y,z,w] to four-vector [-x,-y,-z,-w].
			 * \return negated four-vector.
			 */
			Vector4<T>& operator - ();
			
			/**
			 * Adds a given four-vector by coordinates to current four-vector.
			 * \param v four-vector to be added to current one.
			 * \return four-vector being a sum of current four-vector and the given one.
			 */
			Vector4<T>& operator += (const Vector4<T>& v);
			
			/**
			 * Substracts a given four-vector by coordinates from current four-vector.
			 * \param v four-vector to be substracted from current one.
			 * \return four-vector being result of substraction four-vector v from current one.
			 */
			Vector4<T>& operator -= (const Vector4<T>& v);
			
			/**
			 * Multiplies current four-vector by coordinates by given scalar value c.
			 * \param c a skalar value.
			 * \return current four-vector multiplied by scalar value.
			 */
			Vector4<T>& operator *= (T c);
			
			/**
			 * Divides current four-vector by coordinates by given scalar value c.
			 * Causes error when dividing by zero.
			 * \param a skalar non-zero value.
			 * \return current four-vector divided by scalar value.
			 */
			Vector4<T>& operator /= (T c);

			/**
			 * Adds a given four-vector by coordinates to current four-vector.
			 * \param v four-vector to be added to current one.
			 * \return four-vector being a sum of current four-vector and the given one.
			 */
			Vector4<T> operator + (const Vector4<T>& v);
			
			/**
			 * Substracts a given four-vector by coordinates from current four-vector.
			 * \param v four-vector to be substracted from current one.
			 * \return four-vector being result of substraction four-vector v from current one.
			 */
			Vector4<T> operator - (const Vector4<T>& v);
			
			/**
			 * Multiplies current four-vector by coordinates by given scalar value c.
			 * \param a skalar value.
			 * \return current four-vector multiplied by scalar value.
			 */
			Vector4<T> operator * (T c);
			
			/**
			 * Divides current four-vector by coordinates by given scalar value c.
			 * Causes error when dividing by zero.
			 * \param a skalar non-zero value.
			 * \return current four-vector divided by scalar value.
			 */
			Vector4<T> operator / (T c);
			
			/**
			 * Printing to standard output operator.
			 * \param str output stream - destination for printing process.
			 * \param v a four-vector to print out.
			 */
			friend ostream& operator << (ostream& str, const Vector4<T>& v) {
				str << "(" << v.p.x << ", " << v.p.y << ", " << v.p.z << ", " << v.e << ")";
				return str;
			}
		};
		
		/* Typedefs */
		typedef Vector4<Real64_t> Vector4f;
		
		/* Mutable operators */
		template<typename T>
		Vector4<T>& Vector4<T>::operator - () {
			p = -p; e = -e;
			return *this;
		}

		template<typename T>
		Vector4<T>& Vector4<T>::operator += (const Vector4<T>& v) {
			p += v.p;
			e += v.e;
			return *this;
		}

		template<typename T>
		Vector4<T>& Vector4<T>::operator -= (const Vector4<T>& v) {
			p -= v.p;
			e -= v.e;
			return *this;
		}

		template<typename T>
		Vector4<T>& Vector4<T>::operator *= (T coef) {
			p *= coef;
			e *= coef;
			return *this;
		}

		template<typename T>
		Vector4<T>& Vector4<T>::operator /= (T coef) {
			p /= coef;
			e /= coef;
			return *this;
		}

		/* Static operators */
		template<typename T>
		Vector4<T> Vector4<T>::operator + (const Vector4<T>& r) {
			return Vector4<T>(p + r.p, e + r.e);
		}

		template<typename T>
		Vector4<T> Vector4<T>::operator - (const Vector4<T>& r) {
			return Vector4<T>(p - r.p, e - r.e);
		}

		template<typename T>
		Vector4<T> Vector4<T>::operator * (T coef) {
			return Vector4<T>(p * coef, e * coef);
		}

		template<typename T>
		Vector4<T> Vector4<T>::operator / (T coef) {
			return Vector4<T>(p / coef, e / coef);
		}
		
		/* Comparison operators */
		template<typename T>
		Bool_t Vector4<T>::operator == (const Vector4<T>& v) {
			return p == v.p && e == v.e;
		}
		
		template<typename T>
		Bool_t Vector4<T>::operator != (const Vector4<T>& v) {
			return !(*this == v);
		}
		
	}
}

#endif
