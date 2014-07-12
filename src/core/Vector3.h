#ifndef _AcerDet_Core_Vector3_
#define _AcerDet_Core_Vector3_

#pragma once

#include <cmath>

namespace AcerDet {
	namespace core {

		template<typename T>
		class Vector3 {
		public:
			T x,y,z;

			Vector3() : x(0), y(0), z(0) {};
			Vector3(T X, T Y, T Z) : x(X), y(Y), z(Z) {}

			Vector3<T>& operator - ();
			Vector3<T>& operator += (const Vector3<T>&);
			Vector3<T>& operator -= (const Vector3<T>&);
			Vector3<T>& operator *= (T);
			Vector3<T>& operator /= (T);

			Vector3<T> operator + (const Vector3<T>&, const Vector3<T>&);
			Vector3<T> operator - (const Vector3<T>&, const Vector3<T>&);
			Vector3<T> operator * (const Vector3<T>&, T);
			Vector3<T> operator / (const Vector3<T>&, T);

			inline T length() const;
			
			Vector3<T>& normalize ();
			Vector3<T> getNormalized () const;

			inline T dotProduct (const Vector3<T>&) const;
			inline Vector3<T> crossProduct (const Vector3<T>&) const;

			inline T angle (const Vector3<T>&) const;
			inline T angleOX () const;
			inline T angleOY () const;
			inline T angleOZ () const;
		};

		/* Typedefs */
		typedef Vector3<double> Vector3f;

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
		Vector3<T> vector3<T>::operator + (const Vector3<T>& l, const Vector3<T>& r) {
			return Vector3<T>(l.x + r.x, l.y + r.y, l.z + r.z);
		}

		template<typename T>
		Vector3<T> vector3<T>::operator - (const Vector3<T>& l, const Vector3<T>& r) {
			return Vector3<T>(l.x - r.x, l.y - r.y, l.z - r.z);
		}

		template<typename T>
		Vector3<T> vector3<T>::operator * (const Vector3<T>& v, T coef) {
			return Vector3<T>(v.x * coef, v.y * coef, v.z * coef);
		}

		template<typename T>
		Vector3<T> vector3<T>::operator / (const Vector3<T>& v, T coef) {
			return Vector3<T>(v.x / coef, v.y / coef, v.z / coef);
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
