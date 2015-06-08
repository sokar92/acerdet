#ifndef _AcerDet_Core_ObjectData_
#define _AcerDet_Core_ObjectData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		//! Internal AcerDET representation of object.
		/**
		 * Internal AcerDET representation of object.
		 */
		class ObjectData {
		public:
			/*
			 * K
			 * 0 - num
			 * 1 - status
			 * 2 - barcode
			 * 3 - motherStatus
			 * 4 - 
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - eta
			 * 3 - phi
			 * 4 - pT
			 */
			
			Int32_t num;          /*!< Object unique number in event (sequential). */
			Int32_t barcode;      /*!< Object barcode. */
			Int32_t motherStatus; /*!< Object mother's status (if mother exists). */
			Int32_t status;       /*!< Object status (not converted). */
			Real64_t eta;         /*!< Object eta angle. */
			Real64_t phi;         /*!< Object phi angle. */
			Real64_t pT;          /*!< Object energy. */
			
			/**
			 * Default constructor.
			 * \return new structure describing object in AcerDET.
			 */
			ObjectData();
			
			/**
			 * Sorts a list of objects by energy.
			 * \param vec a vector of objects to be sorted.
			 */
			static void sortBy_pT(vector<ObjectData>& vec);
			
		private:
		
			/**
			 * Object comparator.
			 * \param l first object to compare.
			 * \param r second object to compare.
			 */
			static bool comparator_pT(const ObjectData& l, const ObjectData& r);
		};
	}
}

#endif
