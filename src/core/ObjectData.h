#ifndef _AcerDet_Core_ObjectData_
#define _AcerDet_Core_ObjectData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"
#include "Vector4.h"

namespace AcerDet {
	namespace core {
		
		//! Internal AcerDET representation of object.
		/**
		 * Internal AcerDET representation of object.
		 */
		class ObjectData {
		public:
			Int32_t num;          /*!< Object unique number in event (sequential). */
			Int32_t pdg_id;       /*!< Object pdg_id. */
			Vector4f p;           /*!< Object 4-vector. */
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
