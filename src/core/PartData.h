#ifndef _AcerDet_Core_PartData_
#define _AcerDet_Core_PartData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		/**
		 * Internal AcerDET representation of parton.
		 */
		class PartData {
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
			
			Int32_t num;          /*!< Parton unique number in event (sequential). */
			Int32_t barcode;      /*!< Parton barcode. */
			Int32_t motherStatus; /*!< Parton mother's status (if mother exists). */
			Int32_t status;       /*!< Parton status (not converted). */
			Real64_t eta;         /*!< Parton eta angle. */
			Real64_t phi;         /*!< Parton phi angle. */
			Real64_t pT;          /*!< Parton energy. */
			
			/**
			 * Default constructor.
			 * \return new structure describing parton in AcerDET.
			 */
			PartData();
			
			/**
			 * Sorts a list of partons by energy.
			 * \param vec a vector of partons to be sorted.
			 */
			static void sortBy_pT(vector<PartData>& vec);
			
		private:
		
			/**
			 * Parton comparator.
			 * \param l first parton to compare.
			 * \param r second parton to compare.
			 */
			static bool comparator_pT(const PartData& l, const PartData& r);
		};
	}
}

#endif
