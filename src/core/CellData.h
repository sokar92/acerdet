#ifndef _AcerDet_Core_CellData_
#define _AcerDet_Core_CellData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		//! AcerDET internal representation of detector accumulation basic unit - Cell.
		/**
		 * AcerDET internal representation of detector accumulation basic unit - Cell.
		 */
		class CellData {
		public:
			/*
			 * K
			 * 0 -
			 * 1 -
			 * 2 - cellID
			 * 3 - hits
			 * 4 - status
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 -
			 * 3 -
			 * 4 - pT
			 */
			Int32_t cellID;  /*!< Unique cell id describing it's position in detector. (is computable from eta and phi). */
			Int32_t status;  /*!< Cell status. */
			Int32_t hits;    /*!< Number of energy peeks contained in this cell. */
			Real64_t eta;    /*!< Cell eta coordinate in detector. */
			Real64_t phi;    /*!< Cell phi coordinate in detector. */
			Real64_t pT;     /*!< Total energy accumulated in Cell. */
			
			/**
			 * Default constructor.
			 * \return new structure describing a cell in detector.
			 */
			CellData();
			
			Real64_t pX() const;
			
			Real64_t pY() const;
			
			Real64_t pZ() const;
			
			Real64_t e() const;
			
			/**
			 * Sorst a list of cells by accumulated energy.
			 * \param v list of cells to sort.
			 */
			static void sortBy_pT(vector<CellData>& v);
		
		private:
			static bool comparator_pT(const CellData&, const CellData&);
		};
	}
}

#endif
