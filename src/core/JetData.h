#ifndef _AcerDet_Core_JetData_
#define _AcerDet_Core_JetData_

#pragma once

#include "Typedefs.h"
#include "StlDefs.h"

namespace AcerDet {
	namespace core {
		
		//! Type of partons jet.
		/**
		 * Type of partons jet.
		 */
		enum JetType {
			UNKNOWN = 0,  /*!< Jet not classified as b-jet or c-jet yet. */
			C_JET = 4,    /*!< Jet classified as c-jet. */
			B_JET = 5     /*!< Jet classified as b-jet. */
		};
		
		//! AcerDET internal representation of Jet.
		/**
		 * AcerDET internal representation of Jet.
		 */
		class JetData {
		public:
			/*
			 * K
			 * 0 -
			 * 1 - type
			 * 2 - 
			 * 3 - 
			 * 4 - 
			 */
			
			/*
			 * P
			 * 0 - eta
			 * 1 - phi
			 * 2 - eta_rec
			 * 3 - phi_rec 
			 * 4 - eT
			 */
			
			JetType type;     /*!< A type of jet (if recognised) */
			
			Real64_t eta;     /*!< Jet eta angle. */
			Real64_t phi;     /*!< Jet phi angle. */
			Real64_t eta_rec; /*!< Custom total accumulated eta angle. */
			Real64_t phi_rec; /*!< Custom total accumulated phi angle. */
			Real64_t pT;      /*!< Total energy accumulated in jet. */
			
			/**
			 * Default constructior.
			 * \return new structure describing jet.
			 */
			JetData();
			
			/**
			 * Checks whethear current jet is a b-jet.
			 * \return true if current jet is a b-jet, false otherwise.
			 */
			Bool_t isBJet() const;
			
			/**
			 * Checks whethear current jet is a c-jet.
			 * \return true if current jet is a c-jet, false otherwise.
			 */
			Bool_t isCJet() const;
			
			/**
			 * Sorts a list of jets by accumulated energy.
			 * \param v list of jets to sort.
			 */
			static void sortBy_pT(vector<JetData>& v);
			
		private:
			static bool comparator_pT(const JetData&, const JetData&);
		};
	}
}

#endif
