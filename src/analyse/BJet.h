#ifndef _AcerDet_Analyse_B_Jets_
#define _AcerDet_Analyse_B_Jets_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		//! BJet analyse algorithm class
		/*!
		 * An algorithm recognising Jets as BJets.
		 */	
		class BJet {
		private:
			Real64_t ETJET;  /*!< see Configuration */
			Real64_t RCONE;  /*!< see Configuration */
			Real64_t PTBMIN; /*!< see Configuration */
			Real64_t ETBMAX; /*!< see Configuration */
			Real64_t RJB;    /*!< see Configuration */

			Int32_t KEYHID;  /*!< see Configuration */
			Bool_t KEYBCL;   /*!< see Configuration */

			Int32_t IEVENT;  /*!< see Configuration */
			
			IHistogramManager *histoManager; /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;          /*!< Indicates whether histograms are already registered and might be used now. */

		public:
		
			/**
			 * Constructor.
			 * Creates new instance of BJet analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param hist a global instance of IHistogramManager
			 */
			BJet( const Configuration&, IHistogramManager* );
			
			/**
			 * Destructor.
			 */
			~BJet();
			
			/**
			 * Print information about BJet algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and find b-jets.
			/*!
			 * \param input record.
			 * \param output record.
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );
			
			/**
			 * Print BJet algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
