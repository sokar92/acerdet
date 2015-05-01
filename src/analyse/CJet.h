#ifndef _AcerDet_Analyse_C_Jets_
#define _AcerDet_Analyse_C_Jets_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace analyse {

		//! CJet analyse algorithm class
		/*!
		 * An algorithm recognising Jets as CJets.
		 */	
		class CJet {
		private:
			Real64_t ETJET;  /*!< see Configuration */
			Real64_t RCONE;  /*!< see Configuration */
			Real64_t PTCMIN; /*!< see Configuration */
			Real64_t ETCMAX; /*!< see Configuration */
			Real64_t RJC;    /*!< see Configuration */

			Int32_t KEYHID;  /*!< see Configuration */
			Bool_t KEYBCL;   /*!< see Configuration */

			Int32_t IEVENT;  /*!< Number of events computed by CJet algorithm. */

			IHistogramManager *histoManager; /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;          /*!< Indicates whether histograms are already registered and might be used now. */
			
		public:
			
			/**
			 * Constructor.
			 * Creates new instance of CJet analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param hist a global instance of IHistogramManager
			 */
			CJet( const Configuration& conf, IHistogramManager* hist );
			
			/**
			 * Destructor.
			 */
			~CJet();

			/**
			 * Print information about CJet algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and find c-jets.
			/*!
			 * \param input record.
			 * \param output record.
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );

			/**
			 * Print CJet algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
