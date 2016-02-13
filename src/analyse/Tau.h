#ifndef _AcerDet_Analyse_Tau_
#define _AcerDet_Analyse_Tau_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

#include "./IAnalysePhase.h"

namespace AcerDet {
	namespace analyse {

		//! Tau analyse algorithm class
		/*!
		 * An algorithm for Tau reconstruction.
		 */
		class Tau : public IAnalysePhase {
		private:
			Real64_t ETJET;  /*!< see Configuration */
			Real64_t PTTAU;  /*!< see Configuration */
			Real64_t ETATAU; /*!< see Configuration */
			Real64_t RJTAU;  /*!< see Configuration */
			Real64_t PTFRAC; /*!< see Configuration */

			Int32_t KEYHID;  /*!< see Configuration */
			Bool_t KEYTAU;   /*!< see Configuration */

			Int32_t IEVENT;  /*!< Number of events computed by Tau algorithm. */

			IHistogramManager *histoManager;  /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;           /*!< Indicates whether histograms are already registered and might be used now. */
		public:
			
			/** Constructor.
			 * Creates new instance of Tau analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histMgr a global instance of IHistogramManager
			 */
			Tau( const Configuration& conf, IHistogramManager* histMgr );
			
			/**
			 * Destructor.
			 */
			~Tau();

			/**
			 * Print information about Tau algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and find taus.
			/*!
			 * \param input record.
			 * \param output record.
			 * \param weight event weight.
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output,
				Real64_t weight );

			/**
			 * Print Tau algorithm execution results to standard output.
			 */
			void printResults() const;
		};
		
	}
}

#endif
