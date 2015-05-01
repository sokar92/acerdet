#ifndef _AcerDet_Analyse_Mis_
#define _AcerDet_Analyse_Mis_

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

		//! Missing energy analyse algorithm class
		/*!
		 * An algorithm for collecting missing energy.
		 */
		class Mis {
		private:
			Real64_t PTMUMIN; /*!< see Configuration */
			Real64_t ETAMAX;  /*!< see Configuration */
			Real64_t ETCELL;  /*!< see Configuration */
			Real64_t CALOTH;  /*!< see Configuration */
			
			Int32_t KEYHID;   /*!< see Configuration */
			Bool_t KEYSME;    /*!< see Configuration */
			Int32_t KFINVS;   /*!< see Configuration */
			
			Int32_t IEVENT;   /*!< Number of events computed by Mis algorithm. */
			
			IHistogramManager *histoManager; /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;          /*!< Indicates whether histograms are already registered and might be used now. */
		
		public:
			
			/** Constructor.
			 * Creates new instance of Missing energy analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histMgr a global instance of IHistogramManager
			 */
			Mis( const Configuration& conf, IHistogramManager* histMgr );
			
			/**
			 * Destructor.
			 */
			~Mis();
			
			/**
			 * Print information about Missing energy algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record in order to find missing energy.
			/*!
			 * \param input record.
			 * \param output record.
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );
			
			/**
			 * Print Missing energy algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
