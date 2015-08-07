#ifndef _AcerDet_Analyse_Test_Histograms_
#define _AcerDet_Analyse_Test_Histograms_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../test/Test.h"
using namespace AcerDet::test;

#include "../core/IHistogramManager.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace analyse {

		//! Test Histograms creation algorithm class
		/*!
		 * An algorithm for a creation of test histograms.
		 */
		class Test_Histograms {
		private:
			
			Int32_t KEYHID;    /*!< see Configuration */
			Bool_t TEST;     /*!< see Configuration */

			Int32_t IEVENT;    /*!< Number of events computed by Test Histograms algorithm. */

			IHistogramManager *histoManager; /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;          /*!< Indicates whether histograms are already registered and might be used now. */
			
		public:
			
			/** Constructor.
			 * Creates new instance of Test Histograms creation algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histMgr a global instance of IHistogramManager
			 */
			Test_Histograms( const Configuration& conf, IHistogramManager* histMgr );
			
			/**
			 * Destructor.
			 */
			~Test_Histograms();
			
			/**
			 * Print information about Test Histograms algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input and output record to fill test histograms.
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
			 * Print Test Histograms algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
