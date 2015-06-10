#ifndef _AcerDet_Analyse_Muon_
#define _AcerDet_Analyse_Muon_

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

		//! Muon analyse algorithm class
		/*!
		 * An algorithm for Muon reconstruction.
		 */
		class Muon {
		private:
			Real64_t ETCLU;    /*!< see Configuration */
			Real64_t RCONE;    /*!< see Configuration */
			Real64_t PTMUMIN;  /*!< see Configuration */
			Real64_t ETAMAX;   /*!< see Configuration */
			Real64_t RISOLJ;   /*!< see Configuration */
			Real64_t RDEP;     /*!< see Configuration */
			Real64_t EDMAX;    /*!< see Configuration */
			Real64_t PTMUMINT; /*!< see Configuration */

			Int32_t KEYHID;    /*!< see Configuration */
			Bool_t KEYSME;     /*!< see Configuration */

			Int32_t IEVENT;    /*!< Number of events computed by Muon algorithm. */

			IHistogramManager *histoManager; /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;          /*!< Indicates whether histograms are already registered and might be used now. */
			
		public:
			
			/** Constructor.
			 * Creates new instance of Muon analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histMgr a global instance of IHistogramManager
			 */
			Muon( const Configuration& conf, IHistogramManager* histMgr );
			
			/**
			 * Destructor.
			 */
			~Muon();
			
			/**
			 * Print information about Muon algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and find muons.
			/*!
			 * \param input record.
			 * \param output record.
			 * \param weight event weigth.
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output,
				Real64_t weigth );
			
			/**
			 * Print Muon algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
