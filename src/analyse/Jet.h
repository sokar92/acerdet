#ifndef _AcerDet_Analyse_Jet_
#define _AcerDet_Analyse_Jet_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../core/ParticleDataProvider.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

namespace AcerDet {
	namespace analyse {

		//! Jet analyse algorithm class
		/*!
		 * An algorithm for Jet reconstruction.
		 */
		class Jet {
		private:
			Real64_t ETJET;   /*!< see Configuration */
			Real64_t ETAJET;  /*!< see Configuration */
			Real64_t RCONE;   /*!< see Configuration */
			Real64_t PTMIN;   /*!< see Configuration */
			Real64_t CALOTH;  /*!< see Configuration */

			Int32_t KEYHID;   /*!< see Configuration */
			Bool_t KEYSME;    /*!< see Configuration */
			Bool_t KEYFLD;    /*!< see Configuration */
			Int32_t KFINVS;   /*!< see Configuration */

			Int32_t IEVENT;   /*!< Number of events computed by Jet algorithm. */
			
			IHistogramManager *histoManager;          /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;                   /*!< Indicates whether histograms are already registered and might be used now. */
			const ParticleDataProvider& partProvider; /*!< External provider for particle properties such as charge or name. */
			
		public:
			
			/** Constructor.
			 * Creates new instance of Jet analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histMgr a global instance of IHistogramManager
			 * \param partData a global instance of IParticleDataProvider
			 */
			Jet(
				const Configuration& conf,
				IHistogramManager* histMgr,
				const ParticleDataProvider& partData );

			/**
			 * Destructor.
			 */
			~Jet();
			
			/**
			 * Print information about Jet algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and find jets.
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
			 * Print Jet algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
