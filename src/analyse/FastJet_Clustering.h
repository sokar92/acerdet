#ifndef _AcerDet_Analyse_FastJet_Clustering_
#define _AcerDet_Analyse_FastJet_Clustering_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../core/ParticleDataProvider.h"
using namespace AcerDet::core;

#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"
using namespace AcerDet::io;

#include "./IAnalysePhase.h"

namespace AcerDet {
	namespace analyse {

		//! FastJet_Cell analyse algorithm class
		/*!
		 * An algorithm for grouping particles into cells.
		 */
		class FastJet_Clustering : public IAnalysePhase {
		private:
			Real64_t ETACEL; /*!< see Configuration */
			Real64_t PTMIN;  /*!< see Configuration */
			Real64_t ETTHR;  /*!< see Configuration */
			Real64_t CALOTH; /*!< see Configuration */
			Real64_t DBETA;  /*!< see Configuration */
			Real64_t DBPHI;  /*!< see Configuration */
			
			Int32_t KEYHID; /*!< Histogram unique id */
			Bool_t  KEYFLD; /*!< see Configuration */
			Int32_t KFINVS; /*!< see Configuration */
			
			Int32_t IEVENT; /*!< Number of events computed by Cell algorithm. */
			
			IHistogramManager *histoManager;          /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;                   /*!< Indicates whether histograms are already registered and might be used now. */

			const ParticleDataProvider& partProvider; /*!< External provider for particle properties such as charge or name. */
			
		public:
		
			/**
			 * Constructor.
			 * Creates new instance of Cell analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histMgr a global instance of IHistogramManager
			 * \param partData a global instance of IParticleDataProvider
			 */
			FastJet_Clustering( 
				const Configuration& conf,
				IHistogramManager* histMgr,
				const ParticleDataProvider& partData );

			/**
			 * Destructor.
			 */
			~FastJet_Clustering();
			
			/**
			 * Print information about Cell algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and group cells.
			/*!
			 * \param input record.
			 * \param output record.
			 * \param weight event weight.
			 */
			void analyseRecord(
				const InputRecord& input,
				OutputRecord& output,
				Real64_t weight ); 
			
			/**
			 * Print Cell algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
