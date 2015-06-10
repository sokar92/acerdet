#ifndef _AcerDet_Analyse_Cluster_
#define _AcerDet_Analyse_Cluster_

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

		//! Cluster analyse algorithm class
		/*!
		 * An algorithm for grouping cells into clusters.
		 */
		class Cluster {
		private:
			Real64_t ETCLU;  /*!< see Configuration */
			Real64_t RCONE;  /*!< see Configuration */
			Real64_t ETACLU; /*!< see Configuration */
			Real64_t ETINI;  /*!< see Configuration */
			
			Real64_t PTMIN;  /*!< see Configuration */
			Real64_t CALOTH; /*!< see Configuration */
			
			Int32_t KEYHID;  /*!< see Configuration */
			Bool_t KEYFLD;   /*!< see Configuration */
			Int32_t KFINVS;  /*!< see Configuration */
			
			Int32_t IEVENT;  /*!< Number of events computed by Cluster algorithm. */
			
			IHistogramManager *histoManager;          /*!< A pointer to global histogram manager. */
			Bool_t histoRegistered;                   /*!< Indicates whether histograms are already registered and might be used now. */
			const ParticleDataProvider& partProvider; /*!< External provider for particle properties such as charge or name. */

		public:
		
			/**
			 * Constructor.
			 * Creates new instance of Cell analyse algorithm.
			 * \param conf AcerDET initial Configuration
			 * \param histoMgr a global instance of IHistogramManager
			 * \param partData a global instance of IParticleDataProvider
			 */
			Cluster(
				const Configuration& conf,
				IHistogramManager* histoMgr,
				const ParticleDataProvider& partData );

			/**
			 * Destructor.
			 */
			~Cluster();
			
			/**
			 * Print information about Cluster algorithm class to standard output.
			 */
			void printInfo() const;
			
			//! Analyse input record and group cells into clusters.
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
			 * Print Cluster algorithm execution results to standard output.
			 */
			void printResults() const;
		};

	}
}

#endif
