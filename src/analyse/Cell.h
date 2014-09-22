#ifndef _AcerDet_Analyse_Cell_
#define _AcerDet_Analyse_Cell_

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

		//! Cell algorithm class
		/*!
		 * An algorithm for grouping particles into cells.
		 */
		class Cell {
		private:
			Real64_t ETACEL; /*!< detailed description  */
			Real64_t PTMIN; /*!< detailed description  */
			Real64_t ETTHR; /*!< detailed description  */
			Real64_t CALOTH; /*!< detailed description  */
			Real64_t DBETA; /*!< detailed description  */
			Real64_t DBPHI; /*!< detailed description  */
			
			Int32_t KEYHID; /*!< Histogram unique id */
			Bool_t  KEYFLD; /*!< Should cell mathcing algorithm use folding? */
			Int32_t KFINVS; /*!< Should cell matching algorithm care about invisible particles?  */
			
			Int32_t IEVENT; /*!< Number of events computed by this class */
			
			IHistogramManager *histoManager; /*!< Histogram manager instance */
			Bool_t histoRegistered; /*!< Histogram manager initialization state true if initialized false otherwise */

			const ParticleDataProvider& partProvider; /*!< External provider for particle properties such as charge or name */
			
		public:
			/** A constructor
			 * \param configuration global AcerDET configuration
			 * \param histogramManager histogram manager instance
			 * \param particleData provider of external particle data (such as charge and name) 
			 */
			Cell( 
				const Configuration& configuration,
				IHistogramManager* histogramManager,
				const ParticleDataProvider& particleData );

			/** A destructor
			 */
			~Cell();
			
			//! Print information about Cell algorithm class to standard output
			void printInfo() const;
			
			//! Analyse input record group cells
			/*!
			 * \param input <input data desc>
			 * \param output <output data desc>
			 */
			void analyseRecord(
				const InputRecord& input,
				OutputRecord& output ); 
			
			//! Print Cell algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
