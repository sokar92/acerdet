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

		//! Title
		/*!
		 * Detailed description
		 */
		class Cell {
		private:
			Real64_t ETACEL; /*!< detailed description  */
			Real64_t PTMIN; /*!< detailed description  */
			Real64_t ETTHR; /*!< detailed description  */
			Real64_t CALOTH; /*!< detailed description  */
			Real64_t DBETA; /*!< detailed description  */
			Real64_t DBPHI; /*!< detailed description  */
			
			Int32_t KEYHID; /*!< detailed description  */
			Bool_t  KEYFLD; /*!< detailed description  */
			Int32_t KFINVS; /*!< detailed description  */
			
			Int32_t IEVENT; /*!< detailed description  */
			
			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */

			const ParticleDataProvider& partProvider; /*!< detailed description  */
			
		public:
			Cell( 
				const Configuration&,
				IHistogramManager*,
				const ParticleDataProvider& );

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
