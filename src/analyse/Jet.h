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

		//! Title
		/*!
		 * Detailed description
		 */
		class Jet {
		private:
			Real64_t ETJET; /*!< detailed description  */
			Real64_t ETAJET; /*!< detailed description  */
			Real64_t RCONE; /*!< detailed description  */
			Real64_t PTMIN; /*!< detailed description  */
			Real64_t CALOTH; /*!< detailed description  */

			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYSME; /*!< detailed description  */
			Bool_t KEYFLD; /*!< detailed description  */
			Int32_t KFINVS; /*!< detailed description  */

			Int32_t IEVENT; /*!< detailed description  */
			
			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */
			const ParticleDataProvider& partProvider; /*!< detailed description  */
			
		public:
			Jet(
				const Configuration&,
				IHistogramManager*,
				const ParticleDataProvider& );

			~Jet();
			
			//! Print information about Jet algorithm class to standard output
			void printInfo() const;
			
			//! Analyse input record and find jets
			/*!
			 * \param input <input data desc>
			 * \param output <output data desc>
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );
			
			//! Print Jet algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
