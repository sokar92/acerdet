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

namespace AcerDet {
	namespace analyse {

		//! Title
		/*!
		 * Detailed description
		 */
		class Tau {
		private:
			Real64_t ETJET; /*!< detailed description  */
			Real64_t PTTAU; /*!< detailed description  */
			Real64_t ETATAU; /*!< detailed description  */
			Real64_t RJTAU; /*!< detailed description  */
			Real64_t PTFRAC; /*!< detailed description  */

			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYTAU; /*!< detailed description  */

			Int32_t IEVENT; /*!< detailed description  */

			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */
		public:
			Tau( const Configuration&, IHistogramManager* );
			~Tau();

			//! Print information about Tau algorithm class to standard output
			void printInfo() const;
			
			//! Analyse input record and find taus
			/*!
			 * \param input <input data desc>
			 * \param output <output data desc>
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );

			//! Print Tau algorithm execution results to standard output
			void printResults() const;
		};
		
	}
}

#endif
