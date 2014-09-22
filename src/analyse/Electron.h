#ifndef _AcerDet_Analyse_Electron_
#define _AcerDet_Analyse_Electron_

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
		class Electron {
		private:
			Real64_t ETCLU; /*!< detailed description  */
			Real64_t RCONE; /*!< detailed description  */
			Real64_t PTLMIN; /*!< detailed description  */
			Real64_t ETAMAX; /*!< detailed description  */
			Real64_t RJE; /*!< detailed description  */
			Real64_t RISOLJ; /*!< detailed description  */
			Real64_t RDEP; /*!< detailed description  */
			Real64_t EDMAX; /*!< detailed description  */

			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYSME; /*!< detailed description  */

			Int32_t IEVENT; /*!< detailed description  */
			
			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */

		public:
			Electron( const Configuration&, IHistogramManager* );
			~Electron();
			
			//! Print information about Electron algorithm class to standard output
			void printInfo() const;
			
			//! Analyse input record and find electrons
			/*!
			 * \param input <input data desc>
			 * \param output <output data desc>
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );
			
			//! Print Electron algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
