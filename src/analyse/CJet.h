#ifndef _AcerDet_Analyse_C_Jets_
#define _AcerDet_Analyse_C_Jets_

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
		class CJet {
		private:
			Real64_t ETJET; /*!< detailed description  */
			Real64_t RCONE; /*!< detailed description  */
			Real64_t PTCMIN; /*!< detailed description  */
			Real64_t ETCMAX; /*!< detailed description  */
			Real64_t RJC; /*!< detailed description  */

			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYBCL; /*!< detailed description  */

			Int32_t IEVENT; /*!< detailed description  */

			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */
			
		public:
			CJet( const Configuration&, IHistogramManager* );
			~CJet();

			//! Print information about CJet algorithm class to standard output
			void printInfo() const;
			
			//! Analyse input record and find c-jets
			/*!
			 * \param input <input data desc>
			 * \param output <output data desc>
			 */
			void analyseRecord(
				const io::InputRecord& input,
				io::OutputRecord& output );

			//! Print CJet algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
