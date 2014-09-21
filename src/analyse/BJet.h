#ifndef _AcerDet_Analyse_B_Jets_
#define _AcerDet_Analyse_B_Jets_

#pragma once

#include "../conf/Configuration.h"
using namespace AcerDet::conf;

#include "../core/IHistogramManager.h"
#include "../io/InputRecord.h"
#include "../io/OutputRecord.h"

namespace AcerDet {
	namespace analyse {

		//! Title
		/*!
		 * Detailed description
		 */		
		class BJet {
		private:
			Real64_t ETJET; /*!< detailed description  */
			Real64_t RCONE; /*!< detailed description  */
			Real64_t PTBMIN; /*!< detailed description  */
			Real64_t ETBMAX; /*!< detailed description  */
			Real64_t RJB; /*!< detailed description  */

			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYBCL; /*!< detailed description  */

			Int32_t IEVENT; /*!< detailed description  */
			
			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */

		public:
			BJet( const Configuration&, IHistogramManager* );
			~BJet();
			
			//! Print information about BJet algorithm class to standard output
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			//! Print BJet algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
