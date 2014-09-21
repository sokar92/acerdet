#ifndef _AcerDet_Analyse_Muon_
#define _AcerDet_Analyse_Muon_

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
		class Muon {
		private:
			Real64_t ETCLU; /*!< detailed description  */
			Real64_t RCONE; /*!< detailed description  */
			Real64_t PTMUMIN; /*!< detailed description  */
			Real64_t ETAMAX; /*!< detailed description  */
			Real64_t RISOLJ; /*!< detailed description  */
			Real64_t RDEP; /*!< detailed description  */
			Real64_t EDMAX; /*!< detailed description  */
			Real64_t PTMUMINT; /*!< detailed description  */

			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYSME; /*!< detailed description  */

			Int32_t IEVENT; /*!< detailed description  */

			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */
			
		public:
			Muon( const Configuration&, IHistogramManager* );
			~Muon();
			
			//! Print information about Muon algorithm class to standard output
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			//! Print Muon algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
