#ifndef _AcerDet_Analyse_Mis_
#define _AcerDet_Analyse_Mis_

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
		class Mis {
		private:
			Real64_t PTMUMIN; /*!< detailed description  */
			Real64_t ETAMAX; /*!< detailed description  */
			Real64_t ETCELL; /*!< detailed description  */
			Real64_t CALOTH; /*!< detailed description  */
			
			Int32_t KEYHID; /*!< detailed description  */
			Bool_t KEYSME; /*!< detailed description  */
			Int32_t KFINVS; /*!< detailed description  */
			
			Int32_t IEVENT; /*!< detailed description  */
			
			IHistogramManager *histoManager; /*!< detailed description  */
			Bool_t histoRegistered; /*!< detailed description  */
		
		public:
			Mis( const Configuration&, IHistogramManager* );
			~Mis();
			
			//! Print information about Mis algorithm class to standard output
			void printInfo() const;
			
			void analyseRecord( const io::InputRecord&, io::OutputRecord& );
			
			//! Print Mis algorithm execution results to standard output
			void printResults() const;
		};

	}
}

#endif
